################################################################################
# Network Propagation Score Cutoff — 병렬화 버전 (서버용)
# Method: Degree-preserving permutation based empirical p-value
#
# 사용법 (서버):
#   Rscript np_cutoff_permutation_parallel.R
#
# 코어 수 자동 감지 → 전체 코어의 75% 사용 (안전 마진)
# 수동 지정: N_CORES <- 8 처럼 숫자로 설정
#
# 절차:
# 1. padj == 0인 유전자만 분석 대상 (seed 제외)
# 2. seed와 동일한 수 + nearest-neighbor degree matched random seeds로 N_PERM회 permutation
#    → 매번 page_rank 재실행 → 유전자별 null distribution 생성
#    → degree matching 방식: 각 seed에 대해 degree 차이 0부터 점진적으로 허용
#      (quantile bin 방식 → nearest-neighbor 방식으로 변경, 박사님 코드 기반)
# 3. empirical p-value: p_i = (#{null >= obs} + 1) / (N_PERM + 1)
# 4. seed 수 <= LOW_SEED_THRESHOLD 질환은 low_confidence flag
# ※ BH FDR 보정 미적용: N_PERM=1000 기준 최솟값 p=0.001 x 18408 = 18.4로
#    수학적으로 통과 불가. empirical p-value 자체를 cutoff로 사용. -> 우선은 qval없음.
#
################################################################################

library(igraph)
library(parallel)

# == 파라미터 ==================================================================
DATA_DIR           <- "result_network_propagation"
NETWORK_FILE       <- "tables_expansion/Combined_STRINGv11_OTAR281119_FILTER.rds"
OUTPUT_FILE        <- "result_np_cutoff/np_cutoff_results.csv"
N_PERM             <- 1000
LOW_SEED_THRESHOLD <- 2
RANDOM_SEED        <- 42
FDR_ALPHA          <- 0.05

# 코어 수: NULL이면 자동 감지 (전체의 75%), 숫자로 직접 지정 가능
# 예: N_CORES <- 8
N_CORES <- NULL

# 전체 실행 시 TEST_SPECIFIC <- NULL 로 변경
TEST_SPECIFIC <- c(
  # old datasets
  # "nodes.finngen_R12_ALCOPANCCHRON.rds",
  # "nodes.finngen_R12_ABDOM_HERNIA.rds",
  # "nodes.finngen_R12_L12_ATOPIC.rds",
  # "nodes.finngen_R12_T1D.rds",
  # "nodes.finngen_R12_AUTOIMMUNE_NONTHYROID.rds",
  # "nodes.finngen_R12_T2D_WIDE.rds",
  # "nodes.finngen_R12_I9_CHD.rds",
  # "nodes.finngen_R12_AUTOIMMUNE.rds",
  # "nodes.finngen_R12_K11_IBD_STRICT.rds",
  # "nodes.finngen_R12_I9_HYPTENS.rds"

  # seed >= 5개인 질환들
  "nodes.finngen_R12_I9_HYPTENS.rds",                    # seed 25
  "nodes.finngen_R12_K11_IBD_STRICT.rds",                # seed 16
  "nodes.finngen_R12_AUTOIMMUNE.rds",                    # seed 13
  "nodes.finngen_R12_I9_CHD.rds",                        # seed 11
  "nodes.finngen_R12_T2D_WIDE.rds",                      # seed 9
  "nodes.finngen_R12_C3_BASAL_CELL_CARCINOMA_EXALLC.rds",# seed 12
  "nodes.finngen_R12_C3_PROSTATE_EXALLC.rds",            # seed 8
  "nodes.finngen_R12_C3_SKIN_EXALLC.rds",                # seed 8
  "nodes.finngen_R12_CARDIAC_ARRHYTM.rds",               # seed 9
  "nodes.finngen_R12_T1D.rds",                           # seed 5
  "nodes.finngen_R12_ASTHMMA_ACUTE_RESPIRATORY_INFECTIONS.rds", # seed 7
  "nodes.finngen_R12_AD_EO_EXMORE.rds",                  # seed 5
  "nodes.finngen_R12_ALLERG_ASTHMA.rds",                 # seed 6
  "nodes.finngen_R12_ASTHMA_CHILD_EXMORE.rds"          # seed 5
)
# =============================================================================

# 코어 수 결정
if (is.null(N_CORES)) {
  detected <- detectCores(logical = FALSE)   # 물리 코어
  N_CORES  <- max(1L, floor(detected * 0.75))
}
cat(sprintf("사용 코어 수: %d / %d\n", N_CORES, detectCores(logical = FALSE)))

dir.create(dirname(OUTPUT_FILE), showWarnings = FALSE, recursive = TRUE)

# == Nearest-neighbor degree matching 함수 =====================================
make_nearest_match_sampler <- function(
    deg_all,
    tier1_genes,
    exclude_genes = tier1_genes,
    expand_seq = c(0, 1, 2, 3, 5, 8, 10, 12, 15,
                   20, 30, 40, 50, 60, 70, 80, 90, 100, 200, Inf)
) {
  stopifnot(is.numeric(deg_all), !is.null(names(deg_all)))
  tier1 <- intersect(tier1_genes, names(deg_all))
  pool  <- setdiff(names(deg_all), exclude_genes)

  sampler <- function() {
    picked <- character(0)
    for (g in tier1) {
      d   <- deg_all[g]
      got <- NA
      for (tol in expand_seq) {
        cand <- pool[abs(deg_all[pool] - d) <= tol]
        cand <- setdiff(cand, picked)
        if (length(cand) > 0) { got <- sample(cand, 1); break }
      }
      if (is.na(got)) stop(sprintf("근접 매칭 실패: %s (degree=%s)", g, d))
      picked <- c(picked, got)
    }
    picked
  }
  sampler
}

# == 네트워크 로드 (메인 프로세스에서 한 번만) =================================
cat("Loading PPI network...\n")
string <- as.data.frame(readRDS(NETWORK_FILE), stringsAsFactors = FALSE)
string$combined_score <- as.numeric(string$combined_score)

net_full <- graph_from_data_frame(
  d        = string[, c("ENSG_A", "ENSG_B")],
  directed = FALSE
)
E(net_full)$weight <- string$combined_score
net_full <- igraph::simplify(net_full,
                              remove.loops    = TRUE,
                              remove.multiple = TRUE,
                              edge.attr.comb  = c(weight = "max", "ignore"))
all_nodes   <- V(net_full)$name
net_deg_vec <- igraph::degree(net_full)
base_pv     <- setNames(rep(0.0, length(all_nodes)), all_nodes)
cat(sprintf("Network: %d nodes, %d edges\n\n", length(all_nodes), ecount(net_full)))

# == 단일 질환 처리 함수 =======================================================
# mclapply가 각 worker에서 이 함수를 실행
# net_full, net_deg_vec, base_pv, all_nodes는 fork로 공유됨 (Linux)
process_disease <- function(filepath) {

  node <- as.data.frame(readRDS(filepath), stringsAsFactors = FALSE)
  node$padj      <- as.numeric(node$padj)
  node$page.rank <- as.numeric(node$page.rank)

  trait     <- unique(node$Trait)[1]
  seed_df   <- node[!is.na(node$padj) & node$padj != 0, ]
  target_df <- node[!is.na(node$padj) & node$padj == 0, ]
  n_seeds   <- nrow(seed_df)

  if (n_seeds == 0) return(NULL)

  obs_scores   <- target_df$page.rank
  target_nodes <- target_df$ENSG

  sampler <- make_nearest_match_sampler(
    deg_all       = net_deg_vec,
    tier1_genes   = seed_df$ENSG,
    exclude_genes = seed_df$ENSG
  )

  null_matrix <- matrix(NA_real_, nrow = length(target_nodes), ncol = N_PERM)

  for (perm_idx in seq_len(N_PERM)) {
    rand_seeds <- sampler()
    pv    <- base_pv
    valid <- rand_seeds %in% names(pv)
    pv[rand_seeds[valid]] <- seed_df$padj[valid]
    pr_vec <- page_rank(net_full,
                         personalized = pv,
                         weights      = E(net_full)$weight)$vector
    null_matrix[, perm_idx] <- pr_vec[target_nodes]
  }

  count_ge <- rowSums(null_matrix >= obs_scores, na.rm = TRUE)
  emp_pval <- (count_ge + 1) / (N_PERM + 1)

  target_df$emp_pval       <- emp_pval
  target_df$significant    <- emp_pval < FDR_ALPHA
  target_df$n_seeds        <- n_seeds
  target_df$low_confidence <- (n_seeds <= LOW_SEED_THRESHOLD)

  target_df[, c("ENSG", "gene", "Trait", "page.rank", "degree",
                 "emp_pval", "significant", "n_seeds", "low_confidence")]
}

# == 병렬 실행 =================================================================
files <- sort(list.files(DATA_DIR, pattern = "\\.rds$", full.names = TRUE))
if (!is.null(TEST_SPECIFIC)) files <- files[basename(files) %in% TEST_SPECIFIC]
n_files <- length(files)

cat(sprintf("Processing %d diseases across %d cores...\n", n_files, N_CORES))

# 시작할 질환들 표시
diseases <- gsub("nodes\\.finngen_R12_|\\.rds", "", basename(files))
cat(sprintf("Diseases: %s\n\n", paste(diseases, collapse=", ")))

t0 <- proc.time()["elapsed"]

# 병렬 처리용 wrapper 함수 (진행상황 출력 포함)
process_with_progress <- function(f) {
  fname <- gsub("nodes\\.finngen_R12_|\\.rds", "", basename(f))
  cat(sprintf("Starting: %s\n", fname))

  tryCatch({
    result <- process_disease(f)
    if (!is.null(result)) {
      cat(sprintf("✓ Completed: %s (seeds=%d)\n", fname, result$n_seeds[1]))
      return(result)
    } else {
      cat(sprintf("✗ Skipped: %s (no seeds)\n", fname))
      return(NULL)
    }
  }, error = function(e) {
    cat(sprintf("✗ Error: %s (%s)\n", fname, conditionMessage(e)))
    return(NULL)
  })
}

# mclapply: Linux/Mac에서 fork 기반 병렬화
set.seed(RANDOM_SEED)
all_results <- mclapply(
  files,
  process_with_progress,
  mc.cores    = N_CORES,
  mc.set.seed = TRUE
)

elapsed <- proc.time()["elapsed"] - t0
cat(sprintf("\nParallel processing completed in %.0f seconds (%.1f minutes)\n", elapsed, elapsed/60))

# == 결과 저장 =================================================================
all_results <- Filter(Negate(is.null), all_results)
if (length(all_results) == 0) stop("No results generated.")

final <- do.call(rbind, all_results)
write.csv(final, OUTPUT_FILE, row.names = FALSE)

cat(sprintf("Results saved: %s\n", OUTPUT_FILE))
cat(sprintf("Successfully processed: %d/%d diseases\n", length(all_results), n_files))

cat("\n", strrep("=", 55), "\n", sep = "")
cat(sprintf("완료. 총 %d개 질환 처리\n",    length(unique(final$Trait))))
cat(sprintf("전체 유전자-질환 쌍: %s\n",     format(nrow(final), big.mark = ",")))
cat(sprintf("emp_pval < 0.05 유전자-질환 쌍: %s\n",
            format(sum(final$emp_pval < 0.05), big.mark = ",")))
cat(sprintf("Low confidence 질환 (seed<=%d개): %d개 질환\n",
            LOW_SEED_THRESHOLD,
            length(unique(final$Trait[final$low_confidence]))))
if (length(skipped) > 0)
  cat("Skipped:", paste(skipped, collapse = ", "), "\n")
cat(sprintf("결과 저장: %s\n", OUTPUT_FILE))