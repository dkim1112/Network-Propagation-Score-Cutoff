################################################################################
# Network Propagation Score Cutoff
# Method: Degree-preserving permutation based empirical p-value
#
# 사용법 (김동은 PC기준):
#   & "C:\Program Files\R\R-4.6.0\bin\Rscript.exe" "C:\Users\de7da\OneDrive\Desktop\NPCutOff\np_cutoff_permutation.R"
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

user_lib <- file.path(Sys.getenv("APPDATA"), "R", "win-library", "4.6")
if (dir.exists(user_lib)) {
  .libPaths(unique(c(user_lib, .libPaths())))
}

library(igraph)

# == 파라미터 ==================================================================
DATA_DIR           <- "result_network_propagation"
NETWORK_FILE       <- "tables_expansion/Combined_STRINGv11_OTAR281119_FILTER.rds"
OUTPUT_FILE        <- "result_np_cutoff/np_cutoff_results.csv"
N_PERM             <- 1000
LOW_SEED_THRESHOLD <- 2
RANDOM_SEED        <- 42
FDR_ALPHA          <- 0.05
# seed 수를 다양하게 가져간 10개 질환 (seed 1~25)
# 전체 실행 시 TEST_SPECIFIC <- NULL 로 변경
TEST_SPECIFIC <- c(
  "nodes.finngen_R12_ALCOPANCCHRON.rds",          # seed 1
  "nodes.finngen_R12_ABDOM_HERNIA.rds",            # seed 2
  "nodes.finngen_R12_L12_ATOPIC.rds",              # seed 3
  "nodes.finngen_R12_T1D.rds",                     # seed 5
  "nodes.finngen_R12_AUTOIMMUNE_NONTHYROID.rds",   # seed 7
  "nodes.finngen_R12_T2D_WIDE.rds",                # seed 9
  "nodes.finngen_R12_I9_CHD.rds",                  # seed 11
  "nodes.finngen_R12_AUTOIMMUNE.rds",              # seed 13
  "nodes.finngen_R12_K11_IBD_STRICT.rds",          # seed 16
  "nodes.finngen_R12_I9_HYPTENS.rds"               # seed 25
)
# =============================================================================

set.seed(RANDOM_SEED)

dir.create(dirname(OUTPUT_FILE), showWarnings = FALSE, recursive = TRUE)

# == Nearest-neighbor degree matching 함수 =====================================
# 박사님 코드(AUC_calculation.R)의 make_nearest_match_sampler 기반
# 각 seed에 대해 degree 차이 0부터 점진적으로 허용 범위를 늘려가며
# 가장 가까운 degree의 random 유전자를 선택.
# -> quantile bin 방식 대비 더 정밀하고 임의적 bin 경계 없음.
make_nearest_match_sampler <- function(
    deg_all,                         # named numeric: 모든 노드의 degree
    tier1_genes,                     # seed 유전자 (degree matching 기준)
    exclude_genes = tier1_genes,     # 샘플링 pool에서 제외할 유전자 (기본: seed 자신)
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
        cand <- setdiff(cand, picked)   # 중복 방지
        if (length(cand) > 0) { got <- sample(cand, 1); break }
      }
      if (is.na(got)) stop(sprintf("근접 매칭 실패: %s (degree=%s)", g, d))
      picked <- c(picked, got)
    }
    picked
  }
  sampler
}

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
all_nodes <- V(net_full)$name

# 전체 네트워크 degree 사전 계산 (sampler에 전달)
net_deg_vec <- igraph::degree(net_full)   # named numeric vector

# base personalization vector
base_pv <- setNames(rep(0.0, length(all_nodes)), all_nodes)

# == 3. 단일 질환 처리 함수 ====================================================
process_disease <- function(filepath) {

  node <- as.data.frame(readRDS(filepath), stringsAsFactors = FALSE)
  node$padj      <- as.numeric(node$padj)
  node$page.rank <- as.numeric(node$page.rank)

  trait     <- unique(node$Trait)[1]
  seed_df   <- node[!is.na(node$padj) & node$padj != 0, ]
  target_df <- node[!is.na(node$padj) & node$padj == 0, ] # seed 제외, padj==0인 유전자만 분석 대상
  n_seeds   <- nrow(seed_df)

  if (n_seeds == 0) return(NULL)

  obs_scores   <- target_df$page.rank
  target_nodes <- target_df$ENSG

  # Nearest-neighbor degree matched sampler 생성
  # tier1_genes = 실제 seed, exclude_genes = seed (샘플링 pool에서 제외)
  sampler <- make_nearest_match_sampler(
    deg_all      = net_deg_vec,
    tier1_genes  = seed_df$ENSG,
    exclude_genes = seed_df$ENSG
  )

  # == Permutation =============================================================
  null_matrix <- matrix(NA_real_, nrow = length(target_nodes), ncol = N_PERM)

  for (perm_idx in seq_len(N_PERM)) {

    # nearest-neighbor degree matched random seed 선택
    rand_seeds <- sampler()

    # personalization vector
    pv    <- base_pv
    valid <- rand_seeds %in% names(pv)
    pv[rand_seeds[valid]] <- seed_df$padj[valid]

    # page_rank 재실행 = RWR 재실행
    pr_vec <- page_rank(net_full,
                         personalized = pv,
                         weights      = E(net_full)$weight)$vector

    # 각 유전자별 null score을 null_matrix에 저장
    null_matrix[, perm_idx] <- pr_vec[target_nodes]
  }

  # == Empirical p-value =======================================================
  count_ge <- rowSums(null_matrix >= obs_scores, na.rm = TRUE)
  emp_pval <- (count_ge + 1) / (N_PERM + 1)

  target_df$emp_pval       <- emp_pval
  target_df$significant    <- emp_pval < FDR_ALPHA
  target_df$n_seeds        <- n_seeds
  target_df$low_confidence <- (n_seeds <= LOW_SEED_THRESHOLD)

  target_df[, c("ENSG", "gene", "Trait", "page.rank", "degree",
                 "emp_pval", "significant", "n_seeds", "low_confidence")]
}

# == 4. 전체 질환 처리 =========================================================
files <- sort(list.files(DATA_DIR, pattern = "\\.rds$", full.names = TRUE))
if (!is.null(TEST_SPECIFIC)) files <- files[basename(files) %in% TEST_SPECIFIC]
n_files <- length(files)

cat(sprintf("Processing %d diseases...\n", n_files))

all_results <- vector("list", n_files)
skipped     <- character()

for (i in seq_along(files)) {
  fname <- basename(files[i])
  cat(sprintf("[%d/%d] %s ", i, n_files, gsub("nodes\\.finngen_R12_|\\.rds", "", fname)))

  tryCatch({
    res <- process_disease(files[i])
    if (!is.null(res)) {
      all_results[[i]] <- res
      cat(sprintf("✓ (seeds=%d)\n", res$n_seeds[1]))
    } else {
      skipped <- c(skipped, fname)
      cat("✗ (no seeds)\n")
    }
  }, error = function(e) {
    skipped <<- c(skipped, fname)
    cat(sprintf("✗ (error)\n"))
  })
}

# == 5. 결과 저장 ==============================================================
all_results <- Filter(Negate(is.null), all_results)
if (length(all_results) == 0) stop("No results generated.")

final <- do.call(rbind, all_results)
write.csv(final, OUTPUT_FILE, row.names = FALSE)

cat(sprintf("\nCompleted: %d diseases processed, %d successful\n", n_files, length(all_results)))
cat(sprintf("Results saved to: %s\n", OUTPUT_FILE))

# == 요약 ======================================================================
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