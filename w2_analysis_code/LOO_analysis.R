################################################################################
# Leave-One-Out (LOO) Analysis
# 목적: permutation 기반 empirical p-value 방법론 검증
#
# 절차:
# 1. 각 질환의 seed 중 하나를 제외 (leave-one-out)
# 2. 나머지 seed로 propagation + permutation 수행
# 3. 제외한 seed의 emp_pval 확인
#    → 유의(p < 0.05): 나머지 seed만으로도 해당 단백질 recovery → 방법론 유효
#    → 비유의:  permutation 수 부족하거나 해당 seed가 독립적 위치
#
# 조건: seed 수 >= 2인 질환만 (seed가 1개면 LOO 후 남는 seed가 없음)
# + Assuming we have access to server - 병렬처리 할 것임.
################################################################################

library(igraph)
library(parallel)
if (!requireNamespace("pbmcapply", quietly = TRUE)) install.packages("pbmcapply")
library(pbmcapply) # 실시간 프로그레스 바를 위한 패키지 추가

# == 파라미터 ==================================================================
# Set absolute paths based on project structure
base_dir <- "/Users/kde/Documents/Network-Propagation-Score-Cutoff"

DATA_DIR           <- file.path(base_dir, "result_network_propagation")
NETWORK_FILE       <- file.path(base_dir, "tables_expansion/Combined_STRINGv11_OTAR281119_FILTER.rds")
OUTPUT_FILE        <- file.path(base_dir, "result_np_cutoff/loo_results.csv")
N_PERM             <- 10000
FDR_ALPHA          <- 0.05
N_CORES            <- 48   # 서버 코어 수 직접 지정

# 테스트할 질환 (seed >= 2인 것 위주)
# 전체 실행: TEST_SPECIFIC <- NULL
TEST_SPECIFIC <- c(
  # "nodes.finngen_R12_ABDOM_HERNIA.rds",            # seed 2
  # "nodes.finngen_R12_L12_ATOPIC.rds",              # seed 3
  # "nodes.finngen_R12_T1D.rds",                     # seed 5
  # "nodes.finngen_R12_AUTOIMMUNE_NONTHYROID.rds"   # seed 7
  # "nodes.finngen_R12_T2D_WIDE.rds",                # seed 9
  "nodes.finngen_R12_I9_CHD.rds",                  # seed 11
  "nodes.finngen_R12_AUTOIMMUNE.rds",              # seed 13
  "nodes.finngen_R12_K11_IBD_STRICT.rds"          # seed 16
  # "nodes.finngen_R12_I9_HYPTENS.rds"               # seed 25
)
# =============================================================================

if (is.null(N_CORES)) {
  detected <- detectCores(logical = FALSE)
  N_CORES  <- max(1L, floor(detected * 0.75))
}
cat(sprintf("사용 코어 수: %d\n", N_CORES))

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

# == 네트워크 로드 =============================================================
cat("Loading PPI network...\n")
string <- as.data.frame(readRDS(NETWORK_FILE), stringsAsFactors = FALSE)
string$combined_score <- as.numeric(string$combined_score)
net_full <- graph_from_data_frame(d = string[, c("ENSG_A","ENSG_B")],
                                   directed = FALSE)
E(net_full)$weight <- string$combined_score
net_full <- igraph::simplify(net_full, remove.loops = TRUE,
                              remove.multiple = TRUE,
                              edge.attr.comb = c(weight="max","ignore"))
all_nodes   <- V(net_full)$name
net_deg_vec <- igraph::degree(net_full)
base_pv     <- setNames(rep(0.0, length(all_nodes)), all_nodes)
cat(sprintf("Network: %d nodes, %d edges\n\n", length(all_nodes), ecount(net_full)))

# == LOO 단일 실행 함수 ========================================================
# filepath: 질환 rds 파일
# left_out_idx: seed_df에서 제외할 seed의 인덱스
run_loo_single <- function(filepath, left_out_idx) {

  node <- as.data.frame(readRDS(filepath), stringsAsFactors = FALSE)
  node$padj      <- as.numeric(node$padj)
  node$page.rank <- as.numeric(node$page.rank)

  trait       <- unique(node$Trait)[1]
  all_seed_df <- node[!is.na(node$padj) & node$padj != 0, ]
  n_seeds     <- nrow(all_seed_df)

  if (n_seeds < 2) return(NULL)  # LOO 불가

  # 제외할 seed
  left_out_gene  <- all_seed_df$ENSG[left_out_idx]
  left_out_name  <- all_seed_df$gene[left_out_idx]
  left_out_padj  <- all_seed_df$padj[left_out_idx]

  # 나머지 seed
  remaining_seed_df <- all_seed_df[-left_out_idx, ]

  # 진행 상황 출력 (더 간결하게)
  if (left_out_idx == 1) {
    cat(sprintf("Processing %s (%d seeds)...\n",
                gsub("nodes.finngen_R12_", "", gsub(".rds", "", basename(filepath))),
                n_seeds))
  }

  # 분석 대상: seed 제외 유전자 전체 + left_out 포함
  target_df    <- node[!is.na(node$padj) & node$padj == 0, ]
  left_out_row <- all_seed_df[left_out_idx, ]
  left_out_row$padj <- 0
  target_df    <- rbind(target_df, left_out_row)
  target_nodes <- target_df$ENSG

  # 나머지 seed로 sampler 구성
  sampler <- make_nearest_match_sampler(
    deg_all       = net_deg_vec,
    tier1_genes   = remaining_seed_df$ENSG,
    exclude_genes = remaining_seed_df$ENSG
  )

  # obs_scores: n-1 seed로 page_rank 재실행
  # (원본 파일 page.rank는 left_out이 seed에 포함된 상태 → 사용 불가)
  pv_obs <- base_pv
  valid_obs <- remaining_seed_df$ENSG %in% names(pv_obs)
  pv_obs[remaining_seed_df$ENSG[valid_obs]] <- remaining_seed_df$padj[valid_obs]
  pr_obs     <- page_rank(net_full, personalized = pv_obs,
                           weights = E(net_full)$weight)$vector
  obs_scores <- pr_obs[target_nodes]

  # Permutation
  null_matrix <- matrix(NA_real_, nrow = length(target_nodes), ncol = N_PERM)

  for (perm_idx in seq_len(N_PERM)) {
    rand_seeds <- sampler()
    pv    <- base_pv
    valid <- rand_seeds %in% names(pv)
    pv[rand_seeds[valid]] <- remaining_seed_df$padj[valid]
    pr_vec <- page_rank(net_full, personalized = pv,
                         weights = E(net_full)$weight)$vector
    null_matrix[, perm_idx] <- pr_vec[target_nodes]
  }

  count_ge <- rowSums(null_matrix >= obs_scores, na.rm = TRUE)
  emp_pval <- (count_ge + 1) / (N_PERM + 1)

  left_out_pos  <- which(target_nodes == left_out_gene)
  left_out_pval <- emp_pval[left_out_pos]
  left_out_pr   <- obs_scores[left_out_pos]
  recovered     <- left_out_pval < FDR_ALPHA

  data.frame(
    Trait             = trait,
    n_seeds_total     = n_seeds,
    n_seeds_used      = nrow(remaining_seed_df),
    left_out_ENSG     = left_out_gene,
    left_out_gene     = left_out_name,
    left_out_padj     = left_out_padj,
    left_out_pagerank = left_out_pr,
    emp_pval_loo      = left_out_pval,
    recovered         = recovered,
    stringsAsFactors  = FALSE
  )
}

# == 병렬 실행: (disease, seed_idx) 쌍 단위 ====================================
# 기존: mclapply(diseases) → 내부 lapply(seed_idx)
#   문제: disease 수만큼만 코어 사용 (테스트 1개 질환 → 코어 1개만 사용)
#
# 수정: 모든 (disease, seed_idx) 쌍을 미리 목록으로 만들고 쌍 단위로 병렬화
#   효과: 전체 LOO 케이스 수만큼 코어 활용 가능
#   예) 225개 질환 × 평균 5 seed = 1,124개 작업 → 48코어 동시 처리

files <- sort(list.files(DATA_DIR, pattern = "\\.rds$", full.names = TRUE))
if (!is.null(TEST_SPECIFIC))
  files <- files[basename(files) %in% TEST_SPECIFIC]

# (filepath, seed_idx) 쌍 목록 생성
loo_jobs <- do.call(rbind, lapply(files, function(f) {
  node    <- as.data.frame(readRDS(f), stringsAsFactors = FALSE)
  node$padj <- as.numeric(node$padj)
  n_seeds <- sum(!is.na(node$padj) & node$padj != 0)
  if (n_seeds < 2) return(NULL)
  data.frame(filepath = f, seed_idx = seq_len(n_seeds),
             stringsAsFactors = FALSE)
}))

cat(sprintf("총 LOO 작업 수: %d  (%d diseases × seed 수)\n",
            nrow(loo_jobs), length(files)))

t0 <- proc.time()["elapsed"]
set.seed(42)

cat(sprintf("Starting parallel LOO analysis with %d cores...\n", N_CORES))

all_results <- pbmclapply(
  seq_len(nrow(loo_jobs)),
  function(i) {
    tryCatch(
      run_loo_single(loo_jobs$filepath[i], loo_jobs$seed_idx[i]),
      error = function(e) {
        message(sprintf("[ERROR] Task %d/%d - %s seed_idx=%d: %s",
                        i, nrow(loo_jobs),
                        basename(loo_jobs$filepath[i]),
                        loo_jobs$seed_idx[i],
                        conditionMessage(e)))
        NULL
      }
    )
  },
  mc.cores          = N_CORES,
  mc.set.seed       = TRUE,
  mc.style          = "ETA",     # Show ETA in progress bar
  mc.substyle       = 3,         # More detailed progress bar
  mc.silent         = FALSE      # Show progress messages
)

elapsed <- proc.time()["elapsed"] - t0
cat(sprintf("완료. 소요시간: %.0f초 (%.1f분)\n\n", elapsed, elapsed/60))

# == 결과 저장 =================================================================
all_results <- Filter(Negate(is.null), all_results)
if (length(all_results) == 0) stop("No results generated.")

final <- do.call(rbind, all_results)
write.csv(final, OUTPUT_FILE, row.names = FALSE)

# == 요약 =====================================================================
cat(strrep("=", 60), "\n")
cat(sprintf("총 %d개 질환, %d개 LOO 케이스\n",
            length(unique(final$Trait)), nrow(final)))
cat(sprintf("Recovery (emp_pval < %.2f): %d / %d (%.1f%%)\n",
            FDR_ALPHA,
            sum(final$recovered),
            nrow(final),
            sum(final$recovered) / nrow(final) * 100))

cat("\n질환별 recovery 요약:\n")
for (trait in unique(final$Trait)) {
  sub     <- final[final$Trait == trait, ]
  n_rec   <- sum(sub$recovered)
  n_total <- nrow(sub)
  cat(sprintf("  %-40s %d / %d recovered\n",
              gsub("finngen_R12_","",trait), n_rec, n_total))
}
cat(sprintf("\n결과 저장: %s\n", OUTPUT_FILE))