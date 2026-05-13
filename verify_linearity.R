################################################################################
# 가설 검증: 여러 seed 동시 실행 vs 개별 seed 따로 실행 후 합산
#
# 가설: NP(seed_A + seed_B + ...) ≈ NP(seed_A) + NP(seed_B) + ...
# 이때, seed node: 질환관련유전자
#
# 만약 가설이 맞다면:
#   - 질환마다 NP를 돌리는 대신 seed 단백질별로 한 번씩만 돌리고
#   - 질환별로 해당 seed들의 결과를 합산하면 됨
#   - 329개 고유 seed만 돌리면 되므로 1,124회 → 329회로 절감 (70%)
#
# 검증 방법:
#   1. 질환별로 combined 방식(현재)과 individual 합산 방식을 모두 실행
#   2. 두 결과의 pearson correlation 및 rank correlation 비교
#   3. 상위 N개 유전자 overlap 비교 (실제 분석에서 중요한 지표)
#
# Run with: Rscript verify_linearity.R
################################################################################

user_lib <- file.path(Sys.getenv("APPDATA"), "R", "win-library", "4.6")
if (dir.exists(user_lib)) {
  .libPaths(unique(c(user_lib, .libPaths())))
}

library(igraph)

# == 파라미터 ==================================================================
DATA_DIR     <- "result_network_propagation"
NETWORK_FILE <- "tables_expansion/Combined_STRINGv11_OTAR281119_FILTER.rds"
OUTPUT_FILE  <- "result_np_cutoff/linearity_verification.csv"
TOP_N        <- c(50, 100, 200, 500)   # overlap 비교할 상위 N개

# 검증할 질환 (seed 수 다양하게)
TEST_DISEASES <- c(
  "nodes.finngen_R12_ABDOM_HERNIA.rds",           # seed 2
  "nodes.finngen_R12_L12_ATOPIC.rds",             # seed 3
  "nodes.finngen_R12_T1D.rds",                    # seed 5
  "nodes.finngen_R12_AUTOIMMUNE_NONTHYROID.rds",  # seed 7
  "nodes.finngen_R12_T2D_WIDE.rds",               # seed 9
  "nodes.finngen_R12_AUTOIMMUNE.rds"              # seed 13
)
RANDOM_SEED <- 42
# =============================================================================

set.seed(RANDOM_SEED)
dir.create(dirname(OUTPUT_FILE), showWarnings = FALSE, recursive = TRUE)

cat("=== 가설 검증: Combined vs Individual seed propagation ===\n\n")

# == 1. 네트워크 로드 ==========================================================
cat("Loading PPI network...\n")
string   <- as.data.frame(readRDS(NETWORK_FILE), stringsAsFactors = FALSE)
string$combined_score <- as.numeric(string$combined_score)
net_full <- graph_from_data_frame(d = string[, c("ENSG_A","ENSG_B")],
                                   directed = FALSE)
E(net_full)$weight <- string$combined_score
net_full <- igraph::simplify(net_full, remove.loops = TRUE,
                              remove.multiple = TRUE,
                              edge.attr.comb = c(weight="max","ignore"))
all_nodes <- V(net_full)$name
base_pv   <- setNames(rep(0.0, length(all_nodes)), all_nodes)
cat(sprintf("Network: %d nodes, %d edges\n\n", length(all_nodes), ecount(net_full)))

# == 2. 질환별 검증 ============================================================
results <- list()

for (fname in TEST_DISEASES) {
  fpath <- file.path(DATA_DIR, fname)
  cat(sprintf("--- %s ---\n", fname))

  node     <- as.data.frame(readRDS(fpath), stringsAsFactors = FALSE)
  node$padj      <- as.numeric(node$padj)
  node$page.rank <- as.numeric(node$page.rank)

  trait    <- unique(node$Trait)[1]
  seed_df  <- node[!is.na(node$padj) & node$padj != 0, ]
  target_df <- node[!is.na(node$padj) & node$padj == 0, ]
  n_seeds  <- nrow(seed_df)
  target_ensg <- target_df$ENSG[target_df$ENSG %in% all_nodes]

  cat(sprintf("  Seeds: %d, Targets: %d\n", n_seeds, length(target_ensg)))

  # 방법 A: Combined (모든 seed 동시에)
  pv_all <- base_pv
  pv_all[seed_df$ENSG[seed_df$ENSG %in% all_nodes]] <-
    seed_df$padj[seed_df$ENSG %in% all_nodes]
  pr_combined <- page_rank(net_full, personalized = pv_all,
                            weights = E(net_full)$weight)$vector
  scores_combined <- pr_combined[target_ensg]

  # 방법 B: Individual 합산/평균/중앙값 (seed 하나씩 따로 실행)
  pr_list <- list()
  for (i in seq_len(n_seeds)) {
    ensg_i <- seed_df$ENSG[i]
    if (!ensg_i %in% all_nodes) next
    pv_single <- base_pv
    pv_single[ensg_i] <- seed_df$padj[i]
    pr_i <- page_rank(net_full, personalized = pv_single,
                       weights = E(net_full)$weight)$vector
    pr_list[[ensg_i]] <- pr_i
  }

  # sum / mean / median 계산
  pr_mat       <- do.call(rbind, lapply(pr_list, function(x) x[target_ensg]))
  scores_sum    <- colSums(pr_mat)
  scores_mean   <- colMeans(pr_mat)
  scores_median <- apply(pr_mat, 2, median)

  # == 값 자체 비교 (박사님 가설의 핵심) ======================================
  # combined vs sum/mean/median 이 얼마나 같은지
  # 1) 정규화 후 비교 (스케일이 다르므로 0-1로 맞춤)
  norm <- function(x) (x - min(x)) / (max(x) - min(x))
  nc  <- norm(scores_combined)
  ns  <- norm(scores_sum)
  nm  <- norm(scores_mean)
  nmd <- norm(scores_median)

  mad_sum    <- mean(abs(nc - ns))    # mean absolute difference
  mad_mean   <- mean(abs(nc - nm))
  mad_median <- mean(abs(nc - nmd))

  cat(sprintf("  [값 비교 - 정규화 후 MAD]\n"))
  cat(sprintf("    combined vs individual_sum:    %.4f\n", mad_sum))
  cat(sprintf("    combined vs individual_mean:   %.4f\n", mad_mean))
  cat(sprintf("    combined vs individual_median: %.4f\n", mad_median))

  # 2) rank correlation
  r_rank_sum    <- cor(rank(-scores_combined), rank(-scores_sum),    method="spearman")
  r_rank_mean   <- cor(rank(-scores_combined), rank(-scores_mean),   method="spearman")
  r_rank_median <- cor(rank(-scores_combined), rank(-scores_median), method="spearman")

  cat(sprintf("  [순위 비교 - r]\n"))
  cat(sprintf("    combined vs individual_sum:    %.4f\n", r_rank_sum))
  cat(sprintf("    combined vs individual_mean:   %.4f\n", r_rank_mean))
  cat(sprintf("    combined vs individual_median: %.4f\n", r_rank_median))

  r_pearson <- cor(scores_combined, scores_sum)
  r_rank    <- r_rank_sum  # 대표값으로 sum 사용

  # 상위 N개 overlap (sum/mean/median 각각)
  overlap_results <- list()
  for (n in TOP_N) {
    top_combined <- names(sort(scores_combined, decreasing=TRUE))[1:n]
    for (method in c("sum","mean","median")) {
      sc <- switch(method, sum=scores_sum, mean=scores_mean, median=scores_median)
      top_m   <- names(sort(sc, decreasing=TRUE))[1:n]
      overlap <- length(intersect(top_combined, top_m))
      overlap_results[[paste0(method,"_",n)]] <- round(overlap/n*100, 1)
    }
  }

  cat(sprintf("  [Top N overlap vs combined]\n"))
  cat(sprintf("  %10s %8s %8s %8s %8s\n", "", "Top50", "Top100", "Top200", "Top500"))
  for (method in c("sum","mean","median")) {
    cat(sprintf("  %10s %7.1f%% %7.1f%% %7.1f%% %7.1f%%\n", method,
                overlap_results[[paste0(method,"_50")]],
                overlap_results[[paste0(method,"_100")]],
                overlap_results[[paste0(method,"_200")]],
                overlap_results[[paste0(method,"_500")]]))
  }
  cat("\n")

  results[[trait]] <- data.frame(
    Trait              = trait,
    n_seeds            = n_seeds,
    # 값 자체 비교 (MAD, 정규화 후)
    MAD_sum            = round(mad_sum,    4),
    MAD_mean           = round(mad_mean,   4),
    MAD_median         = round(mad_median, 4),
    # 순위 비교
    spearman_sum       = round(r_rank_sum,    4),
    spearman_mean      = round(r_rank_mean,   4),
    spearman_median    = round(r_rank_median, 4),
    # top N overlap (sum 기준)
    overlap_top50_sum  = overlap_results[["sum_50"]],
    overlap_top100_sum = overlap_results[["sum_100"]],
    overlap_top200_sum = overlap_results[["sum_200"]],
    overlap_top500_sum = overlap_results[["sum_500"]]
  )
}

# == 3. 결과 저장 ==============================================================
final <- do.call(rbind, results)
write.csv(final, OUTPUT_FILE, row.names = FALSE)

cat(strrep("=", 60), "\n")
cat(sprintf("\n결과 저장: %s\n", OUTPUT_FILE))