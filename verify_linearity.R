################################################################################
# 가설 검증: 여러 seed 동시 실행 vs 개별 seed 따로 실행 후 합산
#
# 가설: NP(seed_A + seed_B + ...) ≈ NP(seed_A) + NP(seed_B) + ...
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

  # 방법 B: Individual 합산 (seed 하나씩 따로 실행)
  pr_sum <- setNames(rep(0.0, length(all_nodes)), all_nodes)
  for (i in seq_len(n_seeds)) {
    ensg_i <- seed_df$ENSG[i]
    if (!ensg_i %in% all_nodes) next
    pv_single <- base_pv
    pv_single[ensg_i] <- seed_df$padj[i]
    pr_i <- page_rank(net_full, personalized = pv_single,
                       weights = E(net_full)$weight)$vector
    pr_sum <- pr_sum + pr_i
  }
  scores_sum <- pr_sum[target_ensg]

  # == 비교 ===================================================================
  r_pearson <- cor(scores_combined, scores_sum)
  r_rank    <- cor(rank(-scores_combined), rank(-scores_sum),
                   method = "spearman")

  cat(sprintf("  Pearson r (value): %.4f\n", r_pearson))
  cat(sprintf("  Spearman r (rank): %.4f\n", r_rank))

  # 상위 N개 overlap
  overlap_results <- list()
  for (n in TOP_N) {
    top_combined <- names(sort(scores_combined, decreasing=TRUE))[1:n]
    top_sum      <- names(sort(scores_sum,      decreasing=TRUE))[1:n]
    overlap      <- length(intersect(top_combined, top_sum))
    pct          <- overlap / n * 100
    cat(sprintf("  Top %4d overlap: %d / %d (%.1f%%)\n", n, overlap, n, pct))
    overlap_results[[as.character(n)]] <- pct
  }
  cat("\n")

  results[[trait]] <- data.frame(
    Trait          = trait,
    n_seeds        = n_seeds,
    r_pearson      = r_pearson,
    r_rank_spearman = r_rank,
    overlap_top50  = overlap_results[["50"]],
    overlap_top100 = overlap_results[["100"]],
    overlap_top200 = overlap_results[["200"]],
    overlap_top500 = overlap_results[["500"]]
  )
}

# == 3. 결과 저장 ==============================================================
final <- do.call(rbind, results)
write.csv(final, OUTPUT_FILE, row.names = FALSE)

cat(strrep("=", 55), "\n")
cat("요약:\n")
print(final[, c("Trait","n_seeds","r_pearson","r_rank_spearman",
                "overlap_top100","overlap_top500")],
      row.names = FALSE)
cat(sprintf("\n결과 저장: %s\n", OUTPUT_FILE))