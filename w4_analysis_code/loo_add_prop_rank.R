################################################################################
# Propagation Rank 추가 분석 (기존 LOO 결과 활용)
#
# 기존 loo_results.csv에 emp_pval은 이미 있음 → 그대로 사용
# 추가로 필요한 것: 같은 LOO 조건(n-1 seed)에서 propagation_rank (전에 loo_results.csv 돌릴땐 emp_pval 계산하느라 propagation꺼는 안 했음)
#
# Permutation은 안 함 → page_rank 1회만 실행 → 매우 빠름 (전체 ~1분)
#
# Input:  result_np_cutoff/loo_results.csv  (기존 결과)
# Output: w4_analysis_code/loo_add_prop_rank.csv  (기존 + propagation_rank 추가)
################################################################################

library(igraph)

# == 파라미터 ==================================================================
LOO_FILE     <- "result_np_cutoff/loo_results.csv"  # 기존 LOO 결과 (N_PERM = 10000)
DATA_DIR     <- "result_network_propagation"
NETWORK_FILE <- "tables_expansion/Combined_STRINGv11_OTAR281119_FILTER.rds"
OUTPUT_FILE  <- "w4_analysis_code/loo_add_prop_rank.csv"
FDR_ALPHA    <- 0.05
# =============================================================================

dir.create(dirname(OUTPUT_FILE), showWarnings = FALSE, recursive = TRUE)

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
all_nodes <- V(net_full)$name
base_pv   <- setNames(rep(0.0, length(all_nodes)), all_nodes)
cat(sprintf("Network: %d nodes\n\n", length(all_nodes)))

# == 기존 LOO 결과 로드 =======================================================
loo_df <- read.csv(LOO_FILE, stringsAsFactors = FALSE)
cat(sprintf("기존 LOO 결과: %d cases | %d traits\n\n",
            nrow(loo_df), length(unique(loo_df$Trait))))

# == 질환별 rds 파일 캐시 ======================================================
trait_to_filename <- function(trait) {
  sprintf("%s/nodes.%s.rds", DATA_DIR, trait)
}

rds_cache <- list()
for (trait in unique(loo_df$Trait)) {
  fname <- trait_to_filename(trait)
  if (!file.exists(fname)) {
    stop(sprintf("파일 없음: %s", fname))
  }
  node <- as.data.frame(readRDS(fname), stringsAsFactors = FALSE)
  node$padj <- as.numeric(node$padj)
  rds_cache[[trait]] <- node
}

# == 각 LOO case에 propagation_rank 추가 =======================================
t0 <- proc.time()["elapsed"]
prop_ranks       <- numeric(nrow(loo_df))
prop_ranks_total <- numeric(nrow(loo_df))

for (i in seq_len(nrow(loo_df))) {
  trait         <- loo_df$Trait[i]
  left_out_ENSG <- loo_df$left_out_ENSG[i]

  node       <- rds_cache[[trait]]
  all_seed_df <- node[!is.na(node$padj) & node$padj != 0, ]
  remaining_seed_df <- all_seed_df[all_seed_df$ENSG != left_out_ENSG, ]

  # n-1 seed로 propagation 1회 실행
  pv <- base_pv
  valid <- remaining_seed_df$ENSG %in% names(pv)
  pv[remaining_seed_df$ENSG[valid]] <- remaining_seed_df$padj[valid]
  pr <- page_rank(net_full, personalized = pv,
                  weights = E(net_full)$weight)$vector

  # 남은 seed 제외한 유전자들 중에서 left_out의 rank (공정한 비교)
  pr_excluded <- pr[!names(pr) %in% remaining_seed_df$ENSG]
  prop_ranks[i] <- sum(pr_excluded > pr[left_out_ENSG]) + 1

  # 전체 유전자 중 rank (참고)
  prop_ranks_total[i] <- sum(pr > pr[left_out_ENSG]) + 1

  if (i %% 5 == 0) {
    cat(sprintf("  %d / %d done\n", i, nrow(loo_df)))
  }
}

elapsed <- proc.time()["elapsed"] - t0
cat(sprintf("\n완료: %.1f초\n\n", elapsed))

# == 결과 합치기 ===============================================================
loo_df$propagation_rank       <- prop_ranks
loo_df$propagation_rank_total <- prop_ranks_total

# 두 방법 비교용 컬럼
loo_df$recovered_by_emp  <- loo_df$emp_pval_loo < FDR_ALPHA
loo_df$recovered_by_prop <- loo_df$propagation_rank <= loo_df$n_seeds_total

loo_df$both_recovered <- loo_df$recovered_by_emp & loo_df$recovered_by_prop
loo_df$only_emp       <- loo_df$recovered_by_emp & !loo_df$recovered_by_prop
loo_df$only_prop      <- !loo_df$recovered_by_emp & loo_df$recovered_by_prop
loo_df$neither        <- !loo_df$recovered_by_emp & !loo_df$recovered_by_prop

write.csv(loo_df, OUTPUT_FILE, row.names = FALSE)

# == 요약 출력 =================================================================
cat(strrep("=", 70), "\n")
cat("두 방법 recovery 비교:\n\n")

for (trait in unique(loo_df$Trait)) {
  sub <- loo_df[loo_df$Trait == trait, ]
  n   <- nrow(sub)
  cat(sprintf("%s  (n=%d)\n", gsub("finngen_R12_","",trait), n))
  cat(sprintf("  emp_pval recovered (< %.2f):     %2d / %d\n",
              FDR_ALPHA, sum(sub$recovered_by_emp), n))
  cat(sprintf("  prop recovered (rank <= N=%d): %2d / %d\n",
              sub$n_seeds_total[1], sum(sub$recovered_by_prop), n))
  cat(sprintf("  both:      %2d / %d\n", sum(sub$both_recovered), n))
  cat(sprintf("  only_emp:  %2d / %d\n", sum(sub$only_emp),  n))
  cat(sprintf("  only_prop: %2d / %d\n", sum(sub$only_prop), n))
  cat(sprintf("  neither:   %2d / %d\n\n", sum(sub$neither), n))
}

cat(sprintf("결과 저장: %s\n", OUTPUT_FILE))