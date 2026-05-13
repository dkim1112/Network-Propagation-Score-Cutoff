################################################################################
# Network Propagation Score Cutoff
# Method: Degree-preserving permutation based empirical p-value
#
# 사용법 (김동은 PC기준):
#   & "C:\Program Files\R\R-4.6.0\bin\Rscript.exe" "C:\Users\de7da\OneDrive\Desktop\NPCutOff\np_cutoff_permutation.R"
#
# 절차:
# 1. padj == 0인 유전자만 분석 대상 (seed 제외)
# 2. seed와 동일한 수 + degree-bin matched random seeds로 N_PERM회 permutation
#    → 매번 page_rank 재실행 → 유전자별 null distribution 생성
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
N_BINS             <- 10
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

# == 2. Degree bin 사전 계산 ===================================================
net_deg_vec <- igraph::degree(net_full)
net_deg_df  <- data.frame(node   = names(net_deg_vec),
                           degree = as.numeric(net_deg_vec),
                           stringsAsFactors = FALSE)
breaks <- unique(quantile(net_deg_df$degree,
                           probs = seq(0, 1, length.out = N_BINS + 1)))
net_deg_df$bin <- cut(net_deg_df$degree, breaks = breaks,
                       labels = FALSE, include.lowest = TRUE)

# bin별 후보 노드 pre-index
bin_candidates <- split(net_deg_df$node, net_deg_df$bin)

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

  # 각 seed의 degree bin
  seed_bins <- sapply(seed_df$ENSG, function(g) {
    if (g %in% net_deg_df$node) {
      net_deg_df$bin[net_deg_df$node == g]
    } else {
      net_deg_df$bin[which.min(abs(net_deg_df$degree -
                                     median(net_deg_df$degree)))[1]]
    }
  })

  # == Permutation =============================================================
  null_matrix <- matrix(NA_real_, nrow = length(target_nodes), ncol = N_PERM)

  # degree-bin matched random seed 넣고 N_PERM회 실행
  for (perm_idx in seq_len(N_PERM)) {

    # degree-bin matched random seed 선택
    rand_seeds <- character(n_seeds)
    for (s in seq_len(n_seeds)) {
      cands <- setdiff(bin_candidates[[seed_bins[s]]], seed_df$ENSG)
      if (length(cands) == 0) cands <- net_deg_df$node
      rand_seeds[s] <- cands[sample.int(length(cands), 1)]
    }

    # personalization vector
    pv           <- base_pv
    valid        <- rand_seeds %in% names(pv)
    pv[rand_seeds[valid]] <- seed_df$padj[valid]

    # page_rank 재실행 = RWR 재실행
    # net_full = 네트워크 원본 (Combined_STRINGv11_OTAR281119_FILTER.rds)
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

all_results <- vector("list", n_files)
skipped     <- character()

for (i in seq_along(files)) {
  fname <- basename(files[i])
  tryCatch({
    res <- process_disease(files[i])
    if (!is.null(res)) {
      all_results[[i]] <- res
    } else {
      skipped <- c(skipped, fname)
    }
  }, error = function(e) {
    skipped <<- c(skipped, fname)
  })
}

# == 5. 결과 저장 ==============================================================
all_results <- Filter(Negate(is.null), all_results)
if (length(all_results) == 0) stop("No results generated.")

final <- do.call(rbind, all_results)
write.csv(final, OUTPUT_FILE, row.names = FALSE)

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