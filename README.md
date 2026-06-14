# Network Propagation Score Cutoff

A statistical framework for deciding **which genes count as "significant"** after running network propagation (random-walk-with-restart / personalized PageRank) over a protein–protein interaction (PPI) network for disease gene prioritization.

Instead of picking an arbitrary top-N or raw-score threshold, this project assigns each gene a **permutation-based empirical p-value** using degree-matched random seed sets, and then validates that the resulting cutoff is biologically meaningful against external gold standards (OpenTargets, ChEMBL/clinical precedence).

The disease seed sets come from **FinnGen R12** GWAS endpoints; the network is a combined **STRING v11 + OTAR** interactome.

---

## Core idea

For each disease:

1. **Seeds** are known disease genes (rows where `padj != 0`); their `padj` becomes the personalization weight.
2. Propagation (igraph `page_rank`, weighted, personalized) spreads seed signal across the network. Every non-seed gene (`padj == 0`) gets a propagation score (`page.rank`).
3. To ask *"is this gene's score higher than expected by chance?"*, the seeds are replaced `N_PERM` times by **degree-matched random gene sets** (nearest-neighbor degree matching, not quantile bins), re-running propagation each time to build a per-gene null distribution.
4. **Empirical p-value:** `p = (#{null >= observed} + 1) / (N_PERM + 1)`. Genes with `emp_pval < 0.05` are called significant.

> Note on FDR: BH correction is intentionally **not** applied. At `N_PERM = 1000` the minimum possible p-value (0.001) times ~18,400 genes ≈ 18.4, so nothing could survive multiple-testing correction mathematically. The empirical p-value is used directly as the cutoff. Genes from diseases with very few seeds (`n_seeds <= 2`) are flagged `low_confidence`.

---

## Pipeline

```
result_network_propagation/*.rds   (per-disease propagation results: ENSG, gene, Trait, padj, page.rank, degree)
tables_expansion/*.rds             (Combined STRING v11 + OTAR PPI network)
            │
            ▼
   np_cutoff_permutation[_parallel].R   ──►  result_np_cutoff/np_cutoff_results.csv
            │                                  (emp_pval, significant, n_seeds, low_confidence)
            ▼
   w1–w5 analysis / validation
```

### Main scripts (project root)

| Script | Purpose |
|---|---|
| `np_cutoff_permutation.R` | Reference single-threaded implementation (`N_PERM = 1000`). Windows-oriented paths. |
| `np_cutoff_permutation_parallel.R` | Server version. Parallelized over both diseases and permutations via `mclapply` (auto-detects cores, splits disease/permutation cores). `N_PERM = 10000`. Use this for real runs. |

Both share the same logic and the `make_nearest_match_sampler()` degree-matching routine. The list of diseases to process is controlled by the `TEST_SPECIFIC` variable (set to `NULL` to run all ~225 diseases).

### Analysis stages (`w1`–`w5`)

Each `wN_analysis_code/` folder is a sequential validation step. Most consume `result_np_cutoff/np_cutoff_results.csv` (or a downstream CSV) and write their own results + plots in place.

**w1 — Score ↔ significance relationship**
- `correlation_analysis.py` → per-disease Pearson correlation between `page.rank` and `-log10(emp_pval)` (`disease_correlations.csv`).
- `seed_correlation.py` → tests whether low seed counts weaken that correlation (reliability vs. seed count).
- `scatterplot.py` → per-disease scatter of propagation score vs `-log10(emp_pval)`, overlaying OT / ChEMBL gold standards.
- `verify_linearity.R` → tests the linearity hypothesis `NP(A+B+…) ≈ NP(A)+NP(B)+…`. If true, propagation could be run once per unique seed instead of once per disease (large compute saving).

**w2 — Method validation (Leave-One-Out)**
- `LOO_analysis.R` → drops each seed one at a time, re-runs propagation + permutation on the remaining seeds, and checks whether the held-out gene is recovered (`emp_pval < 0.05`). Parallelized per (disease, seed) pair. Output: `result_np_cutoff/loo_results.csv`.
- `pvalue_precision_analysis.py` → compares LOO empirical p-values at 1K vs 10K permutations to quantify how much extra permutations sharpen p-values.

**w3 — Gold-standard cutoff selection & validation**
- `gold_standard.py` → pulls OpenTargets associations (GraphQL API) per disease (EFO/MONDO mapping), builds gold sets from `globalScore` and `clinical_precedence` (ChEMBL), then computes TP / FN / Novel / precision / sensitivity for significant genes. Output: `gold_standard_validation.csv`, `gold_standard_summary.csv`.
- `choose_cutoff.py` → sweeps OT and clinical-precedence cutoff combinations and plots the resulting precision/sensitivity tradeoff + score histograms.
- `topn_vs_emp_comparison.py` → compares the empirical-p-value cutoff against a simple top-N (N = seed count) propagation cutoff, splitting gold genes into both/only-emp/only-topn/missed and comparing false negatives.

**w4 — Propagation rank add-on**
- `loo_add_prop_rank.R` → augments existing LOO results with the held-out gene's propagation rank (single propagation run, no permutation), to compare recovery by empirical p-value vs. recovery by raw propagation rank.

**w5 — Quadrant enrichment**
- `quadrant_enrichment.py` → splits genes into four quadrants by (propagation score quantile × significance) and measures gold-standard density per quadrant. Validates the cutoff: the "low propagation but significant" quadrant should still be gold-enriched, while "low + non-significant" acts as a negative control.
- `top_left_permutation_test.py` → hypergeometric enrichment test (pooled + per-trait) confirming whether each quadrant's gold count exceeds chance.

---

## Data

| Path | Contents | In git? |
|---|---|---|
| `result_network_propagation/` | ~225 per-disease `.rds` files (`nodes.finngen_R12_*.rds`). Each has columns `ENSG`, `gene`, `Trait`, `padj`, `page.rank`, `degree`. | No (gitignored — large) |
| `tables_expansion/Combined_STRINGv11_OTAR281119_FILTER.rds` | Weighted PPI edge list: `ENSG_A`, `ENSG_B`, `combined_score`. | Yes |
| `result_np_cutoff/` | Pipeline outputs: `np_cutoff_results.csv`, `loo_results.csv`, scatter plots, correlations. | Partial |

`result_np_cutoff/np_cutoff_results.csv` is the central output table consumed by nearly every downstream analysis.

---

## Running it

### Prerequisites
- **R** (4.x) with `igraph`, `parallel`, and `pbmcapply` (the LOO script auto-installs `pbmcapply`).
- **Python 3** with `pandas`, `numpy`, `scipy`, `matplotlib`, `requests`. A `venv/` is present in the repo.

### Typical order (run from project root)

```bash
# 1. Compute empirical p-value cutoffs (server, parallel)
Rscript np_cutoff_permutation_parallel.R          # → result_np_cutoff/np_cutoff_results.csv

# 2. Score/significance relationship
python3 w1_analysis_code/correlation_analysis.py
python3 w1_analysis_code/seed_correlation.py

# 3. Method validation
Rscript w2_analysis_code/LOO_analysis.R           # → result_np_cutoff/loo_results.csv
python3 w2_analysis_code/pvalue_precision_analysis.py

# 4. Gold-standard validation (needs internet for OpenTargets API)
python3 w3_analysis_code/gold_standard.py
python3 w3_analysis_code/choose_cutoff.py
python3 w3_analysis_code/topn_vs_emp_comparison.py
python3 w1_analysis_code/scatterplot.py

# 5. Propagation-rank add-on
Rscript w4_analysis_code/loo_add_prop_rank.R

# 6. Quadrant enrichment
python3 w5_analysis_code/quadrant_enrichment.py
python3 w5_analysis_code/top_left_permutation_test.py
```

The R scripts auto-detect whether they're run from the project root or a subdirectory. Long server runs are typically launched with `nohup ... &` (see the header comment in `LOO_analysis.R`).

### Key parameters
- `N_PERM` — permutations per disease (1000 for quick tests, 10000 for final).
- `TEST_SPECIFIC` — restrict to a list of diseases, or `NULL` for all.
- `LOW_SEED_THRESHOLD` (default 2), `FDR_ALPHA` (default 0.05), `RANDOM_SEED` (42).
- In `w3` scripts: `OT_SCORE_CUTOFF`, `CLINICAL_CUTOFF`, `COMBINE_MODE` (OR/AND) define the gold standard.

---

## Notes
- Many in-code comments are in Korean; the analysis logic and outputs are language-independent.
- The OpenTargets-dependent scripts (`w3`) require network access to `api.platform.opentargets.org`.
- Reproducibility relies on `set.seed(42)` plus `mc.set.seed = TRUE` in the parallel paths.
