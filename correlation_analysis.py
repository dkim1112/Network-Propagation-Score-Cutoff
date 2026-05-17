#!/usr/bin/env python3
"""
Per-disease correlation between propagation score and -log10(empirical p-value).

Input : result_np_cutoff/np_cutoff_results.csv  (produced by np_cutoff_permutation.R)
Output: result_np_cutoff/disease_correlations.csv
"""

import csv
import math

INPUT_FILE = "result_np_cutoff/np_cutoff_results.csv"
OUTPUT_FILE = "result_np_cutoff/disease_correlations.csv"


def pearson_correlation(x, y):
    """Pearson correlation coefficient between two equal-length lists."""
    if len(x) != len(y):
        raise ValueError("Lists must have same length")
    n = len(x)
    if n < 2:
        return float("nan")

    mean_x = sum(x) / n
    mean_y = sum(y) / n
    numerator = sum((x[i] - mean_x) * (y[i] - mean_y) for i in range(n))
    sum_sq_x = sum((x[i] - mean_x) ** 2 for i in range(n))
    sum_sq_y = sum((y[i] - mean_y) ** 2 for i in range(n))
    denominator = math.sqrt(sum_sq_x * sum_sq_y)
    if denominator == 0:
        return float("nan")
    return numerator / denominator


def main():
    traits = {}
    with open(INPUT_FILE, "r", newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            page_rank = float(row["page.rank"])
            neg_log10_emp_pval = -math.log10(float(row["emp_pval"]))
            traits.setdefault(row["Trait"], []).append((page_rank, neg_log10_emp_pval))

    results = []
    for trait, pairs in sorted(traits.items()):
        page_ranks = [p for p, _ in pairs]
        neg_log10_emp_pvals = [v for _, v in pairs]
        corr = pearson_correlation(page_ranks, neg_log10_emp_pvals)
        results.append({"Trait": trait, "corr_prop_emp_pval": corr, "n_genes": len(pairs)})

    with open(OUTPUT_FILE, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["Trait", "corr_prop_emp_pval", "n_genes"])
        writer.writeheader()
        for r in results:
            writer.writerow(r)


if __name__ == "__main__":
    main()
