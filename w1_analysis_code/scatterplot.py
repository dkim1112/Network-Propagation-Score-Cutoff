#!/usr/bin/env python3
"""
Scatter plots: propagation score vs -log10(empirical p-value), one panel per
disease, saved as a single combined PNG.

Goal: the paper highlights the top-N proteins by propagation score; this checks
how well that ranking agrees with the newly introduced empirical p-value, and
in particular whether the agreement is tighter at the extreme (top) of the
propagation-score distribution.

Per panel:
  - x = page.rank (propagation score), y = -log10(emp_pval)
  - red    = top TOP_FRAC by propagation score (the proteins the paper emphasizes)
  - gray   = the rest
  - blue solid line  = trend over all genes
  - red dashed line  = trend over the top-TOP_FRAC genes
  - title shows the disease, its seed count, and Pearson r for all genes vs.
    for the top-TOP_FRAC genes.

Console also prints Pearson r at several score cutoffs (top 50/25/10/5/1%) so
the "are the extremes more correlated?" question gets a numeric answer too.

Input : result_np_cutoff/np_cutoff_results.csv  (produced by np_cutoff_permutation.R)
Output: result_np_cutoff/scatter_prop_vs_emp_pval.png
"""

import csv
import math
import os

import matplotlib
matplotlib.use("Agg")  # non-interactive backend
import matplotlib.pyplot as plt
import numpy as np

INPUT_FILE = "result_np_cutoff/np_cutoff_results.csv"
OUTPUT_FILE = "result_np_cutoff/scatter_prop_vs_emp_pval.png"
TOP_FRAC = 0.05                       # "top n" used in the paper
CUTOFFS = [1.00, 0.50, 0.25, 0.10, 0.05, 0.01]  # score quantiles to report r at
NCOLS = 5


def pearson_correlation(x, y):
    n = len(x)
    if n < 2:
        return float("nan")
    mean_x = sum(x) / n
    mean_y = sum(y) / n
    num = sum((x[i] - mean_x) * (y[i] - mean_y) for i in range(n))
    den = math.sqrt(sum((v - mean_x) ** 2 for v in x) * sum((v - mean_y) ** 2 for v in y))
    return num / den if den else float("nan")


def top_frac_indices(scores, frac):
    """Indices of the highest-scoring `frac` fraction of `scores`."""
    k = max(1, int(round(len(scores) * frac)))
    return sorted(range(len(scores)), key=lambda i: scores[i], reverse=True)[:k]


def main():
    traits = {}  # trait -> dict(page_rank=[], neg_log10=[], n_seeds=int)
    with open(INPUT_FILE, "r", newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            t = traits.setdefault(row["Trait"], {"page_rank": [], "neg_log10": [], "n_seeds": None})
            t["page_rank"].append(float(row["page.rank"]))
            t["neg_log10"].append(-math.log10(float(row["emp_pval"])))
            t["n_seeds"] = int(row["n_seeds"])

    n_traits = len(traits)

    ncols = min(NCOLS, n_traits)
    nrows = math.ceil(n_traits / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3.6 * nrows), squeeze=False)
    fig.suptitle(
        f"Propagation score vs -log10(empirical p-value)  "
        f"(red = top {TOP_FRAC*100:g}% by score)",
        fontsize=14,
    )

    r_all_list, r_top_list = [], []
    for idx, (trait, t) in enumerate(sorted(traits.items())):
        ax = axes[idx // ncols][idx % ncols]
        x = t["page_rank"]
        y = t["neg_log10"]
        n = len(x)
        top_idx = set(top_frac_indices(x, TOP_FRAC))

        gx = [x[i] for i in range(n) if i not in top_idx]
        gy = [y[i] for i in range(n) if i not in top_idx]
        tx = [x[i] for i in range(n) if i in top_idx]
        ty = [y[i] for i in range(n) if i in top_idx]

        ax.scatter(gx, gy, s=2, alpha=0.2, c="gray")
        ax.scatter(tx, ty, s=6, alpha=0.7, c="red")

        if n > 1:
            a, b = np.polyfit(x, y, 1)
            xs = np.linspace(min(x), max(x), 100)
            ax.plot(xs, a * xs + b, "b-", lw=1.2, alpha=0.8)
        if len(tx) > 1:
            a, b = np.polyfit(tx, ty, 1)
            xs = np.linspace(min(tx), max(tx), 50)
            ax.plot(xs, a * xs + b, "r--", lw=1.5, alpha=0.9)

        r_all = pearson_correlation(x, y)
        r_top = pearson_correlation(tx, ty)
        r_all_list.append(r_all)
        r_top_list.append(r_top)

        clean = trait.replace("finngen_R12_", "").replace("_", " ")
        ax.set_title(
            f"{clean}\nseeds = {t['n_seeds']},  n = {n}\n"
            f"r(all) = {r_all:.3f},  r(top {TOP_FRAC*100:g}%) = {r_top:.3f}",
            fontsize=8.5,
        )
        ax.set_xlabel("propagation score (page.rank)", fontsize=8)
        ax.set_ylabel("-log10(emp_pval)", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.grid(True, alpha=0.3)

        # Add P-value significance lines
        ax.axhline(y=-math.log10(0.05), color='orange', linestyle='--', alpha=0.7, linewidth=1, label='P=0.05')
        ax.axhline(y=-math.log10(0.01), color='red', linestyle='--', alpha=0.7, linewidth=1, label='P=0.01')


    for idx in range(n_traits, nrows * ncols):
        axes[idx // ncols][idx % ncols].axis("off")

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)
    fig.savefig(OUTPUT_FILE, dpi=200, bbox_inches="tight")
    plt.close(fig)



if __name__ == "__main__":
    main()
