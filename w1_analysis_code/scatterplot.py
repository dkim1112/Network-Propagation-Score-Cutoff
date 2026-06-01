#!/usr/bin/env python3
"""
Scatter plots: propagation score vs -log10(empirical p-value), one panel per
disease, saved as a single combined PNG.

Per panel:
    - gray   = 일반 유전자
    - red    = propagation score 상위 TOP_FRAC (예: 5%)
    - blue solid line   = 전체 trend
    - red dashed line   = top-TOP_FRAC trend

Input:
  - result_np_cutoff/np_cutoff_results.csv  (propagation + emp_pval 결과)
    - result_np_cutoff/loo_results.csv        (optional; ignored)

Output: result_np_cutoff/scatter_prop_vs_emp_pval.png

사용법:
    python3 w3_analysis_code/scatter_prop_vs_emp_pval.py
"""

import csv
import math
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

INPUT_FILE  = "result_np_cutoff/np_cutoff_results.csv"
OUTPUT_FILE = "result_np_cutoff/scatter_prop_vs_emp_pval.png"
TOP_FRAC = 0.05
CUTOFFS  = [1.00, 0.50, 0.25, 0.10, 0.05, 0.01]
NCOLS    = 5

SCRIPT_DIR  = Path(__file__).resolve().parent
REPO_ROOT   = SCRIPT_DIR.parent
INPUT_PATH  = REPO_ROOT / INPUT_FILE
OUTPUT_PATH = REPO_ROOT / OUTPUT_FILE


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
    k = max(1, int(round(len(scores) * frac)))
    return sorted(range(len(scores)), key=lambda i: scores[i], reverse=True)[:k]


def load_loo_recovered():
    """LOO 결과에서 recovered=True인 유전자들을 trait별로 반환.
       파일 없으면 빈 dict 반환 (LOO 없이도 동작)"""
    # LOO highlighting removed — function kept for compatibility but returns empty
    return {}


def main():
    # LOO 결과 로드
    loo_recovered = load_loo_recovered()
    total_loo = sum(len(g) for g in loo_recovered.values())
    print(f"LOO recovered genes: {total_loo} (across {len(loo_recovered)} traits)")

    # ── propagation score + emp_pval + gene name도 같이 저장 ─────────────────
    traits = {}
    with open(INPUT_PATH, "r", newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            t = traits.setdefault(row["Trait"], {
                "page_rank": [], "neg_log10": [], "gene": [], "n_seeds": None
            })
            t["page_rank"].append(float(row["page.rank"]))
            t["neg_log10"].append(-math.log10(float(row["emp_pval"])))
            t["gene"].append(row.get("gene", ""))
            t["n_seeds"] = int(row["n_seeds"])

    n_traits = len(traits)
    ncols = min(NCOLS, n_traits)
    nrows = math.ceil(n_traits / ncols)
    fig, axes = plt.subplots(nrows, ncols,
                              figsize=(4 * ncols, 3.6 * nrows), squeeze=False)
    fig.suptitle(
        f"Propagation score vs -log10(empirical p-value)",
        fontsize=13,
    )

    for idx, (trait, t) in enumerate(sorted(traits.items())):
        ax = axes[idx // ncols][idx % ncols]
        x        = t["page_rank"]
        y        = t["neg_log10"]
        genes    = t["gene"]
        n        = len(x)
        n_seeds  = t["n_seeds"]
        top_idx  = set(top_frac_indices(x, TOP_FRAC))

        # ── Top-N (= seed 수) 인덱스 계산 ────────────────────────────────────
        topN_idx = set(sorted(range(n), key=lambda i: x[i], reverse=True)[:n_seeds])

        # ── 분류 ────────────────────────────────────────────────────────────
        gx = [x[i] for i in range(n) if i not in top_idx]
        gy = [y[i] for i in range(n) if i not in top_idx]
        tx = [x[i] for i in range(n) if i in top_idx]
        ty = [y[i] for i in range(n) if i in top_idx]

        ax.scatter(gx, gy, s=2, alpha=0.2, c="gray")
        ax.scatter(tx, ty, s=6, alpha=0.7, c="red")

        # LOO highlighting removed

        # ── trend line ───────────────────────────────────────────────────────
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

        clean = trait.replace("finngen_R12_", "").replace("_", " ")
        title = (f"{clean}\nseeds = {n_seeds},  n = {n}\n"
                 f"r(all)={r_all:.3f}  r(top {TOP_FRAC*100:g}%)={r_top:.3f}")
        if n_loo > 0:
            title += f"\nLOO: {loo_in_topN}/{n_loo} in top-N"
        ax.set_title(title, fontsize=8)
        ax.set_xlabel("propagation score (page.rank)", fontsize=8)
        ax.set_ylabel("-log10(emp_pval)", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.grid(True, alpha=0.3)

        # 유의성 기준선
        ax.axhline(y=-math.log10(0.05), color="orange", linestyle="--",
                   alpha=0.6, linewidth=1)
        ax.axhline(y=-math.log10(0.01), color="red",    linestyle="--",
                   alpha=0.6, linewidth=1)

    for idx in range(n_traits, nrows * ncols):
        axes[idx // ncols][idx % ncols].axis("off")

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    os.makedirs(OUTPUT_PATH.parent, exist_ok=True)
    fig.savefig(OUTPUT_PATH, dpi=200, bbox_inches="tight")
    plt.close(fig)

    print(f"\nSaved: {OUTPUT_PATH}")


if __name__ == "__main__":
    main()