#!/usr/bin/env python3
"""
Scatter plots: propagation score vs -log10(empirical p-value), one panel per
disease, saved as combined PNG files.

★ 3개의 그래프 생성 ★
  - scatter_prop_vs_emp_pval_ot.png        (OT globalScore 기준 gold standard 표시)
  - scatter_prop_vs_emp_pval_chembl.png    (ChEMBL/clinical_precedence 기준 gold standard 표시)
  - scatter_prop_vs_emp_pval_no_gold.png   (gold standard 없음, gray와 red만)

Per panel (gold standard 포함):
  - gray   = 일반 유전자
  - red    = propagation score 상위 TOP_FRAC (예: 5%)
  - BLUE   = gold standard (각 그래프마다 정의 다름)
  - blue solid line  = 전체 trend
  - red dashed line  = top-TOP_FRAC trend

Per panel (gold standard 없음):
  - gray   = 일반 유전자
  - red    = propagation score 상위 TOP_FRAC (예: 5%)
  - blue solid line  = 전체 trend
  - red dashed line  = top-TOP_FRAC trend

Input:
  - result_np_cutoff/np_cutoff_results.csv
  - w3_analysis_code/gold_standard_validation.csv

Output:
  - result_np_cutoff/scatter_prop_vs_emp_pval_ot.png
  - result_np_cutoff/scatter_prop_vs_emp_pval_chembl.png
  - result_np_cutoff/scatter_prop_vs_emp_pval_no_gold.png

사용법:
    python3 w1_analysis_code/scatterplot.py
"""

import csv
import math
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

INPUT_FILE           = "result_np_cutoff/np_cutoff_results.csv"
GOLD_FILE            = "w3_analysis_code/gold_standard_validation.csv"
OUTPUT_OT_FILE       = "result_np_cutoff/scatter_prop_vs_emp_pval_ot.png"
OUTPUT_CHEMBL_FILE   = "result_np_cutoff/scatter_prop_vs_emp_pval_chembl.png"
OUTPUT_NO_GOLD_FILE  = "result_np_cutoff/scatter_prop_vs_emp_pval_no_gold.png"
TOP_FRAC = 0.05
NCOLS    = 5

SCRIPT_DIR  = Path(__file__).resolve().parent
REPO_ROOT   = SCRIPT_DIR.parent
INPUT_PATH  = REPO_ROOT / INPUT_FILE
GOLD_PATH   = REPO_ROOT / GOLD_FILE
OUTPUT_OT_PATH       = REPO_ROOT / OUTPUT_OT_FILE
OUTPUT_CHEMBL_PATH   = REPO_ROOT / OUTPUT_CHEMBL_FILE
OUTPUT_NO_GOLD_PATH  = REPO_ROOT / OUTPUT_NO_GOLD_FILE


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


def load_gold_standard(gold_column):
    """Gold standard 유전자 집합 로드.
       gold_column: 'ot_gold' 또는 'clinical_gold'
       반환: {trait: set of genes}"""
    if not GOLD_PATH.exists():
        print(f"  [WARNING] Gold standard file not found: {GOLD_PATH}")
        return {}

    gold = {}
    with open(GOLD_PATH, "r", newline="") as fh:
        for row in csv.DictReader(fh):
            # OR mode만 사용
            if row.get('combine_mode') != 'OR':
                continue
            if row.get(gold_column, '').strip().upper() == 'TRUE':
                gold.setdefault(row['Trait'], set()).add(row['gene'])
    return gold


def load_traits_data():
    """propagation 결과 로드"""
    traits = {}
    with open(INPUT_PATH, "r", newline="") as fh:
        for row in csv.DictReader(fh):
            t = traits.setdefault(row["Trait"], {
                "page_rank": [], "neg_log10": [], "gene": [], "n_seeds": None
            })
            t["page_rank"].append(float(row["page.rank"]))
            t["neg_log10"].append(-math.log10(float(row["emp_pval"])))
            t["gene"].append(row.get("gene", ""))
            t["n_seeds"] = int(row["n_seeds"])
    return traits


def make_scatter_plot(traits, gold_genes, output_path, gold_label):
    """scatter plot 생성"""
    n_traits = len(traits)
    ncols = min(NCOLS, n_traits)
    nrows = math.ceil(n_traits / ncols)
    fig, axes = plt.subplots(nrows, ncols,
                              figsize=(4 * ncols, 3.6 * nrows), squeeze=False)
    fig.suptitle(
        f"Propagation score vs -log10(empirical p-value)\n"
        f"red = top {TOP_FRAC*100:g}% by score  |  BLUE = {gold_label}",
        fontsize=13,
    )

    for idx, (trait, t) in enumerate(sorted(traits.items())):
        ax       = axes[idx // ncols][idx % ncols]
        x        = t["page_rank"]
        y        = t["neg_log10"]
        genes    = t["gene"]
        n        = len(x)
        n_seeds  = t["n_seeds"]
        top_idx  = set(top_frac_indices(x, TOP_FRAC))

        # gold standard 인덱스
        gold_set     = gold_genes.get(trait, set())
        gold_idx     = [i for i in range(n) if genes[i] in gold_set]
        gold_set_idx = set(gold_idx)
        top_only_idx = [i for i in top_idx if i not in gold_set_idx]
        gray_idx     = [i for i in range(n)
                        if i not in top_idx and i not in gold_set_idx]

        # 레이어 1: 일반 (gray)
        ax.scatter([x[i] for i in gray_idx], [y[i] for i in gray_idx],
                   s=2, alpha=0.2, c="gray", zorder=1)
        # 레이어 2: top 5% (gold 아닌 것)
        ax.scatter([x[i] for i in top_only_idx], [y[i] for i in top_only_idx],
                   s=6, alpha=0.7, c="red", zorder=2)
        # 레이어 3: gold standard
        ax.scatter([x[i] for i in gold_idx], [y[i] for i in gold_idx],
                   s=16, alpha=0.8, c="dodgerblue",
                   edgecolors="navy", linewidths=0.4, zorder=3)

        # trend lines
        if n > 1:
            a, b = np.polyfit(x, y, 1)
            xs = np.linspace(min(x), max(x), 100)
            ax.plot(xs, a * xs + b, "b-", lw=1.2, alpha=0.8, zorder=4)
        tx = [x[i] for i in top_idx]; ty = [y[i] for i in top_idx]
        if len(tx) > 1:
            a, b = np.polyfit(tx, ty, 1)
            xs = np.linspace(min(tx), max(tx), 50)
            ax.plot(xs, a * xs + b, "r--", lw=1.5, alpha=0.9, zorder=4)

        r_all = pearson_correlation(x, y)
        r_top = pearson_correlation(tx, ty)
        n_gold_in_panel = len(gold_idx)
        n_gold_sig      = sum(1 for i in gold_idx if y[i] > -math.log10(0.05))

        clean = trait.replace("finngen_R12_", "").replace("_", " ")
        title = (f"{clean}\nseeds = {n_seeds},  n = {n}\n"
                 f"r(all)={r_all:.3f}  r(top {TOP_FRAC*100:g}%)={r_top:.3f}\n"
                 f"Gold: {n_gold_in_panel}  |  sig: {n_gold_sig}")
        ax.set_title(title, fontsize=8)
        ax.set_xlabel("propagation score (page.rank)", fontsize=8)
        ax.set_ylabel("-log10(emp_pval)", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.grid(True, alpha=0.3, zorder=0)
        ax.axhline(y=-math.log10(0.05), color="orange", linestyle="--",
                   alpha=0.6, linewidth=1, zorder=0.5)
        ax.axhline(y=-math.log10(0.01), color="red", linestyle="--",
                   alpha=0.6, linewidth=1, zorder=0.5)

    for idx in range(n_traits, nrows * ncols):
        axes[idx // ncols][idx % ncols].axis("off")

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    os.makedirs(output_path.parent, exist_ok=True)
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {output_path}")


def make_scatter_plot_no_gold(traits, output_path):
    """scatter plot 생성 (gold standard 없음)"""
    n_traits = len(traits)
    ncols = min(NCOLS, n_traits)
    nrows = math.ceil(n_traits / ncols)
    fig, axes = plt.subplots(nrows, ncols,
                              figsize=(4 * ncols, 3.6 * nrows), squeeze=False)
    fig.suptitle(
        f"Propagation score vs -log10(empirical p-value)\n"
        f"red = top {TOP_FRAC*100:g}% by score",
        fontsize=13,
    )

    for idx, (trait, t) in enumerate(sorted(traits.items())):
        ax       = axes[idx // ncols][idx % ncols]
        x        = t["page_rank"]
        y        = t["neg_log10"]
        n        = len(x)
        n_seeds  = t["n_seeds"]
        top_idx  = set(top_frac_indices(x, TOP_FRAC))

        # 인덱스 분리
        gray_idx = [i for i in range(n) if i not in top_idx]
        red_idx  = list(top_idx)

        # 레이어 1: 일반 (gray)
        ax.scatter([x[i] for i in gray_idx], [y[i] for i in gray_idx],
                   s=2, alpha=0.2, c="gray", zorder=1)
        # 레이어 2: top 5% (red)
        ax.scatter([x[i] for i in red_idx], [y[i] for i in red_idx],
                   s=6, alpha=0.7, c="red", zorder=2)

        # trend lines
        if n > 1:
            a, b = np.polyfit(x, y, 1)
            xs = np.linspace(min(x), max(x), 100)
            ax.plot(xs, a * xs + b, "b-", lw=1.2, alpha=0.8, zorder=4)
        tx = [x[i] for i in top_idx]; ty = [y[i] for i in top_idx]
        if len(tx) > 1:
            a, b = np.polyfit(tx, ty, 1)
            xs = np.linspace(min(tx), max(tx), 50)
            ax.plot(xs, a * xs + b, "r--", lw=1.5, alpha=0.9, zorder=4)

        r_all = pearson_correlation(x, y)
        r_top = pearson_correlation(tx, ty)

        clean = trait.replace("finngen_R12_", "").replace("_", " ")
        title = (f"{clean}\nseeds = {n_seeds},  n = {n}\n"
                 f"r(all)={r_all:.3f}  r(top {TOP_FRAC*100:g}%)={r_top:.3f}")
        ax.set_title(title, fontsize=8)
        ax.set_xlabel("propagation score (page.rank)", fontsize=8)
        ax.set_ylabel("-log10(emp_pval)", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.grid(True, alpha=0.3, zorder=0)
        ax.axhline(y=-math.log10(0.05), color="orange", linestyle="--",
                   alpha=0.6, linewidth=1, zorder=0.5)
        ax.axhline(y=-math.log10(0.01), color="red", linestyle="--",
                   alpha=0.6, linewidth=1, zorder=0.5)

    for idx in range(n_traits, nrows * ncols):
        axes[idx // ncols][idx % ncols].axis("off")

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    os.makedirs(output_path.parent, exist_ok=True)
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {output_path}")


def main():
    traits = load_traits_data()
    print(f"Loaded propagation data: {len(traits)} traits\n")

    # 1. OT score 기준
    ot_gold  = load_gold_standard("ot_gold")
    total_ot = sum(len(s) for s in ot_gold.values())
    print(f"OT score gold:        {total_ot} genes across {len(ot_gold)} traits")
    make_scatter_plot(traits, ot_gold, OUTPUT_OT_PATH,
                       gold_label="OT score gold")

    # 2. ChEMBL/clinical 기준
    chembl_gold = load_gold_standard("clinical_gold")
    total_cp    = sum(len(s) for s in chembl_gold.values())
    print(f"ChEMBL/clinical gold: {total_cp} genes across {len(chembl_gold)} traits")
    make_scatter_plot(traits, chembl_gold, OUTPUT_CHEMBL_PATH,
                       gold_label="ChEMBL gold (clinical_precedence > 0)")

    # 3. Gold standard 없음
    print(f"No gold standard:     (simple gray + red visualization)")
    make_scatter_plot_no_gold(traits, OUTPUT_NO_GOLD_PATH)


if __name__ == "__main__":
    main()