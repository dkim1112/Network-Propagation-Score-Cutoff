#!/usr/bin/env python3
"""
Top-N Propagation Score vs Empirical P-value 비교 분석

박사님 요청: "network propagation score를 사용했을 때랑 False Negative를 비교해보면 좋겠다"

비교: Gold standard를 두 방법이 각각 어떻게 잡았는지 4가지 카테고리로 분할
  - both_caught:      둘 다 잡음
  - only_emp_caught:  empirical p-value만 잡음
  - only_topn_caught: top-N만 잡음
  - both_missed:      둘 다 놓침

Top-N 정의: N = 각 질환의 seed 수 (박사님 지정)
            page.rank 상위 N개 유전자

Input:  gold_standard_validation.csv  (gold_standard.py 결과)
Output: topn_vs_emp_comparison.csv, topn_vs_emp_fn_comparison.csv,
        topn_vs_emp_overlap.csv, topn_vs_emp_plot.png, topn_vs_emp_fn_plot.png

사용법:
    python3 w3_analysis_code/topn_vs_emp_comparison.py [gold_standard_validation.csv]
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys, os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
INPUT_CSV  = os.path.join(SCRIPT_DIR, "gold_standard_validation.csv")
OUTPUT_DIR = SCRIPT_DIR

if len(sys.argv) > 1:
    INPUT_CSV = sys.argv[1]
if not os.path.exists(INPUT_CSV):
    sys.exit(1)

detail = pd.read_csv(INPUT_CSV)
detail = detail[detail['combine_mode'] == 'OR'].copy()

results, fn_rows, overlap_rows = [], [], []

for trait in sorted(detail['Trait'].unique()):
    sub = detail[detail['Trait'] == trait].copy()
    n_seeds = int(sub['n_seeds'].iloc[0])
    short   = trait.replace("finngen_R12_", "")

    # 두 방법의 결과 set (Fair: 둘 다 M개)
    emp_genes  = set(sub[sub['significant']]['gene'])              # M = significant 전체
    M          = len(emp_genes)
    topn_genes = set(sub.nlargest(M, 'page.rank')['gene']) if M > 0 else set()
    gold_genes = set(sub[sub['gold']]['gene'])
    


    # Gold standard 4 카테고리 분할
    both_caught   = gold_genes & topn_genes & emp_genes
    only_emp      = (gold_genes & emp_genes) - topn_genes
    only_topn     = (gold_genes & topn_genes) - emp_genes
    both_missed   = gold_genes - topn_genes - emp_genes

    topn_fn = len(gold_genes - topn_genes)
    emp_fn  = len(gold_genes - emp_genes)

    n_gold = len(gold_genes)

    results.append({
        'Trait': short, 'n_seeds': n_seeds, 'n_gold': n_gold,
        'topn_FN': topn_fn,
        'emp_FN': emp_fn,
        'both_caught':      len(both_caught),
        'only_emp_caught':  len(only_emp),
        'only_topn_caught': len(only_topn),
        'both_missed':      len(both_missed),
        'both_caught_pct':      round(len(both_caught) / n_gold * 100, 2) if n_gold > 0 else 0,
        'only_emp_pct':         round(len(only_emp)    / n_gold * 100, 2) if n_gold > 0 else 0,
        'only_topn_pct':        round(len(only_topn)   / n_gold * 100, 2) if n_gold > 0 else 0,
        'both_missed_pct':      round(len(both_missed) / n_gold * 100, 2) if n_gold > 0 else 0,
    })

    fn_rows.append({
        'Trait': short,
        'n_seeds': n_seeds,
        'n_gold': n_gold,
        'topn_FN': topn_fn,
        'emp_FN': emp_fn,
        'topn_FN_pct': round(topn_fn / n_gold * 100, 2) if n_gold > 0 else 0,
        'emp_FN_pct': round(emp_fn / n_gold * 100, 2) if n_gold > 0 else 0,
    })

    for gene in sorted(only_emp):
        overlap_rows.append({'Trait': short, 'gene': gene, 'category': 'only_emp_caught'})
    for gene in sorted(only_topn):
        overlap_rows.append({'Trait': short, 'gene': gene, 'category': 'only_topn_caught'})
    for gene in sorted(both_caught):
        overlap_rows.append({'Trait': short, 'gene': gene, 'category': 'both_caught'})
    for gene in sorted(both_missed):
        overlap_rows.append({'Trait': short, 'gene': gene, 'category': 'both_missed'})

# == 저장 =====================================================================
res_df     = pd.DataFrame(results)
fn_df      = pd.DataFrame(fn_rows)
overlap_df = pd.DataFrame(overlap_rows)
res_df.to_csv(    os.path.join(OUTPUT_DIR, "topn_vs_emp_comparison.csv"), index=False)
fn_df.to_csv(     os.path.join(OUTPUT_DIR, "topn_vs_emp_fn_comparison.csv"), index=False)
overlap_df.to_csv(os.path.join(OUTPUT_DIR, "topn_vs_emp_overlap.csv"),    index=False)

# == 플롯 ======================================================================
fig, axes = plt.subplots(1, 2, figsize=(15, 6))
fig.suptitle("Top-N Propagation vs Empirical P-value: Gold Standard", fontsize=14)

x = np.arange(len(res_df))
labels = res_df['Trait'].str[:20]

# (a) 절대 개수 (stacked)
ax = axes[0]
bottom = np.zeros(len(res_df))
for col, color, label in [
    ('both_caught',      '#2ca02c', 'Both caught'),
    ('only_emp_caught',  '#ff7f0e', 'Only Emp caught'),
    ('only_topn_caught', '#1f77b4', 'Only TopN caught'),
    ('both_missed',      '#d62728', 'Both missed'),
]:
    ax.bar(x, res_df[col], bottom=bottom, label=label, color=color)
    bottom += res_df[col].values
ax.set_ylabel("Gold standard count")
ax.set_title("(a) Absolute Counts")
ax.set_xticks(x); ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=8)
ax.legend(fontsize=9); ax.grid(True, alpha=0.3, axis='y')

# (b) 비율 (%)
ax = axes[1]
bottom = np.zeros(len(res_df))
for col, color, label in [
    ('both_caught_pct',  '#2ca02c', 'Both caught'),
    ('only_emp_pct',     '#ff7f0e', 'Only Emp caught'),
    ('only_topn_pct',    '#1f77b4', 'Only TopN caught'),
    ('both_missed_pct',  '#d62728', 'Both missed'),
]:
    ax.bar(x, res_df[col], bottom=bottom, label=label, color=color)
    bottom += res_df[col].values
ax.set_ylabel("Gold standard ratio (%)")
ax.set_title("(b) Ratio")
ax.set_xticks(x); ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=8)
ax.set_ylim(0, 100)
ax.legend(fontsize=9); ax.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "topn_vs_emp_plot.png"), dpi=150, bbox_inches='tight')
plt.close()

# == FN 전용 플롯 ==============================================================
fig, ax = plt.subplots(figsize=(15, 6))
fig.suptitle("Top-N vs Empirical P-value: False Negative Counts by Disease", fontsize=14)

x = np.arange(len(fn_df))
labels = fn_df['Trait'].str[:20]
width = 0.38
ax.bar(x - width/2, fn_df['topn_FN'], width, label='Top-N FN', color='steelblue')
ax.bar(x + width/2, fn_df['emp_FN'],  width, label='Emp p-value FN', color='tomato')
ax.set_ylabel("False Negative count")
ax.set_title("False Negatives by Disease")
ax.set_xticks(x)
ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=8)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "topn_vs_emp_fn_plot.png"), dpi=150, bbox_inches='tight')
plt.close()