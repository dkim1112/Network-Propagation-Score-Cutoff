#!/usr/bin/env python3
"""
Quadrant Enrichment Test (Hypergeometric)

목적: 각 사분면에 관찰된 gold standard 수가 우연 대비 의미있게 많은지 검정

박사님 제안:
  "좌상단 26개의 gold standard gene이 많은건지 적은건지?
   좌상단 크기와 동일한 수의 유전자를 랜덤하게 뽑고, gold가 얼마나
   포함되는지 1,000회 반복해서 비교"

→ 이건 수학적으로 hypergeometric distribution 그 자체임.
   permutation으로 1,000회 sampling 하나, scipy.stats.hypergeom으로
   해석해를 구하나 결과 동일. hypergeometric이 더 정확하고 빠름.

원리:
  - 전체 풀 크기:    K
  - 풀 안의 gold:    G
  - 사분면 크기:     N
  - 관찰된 gold:     X
  - P(>= X)  = hypergeom.sf(X - 1, K, G, N)

두 가지 분석:
  1. Top row: Pooled (전체 합산) — 모든 trait 합쳐서 한 번의 hypergeometric test
  2. Bottom row: per-trait p-value distribution (box plot + strip overlay)

Input:  w5_analysis_code/quadrant_enrichment.csv
        (Trait, gold_type, quadrant, n_total, n_gold 컬럼 필요)

Output: top_left_permutation_test.csv
        top_left_permutation_test_plot.png

사용법:
    python3 w5_analysis_code/top_left_permutation_test.py
"""
import pandas as pd
import numpy as np
from scipy.stats import hypergeom
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os, sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
INPUT_CSV  = os.path.normpath(
    os.path.join(SCRIPT_DIR, "quadrant_enrichment.csv"))
OUTPUT_DIR = SCRIPT_DIR
ALPHA = 0.05

if len(sys.argv) > 1:
    INPUT_CSV = sys.argv[1]
if not os.path.exists(INPUT_CSV):
    print(f"ERROR: {INPUT_CSV} not found"); sys.exit(1)

print(f"Reading {INPUT_CSV}...")
df = pd.read_csv(INPUT_CSV)
print(f"  {len(df):,} rows | {df['Trait'].nunique()} traits\n")


def hypergeom_test(n_total, n_gold, pool_size, pool_gold):
    if n_total == 0 or pool_size == 0:
        return float('nan'), float('nan')
    expected = n_total * pool_gold / pool_size
    pval = hypergeom.sf(n_gold - 1, pool_size, pool_gold, n_total)
    return float(pval), float(expected)


QUADRANTS = ['UL_low_prop_sig', 'UR_high_prop_sig',
             'LL_low_prop_nonsig', 'LR_high_prop_nonsig']
quad_labels = {
    'UL_low_prop_sig':      'UL\n(low+sig)',
    'UR_high_prop_sig':     'UR\n(high+sig)',
    'LL_low_prop_nonsig':   'LL\n(low+nonsig)',
    'LR_high_prop_nonsig':  'LR\n(high+nonsig)',
}
colors = {'UL_low_prop_sig':     '#ff7f0e',
          'UR_high_prop_sig':    '#2ca02c',
          'LL_low_prop_nonsig':  '#d62728',
          'LR_high_prop_nonsig': '#7f7f7f'}


# == 1. Pooled ===============================================================
pooled_results = []
for gold_type in df['gold_type'].unique():
    sub = df[df['gold_type'] == gold_type]
    pool_size = int(sub['n_total'].sum())
    pool_gold = int(sub['n_gold'].sum())
    for q in QUADRANTS:
        q_sub = sub[sub['quadrant'] == q]
        n_total = int(q_sub['n_total'].sum())
        n_gold  = int(q_sub['n_gold'].sum())
        pval, exp = hypergeom_test(n_total, n_gold, pool_size, pool_gold)
        enrichment = n_gold / exp if exp > 0 else float('nan')
        pooled_results.append({
            'gold_type': gold_type, 'quadrant': q,
            'pool_size': pool_size, 'pool_gold': pool_gold,
            'n_total': n_total, 'observed': n_gold,
            'expected': round(exp, 2), 'enrichment': round(enrichment, 3),
            'p_value': pval, 'significant': pval < ALPHA,
        })
pooled_df = pd.DataFrame(pooled_results)

# == 2. Per-trait ===========================================================
per_trait_results = []
for gold_type in df['gold_type'].unique():
    sub = df[df['gold_type'] == gold_type]
    for trait in sub['Trait'].unique():
        trait_sub = sub[sub['Trait'] == trait]
        pool_size = int(trait_sub['n_total'].sum())
        pool_gold = int(trait_sub['n_gold'].sum())
        if pool_gold == 0:
            continue
        for _, row in trait_sub.iterrows():
            pval, exp = hypergeom_test(
                int(row['n_total']), int(row['n_gold']),
                pool_size, pool_gold)
            enrichment = row['n_gold'] / exp if exp > 0 else float('nan')
            per_trait_results.append({
                'Trait': trait, 'gold_type': gold_type,
                'quadrant': row['quadrant'],
                'pool_size': pool_size, 'pool_gold': pool_gold,
                'n_total': int(row['n_total']), 'observed': int(row['n_gold']),
                'expected': round(exp, 2),
                'enrichment': round(enrichment, 3) if not np.isnan(enrichment) else float('nan'),
                'p_value': pval, 'significant': pval < ALPHA,
            })
per_trait_df = pd.DataFrame(per_trait_results)

pooled_df.to_csv(os.path.join(OUTPUT_DIR, "quadrant_enrichment_test_pooled.csv"), index=False)
per_trait_df.to_csv(os.path.join(OUTPUT_DIR, "quadrant_enrichment_test_per_trait.csv"), index=False)


# == 시각화 ===================================================================
gold_labels = list(pooled_df['gold_type'].unique())

fig, axes = plt.subplots(2, len(gold_labels),
                          figsize=(4.5 * len(gold_labels), 9))

fig.suptitle(
    "Quadrant Gold Enrichment Test (Hypergeometric)\n"
    "Top row: pooled (observed vs expected, enrichment ratio)\n"
    "Bottom row: per-trait p-value distribution (box plot + each trait as a point)",
    fontsize=12,
)

# Row 1: pooled enrichment ratio
for col, gold_lbl in enumerate(gold_labels):
    ax = axes[0][col]
    sub = pooled_df[pooled_df['gold_type'] == gold_lbl]
    sub = sub.set_index('quadrant').reindex(QUADRANTS)

    x = np.arange(len(QUADRANTS))
    ax.bar(x, sub['enrichment'],
           color=[colors[q] for q in QUADRANTS],
           edgecolor='black', linewidth=0.8)
    ax.axhline(1.0, color='black', linestyle='--', alpha=0.5, linewidth=1)

    ymax = max(sub['enrichment'].max() if not sub['enrichment'].isna().all() else 1, 1.1)
    for i, q in enumerate(QUADRANTS):
        enr = sub.loc[q, 'enrichment']
        pval = sub.loc[q, 'p_value']
        obs = int(sub.loc[q, 'observed'])
        exp = sub.loc[q, 'expected']
        if np.isnan(enr): continue
        sig_mark = "*" if pval < ALPHA else ""
        ax.text(i, enr + ymax * 0.02,
                f"{enr:.2f}{sig_mark}\nobs={obs}\nexp={exp:.1f}",
                ha='center', fontsize=7)

    ax.set_xticks(x)
    ax.set_xticklabels([quad_labels[q] for q in QUADRANTS], fontsize=8)
    ax.set_ylabel("Enrichment (obs / expected)")
    ax.set_title(f"Pooled: {gold_lbl}")
    ax.grid(True, alpha=0.3, axis='y')
    ax.set_ylim(0, ymax * 1.30)

# Row 2: box plot + strip overlay
rng = np.random.default_rng(42)
for col, gold_lbl in enumerate(gold_labels):
    ax = axes[1][col]
    sub = per_trait_df[per_trait_df['gold_type'] == gold_lbl]

    # 각 quadrant의 p-value array 수집
    box_data = []
    valid_positions = []
    for i, q in enumerate(QUADRANTS):
        q_data = sub[sub['quadrant'] == q]
        pvals = q_data['p_value'].values
        box_data.append(pvals)
        if len(pvals) > 0:
            valid_positions.append(i)

    # box plot 그리기 (개별 color, 반투명)
    positions = np.arange(len(QUADRANTS))
    bp = ax.boxplot(
        box_data,
        positions=positions,
        widths=0.55,
        patch_artist=True,
        showfliers=False,        # outlier는 strip에서 별도로 보여줌
        medianprops=dict(color='black', linewidth=2),
        whiskerprops=dict(color='black', linewidth=1.2),
        capprops=dict(color='black', linewidth=1.2),
        boxprops=dict(linewidth=1.2),
        zorder=2,
    )
    for patch, q in zip(bp['boxes'], QUADRANTS):
        patch.set_facecolor(colors[q])
        patch.set_alpha(0.35)

    # strip overlay (jittered points)
    for i, q in enumerate(QUADRANTS):
        pvals = box_data[i]
        if len(pvals) == 0:
            continue
        jitter = rng.uniform(-0.15, 0.15, size=len(pvals))
        ax.scatter([i + j for j in jitter], pvals,
                   color=colors[q], alpha=0.85, s=45,
                   edgecolors='black', linewidth=0.6, zorder=4)

    ax.axhline(0.05, color='red', linestyle='--', alpha=0.8,
               linewidth=1.3, label='p = 0.05', zorder=3)
    ax.set_xticks(positions)
    ax.set_xticklabels([quad_labels[q] for q in QUADRANTS], fontsize=8)
    ax.set_ylabel("Per-trait hypergeometric p-value")
    ax.set_title(f"Per-trait p-values: {gold_lbl}")
    ax.set_ylim(-0.05, 1.05)
    ax.grid(True, alpha=0.3, axis='y')
    ax.legend(fontsize=8, loc='center right')

plt.tight_layout()
out_png = os.path.join(OUTPUT_DIR, "quadrant_enrichment_test_plot.png")
plt.savefig(out_png, dpi=140, bbox_inches='tight')
plt.close()
print(f"Saved: {out_png}")