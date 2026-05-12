#!/usr/bin/env python3
"""
Seed correlation pattern comparison.
Input: np_cutoff_results CSV file (emp_pval, page.rank, n_seeds columns)
Output: - seed_correlation_summary.csv (trait-level summary)
      - seed_vs_correlation.png (scatter plot)
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
import sys
import os

# == 파라미터 ==================================================================
INPUT_CSV  = "result_np_cutoff/np_cutoff_results.csv"   # 기본 입력 파일
OUTPUT_CSV = "seed_correlation_summary.csv"
OUTPUT_PNG = "seed_vs_correlation.png"
# =============================================================================

# 커맨드라인 인수로 파일명 받을 수 있도록
if len(sys.argv) > 1:
    INPUT_CSV = sys.argv[1]

if not os.path.exists(INPUT_CSV):
    print(f"ERROR: {INPUT_CSV} 파일을 찾을 수 없습니다.")
    print("사용법: python seed_correlation_analysis.py [csv파일명]")
    sys.exit(1)

print(f"Reading {INPUT_CSV}...")
df = pd.read_csv(INPUT_CSV)
df['log10_emp'] = -np.log10(df['emp_pval'].clip(lower=1e-10))

print(f"  총 {len(df):,}행 | {df['Trait'].nunique()}개 질환\n")

# == Correlation by trait ======================================================
rows = []
for trait, sub in df.groupby('Trait'):
    n_seeds    = sub['n_seeds'].iloc[0]
    page_ranks = sub['page.rank'].values
    log_emp    = sub['log10_emp'].values

    # Overall Pearson r
    r_all, _   = stats.pearsonr(page_ranks, log_emp)

    row = {
        'Trait':     trait,
        'n_seeds':   n_seeds,
        'n_genes':   len(sub),
        'r_all_emp': round(r_all, 4),
    }

    rows.append(row)
    print(f"  {trait.replace('finngen_R12_',''):35s} "
          f"seed={n_seeds:2d} | "
          f"r_all={r_all:.3f}")

summary = pd.DataFrame(rows).sort_values('n_seeds')
summary.to_csv(OUTPUT_CSV, index=False)
print(f"\n결과 저장: {OUTPUT_CSV}")

# == 요약 통계 =================================================================
print("\n=== SUMMARY ===")
print(f"r_all (emp_pval)  — mean: {summary['r_all_emp'].mean():.3f}, "
      f"median: {summary['r_all_emp'].median():.3f}, "
      f"range: {summary['r_all_emp'].min():.3f}~{summary['r_all_emp'].max():.3f}")

# seed 수에 따른 패턴 확인
if len(summary) >= 3:
    r_seed_vs_corr, p_seed = stats.pearsonr(summary['n_seeds'], summary['r_all_emp'])
    print(f"\nseed count vs r_all correlation: r={r_seed_vs_corr:.3f} (p={p_seed:.3f})")
    if r_seed_vs_corr > 0.3:
        print("  -> higher seed counts tend to have higher correlation")
    else:
        print("  -> no clear trend between seed count and correlation")

# == Scatter Plot ==============================================================
seed_counts = summary['n_seeds'].values
r_all_vals  = summary['r_all_emp'].values

fig, ax = plt.subplots(1, 1, figsize=(6, 5))
fig.suptitle('Seed Count vs Correlation Pattern', fontsize=14)

ax.scatter(seed_counts, r_all_vals, s=80, color='steelblue', alpha=0.8)
for i, row in summary.iterrows():
    ax.annotate(row['Trait'].replace('finngen_R12_',''),
                (row['n_seeds'], row['r_all_emp']),
                fontsize=7, ha='left', va='bottom',
                xytext=(3, 3), textcoords='offset points')
if len(summary) >= 3:
    z = np.polyfit(seed_counts, r_all_vals, 1)
    x_line = np.linspace(seed_counts.min(), seed_counts.max(), 100)
    ax.plot(x_line, np.poly1d(z)(x_line), 'r--', alpha=0.6)
ax.set_xlabel('Seed count (n_seeds)', fontsize=11)
ax.set_ylabel('Pearson r (overall)', fontsize=11)
ax.set_title(f'Seed count vs r_all\n(emp_pval, r={r_seed_vs_corr:.3f})', fontsize=11)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(OUTPUT_PNG, dpi=150, bbox_inches='tight')
print(f"Plot saved: {OUTPUT_PNG}")