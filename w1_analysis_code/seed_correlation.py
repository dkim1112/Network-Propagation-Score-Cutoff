#!/usr/bin/env python3
"""
Seed correlation pattern comparison.

Research Question: Does the reliability of permutation-based empirical p-values
decrease when seed count is low due to increased variance in null distributions?

Hypothesis: With fewer seeds (1-2), individual random gene network positions have
larger impact on null distribution → higher variance → less reliable emp_pval
→ weaker correlation between propagation score and statistical significance.

Input: np_cutoff_results CSV file (emp_pval, page.rank, n_seeds columns)
Output: - seed_vs_correlation.png (scatter plot)
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
OUTPUT_PNG = "result_np_cutoff/seed_vs_correlation.png"
# =============================================================================

# 커맨드라인 인수로 파일명 받을 수 있도록
if len(sys.argv) > 1:
    INPUT_CSV = sys.argv[1]

if not os.path.exists(INPUT_CSV):
    print(f"ERROR: {INPUT_CSV} 파일을 찾을 수 없습니다.")
    print("사용법: python seed_correlation_analysis.py [csv파일명]")
    sys.exit(1)

df = pd.read_csv(INPUT_CSV)
df['log10_emp'] = -np.log10(df['emp_pval'].clip(lower=1e-10))

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

summary = pd.DataFrame(rows).sort_values('n_seeds')

if len(summary) >= 3:
    r_seed_vs_corr, p_seed = stats.pearsonr(summary['n_seeds'], summary['r_all_emp'])
else:
    r_seed_vs_corr, p_seed = 0, 1

# == Scatter Plot ==============================================================
seed_counts = summary['n_seeds'].values
r_all_vals  = summary['r_all_emp'].values

fig, ax = plt.subplots(1, 1, figsize=(8, 6))
fig.suptitle('Impact of Seed Count on Empirical P-value Reliability', fontsize=14)

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
ax.set_xlabel('Number of seed genes', fontsize=11)
ax.set_ylabel('Correlation: Propagation Score vs -log10(Empirical P-value)', fontsize=11)
ax.set_title(f'Fewer seeds → Less reliable empirical p-values?\n(r={r_seed_vs_corr:.3f}, p={p_seed:.3f})', fontsize=11)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(OUTPUT_PNG, dpi=150, bbox_inches='tight')