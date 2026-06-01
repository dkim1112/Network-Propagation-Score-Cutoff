#!/usr/bin/env python3
"""
P-value Precision Analysis

목적: 1K vs 10K permutation에서 LOO empirical p-value 정밀도 비교
핵심 질문: 더 많은 permutation이 p-value를 더 작게 만드는가? 얼마나?
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path
from scipy.stats import pearsonr

# 설정
BASE_DIR = Path(__file__).parent
TEST_DATA_DIR = BASE_DIR / "test_data"
OUTPUT_DIR = BASE_DIR / "precision_analysis_results"
OUTPUT_DIR.mkdir(exist_ok=True)

def main():
    # 데이터 로드
    df_1k = pd.read_csv(TEST_DATA_DIR / "loo_results_1000.csv")   # 1K permutations
    df_10k = pd.read_csv(TEST_DATA_DIR / "loo_results_10000.csv")  # 10K permutations

    # 공통 케이스 매칭
    df_1k['case_id'] = df_1k['Trait'] + '_' + df_1k['left_out_ENSG']
    df_10k['case_id'] = df_10k['Trait'] + '_' + df_10k['left_out_ENSG']

    common_cases = set(df_1k['case_id']) & set(df_10k['case_id'])
    df_1k_common = df_1k[df_1k['case_id'].isin(common_cases)].set_index('case_id')
    df_10k_common = df_10k[df_10k['case_id'].isin(common_cases)].set_index('case_id')

    # 비교 데이터 생성
    comparison_data = []
    for case in common_cases:
        row_1k = df_1k_common.loc[case]
        row_10k = df_10k_common.loc[case]

        comparison_data.append({
            'case_id': case,
            'trait': row_1k['Trait'].replace('finngen_R12_', ''),
            'gene': row_1k['left_out_gene'],
            'pval_1k': row_1k['emp_pval_loo'],
            'pval_10k': row_10k['emp_pval_loo'],
            'recovered_1k': row_1k['recovered'],
            'recovered_10k': row_10k['recovered']
        })

    df = pd.DataFrame(comparison_data)

    # 핵심 통계 - 방향성 중심
    pearson_r, _ = pearsonr(df['pval_1k'], df['pval_10k'])
    pval_diff = df['pval_1k'] - df['pval_10k']  # 양수 = 10K에서 더 작은 p-value (개선)
    df['pval_diff'] = pval_diff  # Add to dataframe for later use
    recovery_changed = (df['recovered_1k'] != df['recovered_10k']).sum()

    # 방향별 분석
    smaller_in_10k = pval_diff > 0
    larger_in_10k = pval_diff < 0
    same_in_10k = pval_diff == 0

    n_smaller = smaller_in_10k.sum()
    n_larger = larger_in_10k.sum()
    n_same = same_in_10k.sum()

    mean_decrease = pval_diff[smaller_in_10k].mean() if n_smaller > 0 else 0
    mean_increase = pval_diff[larger_in_10k].mean() if n_larger > 0 else 0

    # Removed print statements - all info will be in the plot

    # 방향성 변화에 집중한 시각화 (1x2)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

    # 1. P-value 변화: 1K vs 10K 직접 비교 with case identification
    colors = ['red' if diff > 0 else 'blue' for diff in pval_diff]

    ax1.scatter(df['pval_1k'], df['pval_10k'], c=colors, s=60, alpha=0.7, edgecolors='black', linewidth=0.5)
    ax1.plot([0, 1], [0, 1], 'k--', alpha=0.5, linewidth=1, label='No change line')
    ax1.axhline(0.05, color='gray', linestyle=':', alpha=0.7, label='p=0.05 threshold')
    ax1.axvline(0.05, color='gray', linestyle=':', alpha=0.7)

    ax1.set_xlabel('P-value (1K permutations)')
    ax1.set_ylabel('P-value (10K permutations)')
    ax1.set_title(f'Paired P-value Comparison: 1K → 10K\n({len(df)} matched LOO cases)')
    ax1.grid(True, alpha=0.3)

    # 범례 추가
    red_patch = mpatches.Patch(color='red', label='10K < 1K (improved)')
    blue_patch = mpatches.Patch(color='blue', label='10K > 1K (worse)')
    ax1.legend(handles=[red_patch, blue_patch], loc='upper right')

    # 의미있는 범위로 확대
    max_pval = max(df['pval_1k'].max(), df['pval_10k'].max())
    if max_pval > 0.2:
        ax1.set_xlim(0, 0.2)
        ax1.set_ylim(0, 0.2)

    # Case identification - 몇 개 케이스에 라벨 표시
    extreme_cases = df.nlargest(3, 'pval_diff')[['pval_1k', 'pval_10k', 'trait', 'gene']]
    for _, case in extreme_cases.iterrows():
        ax1.annotate(f"{case['trait'][:8]}\n{case['gene'][:8]}",
                    (case['pval_1k'], case['pval_10k']),
                    xytext=(5, 5), textcoords='offset points',
                    fontsize=8, alpha=0.8,
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.7))

    # 2. 통계 요약을 시각적으로 표시
    categories = ['10K Smaller\n(Improved)', 'No Change', '10K Larger\n(Worse)']
    counts = [n_smaller, n_same, n_larger]
    colors_bar = ['red', 'gray', 'blue']

    bars = ax2.bar(categories, counts, color=colors_bar, alpha=0.7, edgecolor='black')
    ax2.set_ylabel('Number of Cases')
    ax2.set_title('Permutation Effect Summary\n(1K → 10K permutations)')
    ax2.grid(True, alpha=0.3, axis='y')

    # 막대 위에 숫자와 평균 변화량 표시
    for i, (bar, count) in enumerate(zip(bars, counts)):
        height = bar.get_height()
        percentage = 100 * count / len(df)
        if i == 0 and n_smaller > 0:
            ax2.text(bar.get_x() + bar.get_width()/2., height + max(counts)*0.01,
                    f'{count} ({percentage:.1f}%)\nAvg Δ: {mean_decrease:.5f}',
                    ha='center', va='bottom', fontweight='bold', fontsize=10)
        elif i == 2 and n_larger > 0:
            ax2.text(bar.get_x() + bar.get_width()/2., height + max(counts)*0.01,
                    f'{count} ({percentage:.1f}%)\nAvg Δ: {mean_increase:.5f}',
                    ha='center', va='bottom', fontweight='bold', fontsize=10)
        else:
            ax2.text(bar.get_x() + bar.get_width()/2., height + max(counts)*0.01,
                    f'{count} ({percentage:.1f}%)',
                    ha='center', va='bottom', fontweight='bold', fontsize=10)

    # 핵심 통계만 텍스트로 추가
    stats_text = f"""Dataset Details:
• Total matched cases: {len(df)}
• Max |change|: {np.abs(pval_diff).max():.6f}
• Recovery status changes: {recovery_changed} cases"""

    ax2.text(0.02, 0.02, stats_text, transform=ax2.transAxes,
            fontsize=9, verticalalignment='bottom',
            bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.8))

    plt.tight_layout()

    # 저장
    plot_path = OUTPUT_DIR / 'pvalue_precision_comparison.png'
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    # 저장 (silent)
    results_path = OUTPUT_DIR / 'pvalue_precision_results.csv'
    df.to_csv(results_path, index=False)

    summary = {
        'metric': ['total_cases', 'correlation', 'mean_abs_diff', 'max_abs_diff', 'recovery_changes', 'change_rate_percent'],
        'value': [len(df), pearson_r, np.abs(pval_diff).mean(), np.abs(pval_diff).max(), recovery_changed, 100*recovery_changed/len(df)]
    }
    summary_path = OUTPUT_DIR / 'pvalue_precision_summary.csv'
    pd.DataFrame(summary).to_csv(summary_path, index=False)

    plt.show()

if __name__ == "__main__":
    main()