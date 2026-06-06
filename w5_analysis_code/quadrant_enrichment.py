#!/usr/bin/env python3
"""
Quadrant Gold Standard Distribution Analysis

목적: emp_pval 사용 정당성 검증
  - 좌상단(low propagation, significant emp_pval)에 gold standard가 많이 분포해야 함
  - 좌하단(propagation도 낮고 emp_pval도 유의하지 않음)에는 적게 분포해야 함

사분면 정의:
  - UL: propagation < lower-quantile  AND  emp_pval <  0.05
  - UR: propagation >= lower-quantile AND  emp_pval <  0.05
  - LL: propagation < lower-quantile  AND  emp_pval >= 0.05  (Negative control)
  - LR: propagation >= lower-quantile AND  emp_pval >= 0.05

Gold standard (2가지 비교):
  - OT globalScore >= cutoff (0.3 / 0.5 / 0.75)
  - ChEMBL: clinical_precedence > 0 (binary)

Input:  w3_analysis_code/gold_standard_validation.csv
Output: w5_analysis_code/quadrant_enrichment.csv
        w5_analysis_code/quadrant_enrichment_plot.png

"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os, sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# =============================================================================
# ★ 조정 가능 ★
# =============================================================================
INPUT_CSV  = os.path.normpath(
    os.path.join(SCRIPT_DIR, "..", "w3_analysis_code", "gold_standard_validation.csv"))
OUTPUT_DIR = SCRIPT_DIR

OT_CUTOFFS         = [0.3, 0.5, 0.75]   # OT globalScore cutoff
PROP_QUANTILE      = 0.25               # 좌상단을 더 좁게: 하위 25%만 "low propagation"
                                         # (이전: 0.5 = median 기준 → 너무 넓음)
EMP_PVAL_THRESHOLD = 0.05
# =============================================================================

if len(sys.argv) > 1:
    INPUT_CSV = sys.argv[1]
if not os.path.exists(INPUT_CSV):
    print(f"ERROR: {INPUT_CSV} not found")
    sys.exit(1)

print(f"Reading {INPUT_CSV}...")
df = pd.read_csv(INPUT_CSV)
df = df[df['combine_mode'] == 'OR'].copy()
print(f"  {len(df):,} rows | {df['Trait'].nunique()} diseases\n")
print(f"Parameters:")
print(f"  OT cutoffs:         {OT_CUTOFFS}")
print(f"  PROP quantile:      {PROP_QUANTILE} (lower {PROP_QUANTILE*100:.0f}% = low propagation)")
print(f"  EMP_PVAL threshold: {EMP_PVAL_THRESHOLD}\n")


def assign_quadrant(sub, prop_threshold, pval_threshold):
    low_prop = sub['page.rank'] < prop_threshold
    sig      = sub['emp_pval']  < pval_threshold
    sub = sub.copy()
    sub['quadrant'] = np.select(
        [
            low_prop  &  sig,
            ~low_prop &  sig,
            low_prop  & ~sig,
            ~low_prop & ~sig,
        ],
        ['UL_low_prop_sig', 'UR_high_prop_sig',
         'LL_low_prop_nonsig', 'LR_high_prop_nonsig'],
        default=''
    )
    return sub


def quadrant_stats(sub, quadrant_name):
    in_q = sub['quadrant'] == quadrant_name
    n_total = int(in_q.sum())
    n_gold  = int((in_q & sub['gold']).sum())
    return {
        'n_total':  n_total,
        'n_gold':   n_gold,
        'gold_pct': round(n_gold / n_total * 100, 3) if n_total > 0 else 0,
    }


def compute_for_gold(df_in, gold_label):
    """주어진 gold 기준으로 사분면별 통계 계산
       반환: per-trait records, summary dict for plotting"""
    records = []
    for trait in sorted(df_in['Trait'].unique()):
        sub = df_in[df_in['Trait'] == trait].copy()
        if sub['gold'].sum() == 0:
            continue
        prop_thresh = sub['page.rank'].quantile(PROP_QUANTILE)
        sub = assign_quadrant(sub, prop_thresh, EMP_PVAL_THRESHOLD)
        for q in QUADRANTS:
            r = quadrant_stats(sub, q)
            records.append({
                'Trait':         trait.replace('finngen_R12_', ''),
                'gold_type':     gold_label,
                'quadrant':      q,
                'n_total':       r['n_total'],
                'n_gold':        r['n_gold'],
                'gold_pct':      r['gold_pct'],
            })
    return records


# == 메인 분석 ================================================================
results = []
QUADRANTS = ['UL_low_prop_sig', 'UR_high_prop_sig',
             'LL_low_prop_nonsig', 'LR_high_prop_nonsig']

# 1. OT cutoffs
for ot_cut in OT_CUTOFFS:
    print(f"{'='*60}")
    print(f"OT globalScore >= {ot_cut}")
    df_cut = df.copy()
    df_cut['gold'] = df_cut['ot_score'].fillna(0) >= ot_cut
    label = f"OT>={ot_cut}"
    recs = compute_for_gold(df_cut, label)
    results.extend(recs)

    # 합산 출력
    df_all = df_cut.copy()
    df_all['quadrant'] = ''
    for trait in df_all['Trait'].unique():
        mask = df_all['Trait'] == trait
        prop_thresh = df_all.loc[mask, 'page.rank'].quantile(PROP_QUANTILE)
        df_all.loc[mask] = assign_quadrant(df_all[mask].copy(),
                                            prop_thresh, EMP_PVAL_THRESHOLD).values
    print(f"  {'quadrant':30s} {'n_total':>9s} {'n_gold':>7s} {'gold%':>7s}")
    for q in QUADRANTS:
        r = quadrant_stats(df_all, q)
        print(f"  {q:30s} {r['n_total']:>9d} {r['n_gold']:>7d} {r['gold_pct']:>6.2f}%")
    print()

# 2. ChEMBL (clinical_precedence > 0)
print(f"{'='*60}")
print(f"ChEMBL (clinical_precedence > 0)")
df_cm = df.copy()
df_cm['gold'] = df_cm['clinical_gold'] == True
label = "ChEMBL"
recs = compute_for_gold(df_cm, label)
results.extend(recs)

df_all = df_cm.copy()
df_all['quadrant'] = ''
for trait in df_all['Trait'].unique():
    mask = df_all['Trait'] == trait
    prop_thresh = df_all.loc[mask, 'page.rank'].quantile(PROP_QUANTILE)
    df_all.loc[mask] = assign_quadrant(df_all[mask].copy(),
                                        prop_thresh, EMP_PVAL_THRESHOLD).values
print(f"  {'quadrant':30s} {'n_total':>9s} {'n_gold':>7s} {'gold%':>7s}")
for q in QUADRANTS:
    r = quadrant_stats(df_all, q)
    print(f"  {q:30s} {r['n_total']:>9d} {r['n_gold']:>7d} {r['gold_pct']:>6.2f}%")

# == CSV 저장 =================================================================
res_df = pd.DataFrame(results)
out_csv = os.path.join(OUTPUT_DIR, "quadrant_enrichment.csv")
res_df.to_csv(out_csv, index=False)
print(f"\nSaved: {out_csv}")

# == 시각화 (English only, auto-scaled y-axis) ================================
gold_labels_order = [f"OT>={c}" for c in OT_CUTOFFS] + ["ChEMBL"]
n_panels = len(gold_labels_order)
fig, axes = plt.subplots(1, n_panels,
                          figsize=(5 * n_panels, 6),
                          sharey=False)   # ★ 각 panel 자동 scale (차이 보이게)
if n_panels == 1:
    axes = [axes]

fig.suptitle(
    f"Quadrant Gold Standard Distribution\n"
    f"Expected: UL > LL (UL enriched, LL is negative control)\n"
    f"low_prop = bottom {PROP_QUANTILE*100:.0f}% by propagation score",
    fontsize=11,
)

quad_labels = {
    'UL_low_prop_sig':      'UL\n(low prop\n+ sig)',
    'UR_high_prop_sig':     'UR\n(high prop\n+ sig)',
    'LL_low_prop_nonsig':   'LL\n(low prop\n+ non-sig)',
    'LR_high_prop_nonsig':  'LR\n(high prop\n+ non-sig)',
}
colors = {'UL_low_prop_sig':     '#ff7f0e',
          'UR_high_prop_sig':    '#2ca02c',
          'LL_low_prop_nonsig':  '#d62728',
          'LR_high_prop_nonsig': '#7f7f7f'}

for ax, gold_lbl in zip(axes, gold_labels_order):
    sub = res_df[res_df['gold_type'] == gold_lbl]
    stats = sub.groupby('quadrant').agg(
        gold_pct_mean = ('gold_pct', 'mean'),
        n_traits      = ('Trait', 'count'),
        n_gold_sum    = ('n_gold', 'sum'),
    ).reindex(QUADRANTS)

    x = np.arange(len(QUADRANTS))
    bars = ax.bar(x, stats['gold_pct_mean'],
                   color=[colors[q] for q in QUADRANTS],
                   edgecolor='black', linewidth=0.8)

    # value label on top of each bar
    ymax = stats['gold_pct_mean'].max() if stats['gold_pct_mean'].max() > 0 else 1
    for i, q in enumerate(QUADRANTS):
        pct = stats.loc[q, 'gold_pct_mean']
        n_g = int(stats.loc[q, 'n_gold_sum'])
        ax.text(i, pct + ymax * 0.02,
                f"{pct:.2f}%\nn={n_g}",
                ha='center', fontsize=8)

    ax.set_xticks(x)
    ax.set_xticklabels([quad_labels[q] for q in QUADRANTS], fontsize=8)
    ax.set_ylabel("Gold standard ratio (trait mean, %)")
    ax.set_title(gold_lbl)
    ax.grid(True, alpha=0.3, axis='y')
    # add headroom for text labels
    ax.set_ylim(0, ymax * 1.25)

plt.tight_layout()
out_png = os.path.join(OUTPUT_DIR, "quadrant_enrichment_plot.png")
plt.savefig(out_png, dpi=150, bbox_inches='tight')
plt.close()
print(f"Saved: {out_png}")