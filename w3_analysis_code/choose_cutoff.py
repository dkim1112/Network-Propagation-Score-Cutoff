#!/usr/bin/env python3
"""
Gold Standard Cutoff Analysis
목적:
  1. OpenTargets globalScore 분포 histogram → 적합한 cutoff 탐색
  2. clinical_precedence score 분포 histogram
  3. 다양한 cutoff 조합에서 sensitivity/precision 비교

Output (w3_analysis_code/ 폴더 내):
  - choose_cutoff_analysis_ot_histogram.png
  - choose_cutoff_analysis_cp_histogram.png
  - choose_cutoff_analysis_comparison.csv
  - choose_cutoff_analysis_comparison.png

사용법:
    python3 w3_analysis_code/choose_cutoff.py [np_cutoff_results.csv]
"""

import pandas as pd
import numpy as np
import requests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys, os, time

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# =============================================================================
INPUT_CSV = os.path.normpath(
    os.path.join(SCRIPT_DIR, "..", "result_np_cutoff", "np_cutoff_results.csv"))
OUTPUT_DIR = SCRIPT_DIR

# 시험할 cutoff 값들
OT_CUTOFFS = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
CP_CUTOFFS = [0.0, 0.25, 0.5, 0.75]  # 0.0 = clinical_precedence 아예 안 씀
# =============================================================================

if len(sys.argv) > 1:
    INPUT_CSV = sys.argv[1]
if not os.path.exists(INPUT_CSV):
    print(f"ERROR: {INPUT_CSV} not found"); sys.exit(1)

TRAIT_TO_OT_ID = {
    "finngen_R12_AUTOIMMUNE":           "EFO_0005140",
    "finngen_R12_I9_CHD":               "EFO_0001645",
    "finngen_R12_I9_HYPTENS":           "EFO_0000537",
    "finngen_R12_K11_IBD_STRICT":       "EFO_0003767",
    "finngen_R12_L12_ATOPIC":           "EFO_0000274",
    "finngen_R12_T1D":                  "MONDO_0005147",
    "finngen_R12_T2D_WIDE":             "MONDO_0005148",
    "finngen_R12_CARDIAC_ARRHYTM":      "EFO_0004269",
    "finngen_R12_C3_PROSTATE_EXALLC":   "MONDO_0008315",
    "finngen_R12_C3_BASAL_CELL_CARCINOMA_EXALLC": "EFO_0004193",
    "finngen_R12_C3_SKIN_EXALLC":       "MONDO_0002898",
    "finngen_R12_ALLERG_ASTHMA":        "MONDO_0004784",
    "finngen_R12_ASTHMA_CHILD_EXMORE":  "MONDO_0005405",
    "finngen_R12_AD_EO_EXMORE":         "MONDO_0004975",
}

OT_API_URL = "https://api.platform.opentargets.org/api/v4/graphql"
QUERY = """
{
  disease(efoId: "%s") {
    associatedTargets(page: {index: 0, size: 1000}
                      orderByScore: "score"
                      enableIndirect: true) {
      rows {
        target { approvedSymbol }
        score
        datasourceScores { id  score }
      }
    }
  }
}
"""

def fetch_ot(ot_id):
    try:
        resp = requests.post(OT_API_URL,
                             json={"query": QUERY % ot_id}, timeout=30)
        rows = resp.json()["data"]["disease"]["associatedTargets"]["rows"]
        records = []
        for r in rows:
            cp = next((ds["score"] for ds in r["datasourceScores"]
                       if ds["id"] == "clinical_precedence"), 0.0)
            records.append({
                "symbol":       r["target"]["approvedSymbol"],
                "globalScore":  r["score"],
                "cp_score":     cp if cp else 0.0,
            })
        df = pd.DataFrame(records)
        return df.groupby("symbol", as_index=False).agg(
            {"globalScore": "max", "cp_score": "max"})
    except Exception as e:
        print(f"  [API ERROR] {ot_id}: {e}")
        return pd.DataFrame()

# == 1. 데이터 수집 ============================================================
print("Reading input CSV...")
df_np = pd.read_csv(INPUT_CSV)
df_np = df_np[df_np['Trait'].isin(TRAIT_TO_OT_ID.keys())].copy()
print(f"  {len(df_np):,} rows | {df_np['Trait'].nunique()} diseases\n")

ot_cache = {}
for trait, ot_id in TRAIT_TO_OT_ID.items():
    if ot_id not in ot_cache and trait in df_np['Trait'].values:
        print(f"  Fetching {trait.replace('finngen_R12_','')} ({ot_id})...")
        ot_cache[ot_id] = fetch_ot(ot_id)
        time.sleep(0.4)

# 전체 score 모으기 (히스토그램용)
all_ot_scores = np.concatenate([d['globalScore'].values
                                  for d in ot_cache.values() if not d.empty])
all_cp_scores = np.concatenate([d['cp_score'].values
                                  for d in ot_cache.values() if not d.empty])
cp_nonzero = all_cp_scores[all_cp_scores > 0]

# == 2. Histogram: globalScore =================================================
fig, ax = plt.subplots(figsize=(10, 5))
fig.suptitle("OpenTargets globalScore Distribution", fontsize=14)

ax.hist(all_ot_scores, bins=80, color="steelblue", alpha=0.8, edgecolor="white")
for cut in OT_CUTOFFS:
    n_above = (all_ot_scores >= cut).sum()
    ax.axvline(cut, color="red", linestyle="--", linewidth=1, alpha=0.7,
               label=f"{cut} (n={n_above})")
ax.set_xlabel("globalScore", fontsize=11)
ax.set_ylabel("Count", fontsize=11)
ax.set_title("All diseases combined")
ax.legend(title="Cutoffs", fontsize=8)
ax.grid(True, alpha=0.3)

plt.tight_layout()
out1 = os.path.join(OUTPUT_DIR, "cutoff_analysis_ot_histogram.png")
plt.savefig(out1, dpi=150, bbox_inches="tight")
plt.close()
print(f"\nSaved: {out1}")

# == 3. Histogram: clinical_precedence =========================================
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
fig.suptitle("clinical_precedence Score Distribution\n(ChEMBL drug target evidence)", fontsize=13)

ax = axes[0]
ax.hist(cp_nonzero, bins=60, color="tomato", alpha=0.8, edgecolor="white")
for cut in CP_CUTOFFS:
    if cut > 0:
        n_above = (cp_nonzero >= cut).sum()
        ax.axvline(cut, color="navy", linestyle="--", linewidth=1.2,
                   label=f"{cut} (n={n_above})")
ax.set_xlabel("clinical_precedence score (> 0 only)", fontsize=11)
ax.set_ylabel("Count", fontsize=11)
ax.set_title(f"n={len(cp_nonzero)} genes with cp > 0")
ax.legend(title="CP cutoffs", fontsize=8)
ax.grid(True, alpha=0.3)

ax = axes[1]
# 0 포함 전체 -> pie chart
zero_count = int((all_cp_scores == 0).sum())
nonzero_count = int(len(cp_nonzero))
labels = ["cp = 0", "cp > 0"]
sizes = [zero_count, nonzero_count]
colors = ["lightgray", "tomato"]
ax.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%', startangle=90)
ax.axis('equal')
ax.set_title(f"Genes with/without clinical evidence\n(0: {zero_count:,} | >0: {nonzero_count:,})")

plt.tight_layout()
out2 = os.path.join(OUTPUT_DIR, "cutoff_analysis_cp_histogram.png")
plt.savefig(out2, dpi=150, bbox_inches="tight")
plt.close()
print(f"Saved: {out2}")

# == 4. Cutoff 조합별 sensitivity/precision 비교 ================================
print("\nComputing cutoff combinations...")
comparison_rows = []

for ot_cut in OT_CUTOFFS:
    for cp_cut in CP_CUTOFFS:
        tp_total = fn_total = sig_total = gold_total = 0

        for trait, ot_id in TRAIT_TO_OT_ID.items():
            if ot_id not in ot_cache or ot_cache[ot_id].empty:
                continue
            sub = df_np[df_np['Trait'] == trait]
            if len(sub) == 0:
                continue
            ot_df = ot_cache[ot_id]

            ot_gold = set(ot_df[ot_df['globalScore'] >= ot_cut]['symbol'])
            if cp_cut == 0.0:
                cp_gold = set()
            else:
                cp_gold = set(ot_df[ot_df['cp_score'] >= cp_cut]['symbol'])
            gold = ot_gold | cp_gold

            sig  = set(sub[sub['significant']]['gene'])
            tp   = len(sig & gold)
            fn   = len(gold - sig)

            tp_total   += tp
            fn_total   += fn
            sig_total  += len(sig)
            gold_total += len(gold)

        prec = tp_total / sig_total  * 100 if sig_total  > 0 else 0
        sens = tp_total / gold_total * 100 if gold_total > 0 else 0
        comparison_rows.append({
            'ot_cutoff': ot_cut,
            'cp_cutoff': cp_cut,
            'cp_label':  'none' if cp_cut == 0 else f'>={cp_cut}',
            'n_gold':    gold_total,
            'TP':        tp_total,
            'FN':        fn_total,
            'precision_pct':   round(prec, 2),
            'sensitivity_pct': round(sens, 2),
        })

comp_df = pd.DataFrame(comparison_rows)
out3 = os.path.join(OUTPUT_DIR, "cutoff_analysis_comparison.csv")
comp_df.to_csv(out3, index=False)
print(f"Saved: {out3}")
print()
print(comp_df[['ot_cutoff','cp_label','n_gold','TP',
               'precision_pct','sensitivity_pct']].to_string(index=False))

# == 5. 비교 Plot ==============================================================
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
fig.suptitle("Sensitivity & Precision by Cutoff Combination", fontsize=14)

colors = ["#1f77b4","#ff7f0e","#2ca02c","#d62728"]
cp_labels = comp_df['cp_label'].unique()

for ax, metric in zip(axes, ['sensitivity_pct', 'precision_pct']):
    for i, cp_lbl in enumerate(cp_labels):
        sub = comp_df[comp_df['cp_label'] == cp_lbl]
        ax.plot(sub['ot_cutoff'], sub[metric],
                marker='o', label=f'cp {cp_lbl}',
                color=colors[i], linewidth=2)
    ax.set_xlabel("OT globalScore cutoff", fontsize=11)
    ax.set_ylabel(f"{metric.replace('_pct','').capitalize()} (%)", fontsize=11)
    ax.set_title(metric.replace('_pct','').capitalize())
    ax.legend(title="ChEMBL cutoff", fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_xticks(OT_CUTOFFS)

plt.tight_layout()
out4 = os.path.join(OUTPUT_DIR, "cutoff_analysis_comparison.png")
plt.savefig(out4, dpi=150, bbox_inches="tight")
plt.close()
print(f"Saved: {out4}")

print("\n완료. 출력 파일:")
for f in [out1, out2, out3, out4]:
    print(f"  {f}")