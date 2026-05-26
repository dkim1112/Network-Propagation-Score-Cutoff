#!/usr/bin/env python3
"""
Gold Standard Validation
목적: empirical p-value가 유의한 유전자들이 실제로 의미있는지 검증

Gold Standard 기준 (2가지 + AI, 조정 가능):
1. OpenTargets globalScore >= OT_SCORE_CUTOFF
2. clinical_precedence > 0 (ChEMBL 기반 임상/약물 근거)
   ※ OpenTargets API에서 "chembl" datasource가 별도로 존재하지 않음.
      clinical_precedence 컬럼이 ChEMBL approved drug target을 포함하는
      임상 근거 기반 점수임.
3. AI (ChatGPT) known genes: AI_GOLD 딕셔너리에 수동 입력

사용법:
    pip install requests pandas
    python3 w3_analysis_code/gold_standard.py [input_csv]
"""

import pandas as pd
import numpy as np
import requests
import sys, os, time

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# =============================================================================
# ★ 조정 가능한 파라미터 ★
# =============================================================================
INPUT_CSV       = os.path.normpath(os.path.join(SCRIPT_DIR, "..", "result_np_cutoff", "np_cutoff_results.csv"))
OUTPUT_DIR      = SCRIPT_DIR

OT_SCORE_CUTOFF = 0.3   # OpenTargets globalScore 기준 (0~1)
                         # 0.1=상위30%, 0.3=상위12%, 0.5=상위2.6%

COMBINE_MODE    = "OR"  # "OR": OT score OR clinical_precedence
                         # "AND": 둘 다 해당
# =============================================================================

TEST_TRAITS = [
    "finngen_R12_AUTOIMMUNE",
    "finngen_R12_I9_CHD",
    "finngen_R12_I9_HYPTENS",
    "finngen_R12_K11_IBD_STRICT",
    "finngen_R12_L12_ATOPIC",
    "finngen_R12_T1D",
    "finngen_R12_T2D_WIDE",
    "finngen_R12_CARDIAC_ARRHYTM",
    "finngen_R12_C3_PROSTATE_EXALLC",
    "finngen_R12_C3_BASAL_CELL_CARCINOMA_EXALLC",
    "finngen_R12_C3_SKIN_EXALLC",
    "finngen_R12_ALLERG_ASTHMA",
    "finngen_R12_ASTHMA_CHILD_EXMORE",
    "finngen_R12_AD_EO_EXMORE",
]

if len(sys.argv) > 1:
    INPUT_CSV = sys.argv[1]
    if not os.path.isabs(INPUT_CSV):
        search_paths = [
            os.path.normpath(os.path.join(os.getcwd(), INPUT_CSV)),
            os.path.normpath(os.path.join(SCRIPT_DIR, INPUT_CSV)),
            os.path.normpath(os.path.join(SCRIPT_DIR, "..", INPUT_CSV)),
            os.path.normpath(os.path.join(SCRIPT_DIR, "..", "result_np_cutoff", INPUT_CSV)),
        ]
        for candidate in search_paths:
            if os.path.exists(candidate):
                INPUT_CSV = candidate
                break
if not os.path.exists(INPUT_CSV):
    print(f"ERROR: {INPUT_CSV} not found"); sys.exit(1)

# == FinnGen trait → OpenTargets EFO ID =======================================
TRAIT_TO_OT_ID = {
    "finngen_R12_AUTOIMMUNE":                            "EFO_0005140", # autoimmune disease
    "finngen_R12_I9_CHD":                                "EFO_0001645", # coronary heart disease
    "finngen_R12_I9_HYPTENS":                            "EFO_0000537", # hypertension
    "finngen_R12_K11_IBD_STRICT":                        "EFO_0003767", # inflammatory bowel disease
    "finngen_R12_L12_ATOPIC":                            "EFO_0000274", # atopic dermatitis
    "finngen_R12_T1D":                                   "MONDO_0005147", # type 1 diabetes
    "finngen_R12_T2D_WIDE":                              "MONDO_0005148", # type 2 diabetes
    "finngen_R12_CARDIAC_ARRHYTM":                       "EFO_0004269", # cardiac arrhythmia
    "finngen_R12_C3_PROSTATE_EXALLC":                    "MONDO_0008315", # prostate cancer
    "finngen_R12_C3_BASAL_CELL_CARCINOMA_EXALLC":        "EFO_0004193", # basal cell carcinoma
    "finngen_R12_C3_SKIN_EXALLC":                        "MONDO_0002898", # skin cancer
    "finngen_R12_ALLERG_ASTHMA":                         "MONDO_0004784", # allergic asthma
    "finngen_R12_ASTHMA_CHILD_EXMORE":                   "MONDO_0005405", # childhood asthma
    "finngen_R12_AD_EO_EXMORE":                          "MONDO_0004975" # Alzheimer's disease early onset
}

# == AI gold standard (수동 입력 - 추후 교체 가능) ============================
AI_GOLD = {
    "finngen_R12_AUTOIMMUNE":           [],
    "finngen_R12_I9_CHD":               [],
    "finngen_R12_I9_HYPTENS":           [],
    "finngen_R12_K11_IBD_STRICT":       [],
    "finngen_R12_L12_ATOPIC":           [],
    "finngen_R12_T1D":                  [],
    "finngen_R12_T2D_WIDE":             [],
    "finngen_R12_CARDIAC_ARRHYTM":      [],
    "finngen_R12_C3_PROSTATE_EXALLC":   [],
    "finngen_R12_C3_BASAL_CELL_CARCINOMA_EXALLC":        [],
    "finngen_R12_C3_SKIN_EXALLC":       [],
    "finngen_R12_ALLERG_ASTHMA":        [],
    "finngen_R12_ASTHMA_CHILD_EXMORE":  [],
    "finngen_R12_AD_EO_EXMORE":         [],
}

# == OpenTargets GraphQL ======================================================
OT_API_URL = "https://api.platform.opentargets.org/api/v4/graphql"
QUERY = """
{
  disease(efoId: "%s") {
    associatedTargets(page: {index: 0, size: %d}
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

def fetch_ot_targets(ot_id, max_targets=1000):
    try:
        resp = requests.post(
            OT_API_URL,
            json={"query": QUERY % (ot_id, max_targets)},
            timeout=30
        )
        rows = resp.json()["data"]["disease"]["associatedTargets"]["rows"]
        records = []
        for r in rows:
            cp_score = next(
                (ds["score"] for ds in r["datasourceScores"]
                 if ds["id"] == "clinical_precedence"),
                None
            )
            records.append({
                "symbol":          r["target"]["approvedSymbol"],
                "globalScore":     round(r["score"], 4),
                "has_clinical":    cp_score is not None and cp_score > 0,
                "clinical_score":  cp_score if cp_score else 0.0,
            })
        if not records:
            return pd.DataFrame()

        return (
            pd.DataFrame(records)
              .groupby("symbol", as_index=False)
              .agg({
                  "globalScore": "max",
                  "has_clinical": "max",
                  "clinical_score": "max",
              })
        )
    except Exception as e:
        print(f"  [API ERROR] {ot_id}: {e}")
        return pd.DataFrame()

# == 메인 =====================================================================
print(f"Reading {INPUT_CSV}...")
df = pd.read_csv(INPUT_CSV)
df = df[df['Trait'].isin(TEST_TRAITS)].copy()
print(f"  {len(df):,} rows | {df['Trait'].nunique()} diseases")
print(f"\n파라미터: OT_SCORE_CUTOFF={OT_SCORE_CUTOFF} | COMBINE_MODE={COMBINE_MODE}\n")

all_results, summary_rows = [], []
ot_cache = {}

for trait in sorted(df['Trait'].unique()):
    short = trait.replace("finngen_R12_", "")
    print(f"\n{'='*55}\n{short}")

    ot_id = TRAIT_TO_OT_ID.get(trait)
    if not ot_id:
        print("  [SKIP] EFO ID 없음"); continue

    sub = df[df['Trait'] == trait].copy()

    # OpenTargets API
    if ot_id not in ot_cache:
        print(f"  API 호출: {ot_id}")
        ot_cache[ot_id] = fetch_ot_targets(ot_id)
        time.sleep(0.5)
    ot_df = ot_cache[ot_id]

    if ot_df.empty:
        print("  [SKIP] OT 결과 없음"); continue

    # Gold standard 집합
    ot_gold       = set(ot_df[ot_df['globalScore'] >= OT_SCORE_CUTOFF]['symbol'])
    clinical_gold = set(ot_df[ot_df['has_clinical']]['symbol'])
    ai_gold       = set(AI_GOLD.get(trait, []))

    if COMBINE_MODE == "OR":
        db_gold = ot_gold | clinical_gold
    else:
        db_gold = ot_gold & clinical_gold

    any_gold = db_gold | ai_gold

    print(f"  OT (>={OT_SCORE_CUTOFF}): {len(ot_gold):4d} | "
          f"clinical_precedence: {len(clinical_gold):4d} | "
          f"AI: {len(ai_gold):4d} | "
          f"Any: {len(any_gold):4d}")

    # 유전자별 표시
    sub = sub.copy()
    sub['ot_score']       = sub['gene'].map(ot_df.set_index('symbol')['globalScore'])
    sub['ot_gold']        = sub['gene'].isin(ot_gold)
    sub['clinical_gold']  = sub['gene'].isin(clinical_gold)
    sub['ai_gold']        = sub['gene'].isin(ai_gold)
    sub['any_gold']       = sub['gene'].isin(any_gold)

    # 카테고리
    def cat(row):
        s, g = row['significant'], row['any_gold']
        if   s and     g: return "TP"
        elif s and not g: return "Novel_sig"
        elif not s and g: return "FN"
        else:             return "TN"
    sub['category'] = sub.apply(cat, axis=1)
    all_results.append(sub)

    # 통계
    n_sig  = int(sub['significant'].sum())
    n_gold = int(sub['any_gold'].sum())
    tp     = int((sub['category'] == 'TP').sum())
    novel  = int((sub['category'] == 'Novel_sig').sum())
    fn     = int((sub['category'] == 'FN').sum())
    prec   = tp / n_sig  * 100 if n_sig  > 0 else 0
    sens   = tp / n_gold * 100 if n_gold > 0 else 0

    print(f"  Sig={n_sig} | TP={tp} | Novel_sig={novel} | FN={fn}")
    print(f"  Precision={prec:.1f}% | Sensitivity={sens:.1f}%")

    summary_rows.append({
        'Trait': short, 'n_seeds': int(sub['n_seeds'].iloc[0]),
        'n_significant': n_sig,
        'n_ot_gold': len(ot_gold),
        'n_clinical_gold': len(clinical_gold),
        'n_ai_gold': len(ai_gold),
        'n_any_gold': n_gold,
        'TP': tp, 'Novel_sig': novel, 'FN': fn,
        'precision_pct': round(prec, 2),
        'sensitivity_pct': round(sens, 2),
        'ot_cutoff': OT_SCORE_CUTOFF,
        'combine_mode': COMBINE_MODE,
    })

# == 저장 =====================================================================
final   = pd.concat(all_results, ignore_index=True)
summary = pd.DataFrame(summary_rows)
final.to_csv(  os.path.join(OUTPUT_DIR, "gold_standard_validation.csv"), index=False)
summary.to_csv(os.path.join(OUTPUT_DIR, "gold_standard_summary.csv"),    index=False)

print("\n" + "="*65)
print(summary[[
    'Trait','n_seeds','n_significant',
    'n_ot_gold','n_clinical_gold',
    'TP','Novel_sig','FN',
    'precision_pct','sensitivity_pct'
]].to_string(index=False))
print("\n결과 저장: gold_standard_validation.csv, gold_standard_summary.csv")