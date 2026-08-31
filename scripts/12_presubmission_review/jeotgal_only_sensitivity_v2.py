#!/usr/bin/env python3
"""
jeotgal_only_sensitivity_v2.py

food-source 원본 파일(WGS_Summaries.xlsx) 재확인 결과, 이전 젓갈 내부
Bacillus-group 개수(30개)가 5개 누락되어 있었음이 확인됨(HN00226248,
HN00281774 두 배치가 파일에서 통째로 빠져 있었음). 이 스크립트는 정정된
전체 젓갈 시료 목록(LAB 35개, Bacillus-group 35개, 아래 하드코딩)으로
COG/CARD/VFDB/antiSMASH 비교를 처음부터 다시 수행한다.

사용법:
    python3 jeotgal_only_sensitivity_v2.py \
        --cog-ratio eggnog_cog_ratio_wide.tsv \
        --card-summary abricate_out/card_summary.tsv \
        --vfdb-summary abricate_out/vfdb_summary.tsv \
        --antismash-summary antismash_analysis/antismash_bgc_summary.tsv \
        --outdir jeotgal_only_result_v2
"""
import argparse
import os
import re
import pandas as pd
import numpy as np
from scipy import stats

JEOTGAL_LAB = [
    "HN00158226_GS200142", "HN00173995_T10-42", "HN00175004_GS00019", "HN00175004_GS00024",
    "HN00175796_GS201401", "HN00175796_GS201402", "HN00176442_GS40304", "HN00176442_GS41429",
    "HN00180670_MK1-10", "HN00180670_MK1-12", "HN00180670_MK1-13", "HN00180670_MK1-18",
    "HN00180670_MK1-34", "HN00180670_MK1-9", "HN00180670_MK2-12", "HN00180670_MK2-3",
    "HN00180670_MK2-6", "HN00195341_G2-10", "HN00195341_G3-29", "HN00195341_JJ1-135",
    "HN00200515_GS00003", "HN00200515_GS00008", "HN00200515_GS00049", "HN00200515_GS200102",
    "HN00200515_GS200141", "HN00200749_GS00038", "HN00200749_T10-11", "HN00222446_GS40726",
    "HN00223109_MK1-41", "HN00223109_T7-42", "HN00251139_BTH25001", "HN00252244_BMX25010",
    "HN00280011_BEG25026", "HN00280011_C4-50", "HN00280011_GS00031",
]
JEOTGAL_BAC = [
    "HN00163376_FA111", "HN00167634_C4-2", "HN00173995_MK3-11", "HN00173995_MK3-15",
    "HN00175004_GS00026", "HN00175486_GS00022", "HN00175796_GS41424", "HN00175796_GS41425",
    "HN00176442_GS40306", "HN00176442_GSY0003", "HN00176442_GSY0007", "HN00176442_GSY0017",
    "HN00204805_FA0221", "HN00222446_FA0107", "HN00222446_FA0422", "HN00222446_GSY0035",
    "HN00222446_Sea08-36", "HN00222446_T1-19", "HN00223109_TPP6047", "HN00226248_TPP3038",
    "HN00251139_BGIL25015", "HN00251139_BGIL25041", "HN00251139_BJS25031", "HN00251139_BMX25005",
    "HN00280011_BGIL25058", "HN00280011_KCUT25001", "HN00280011_KCUT25011", "HN00280011_KCUT25019",
    "HN00280011_KCUT25026", "HN00280011_L11-28", "HN00280011_Sea9-11", "HN00281774_BMX25006",
    "HN00281774_BSO25019", "HN00281774_BSO25044", "HN00281774_Sea9-10",
]


def derive_sample_id(raw):
    base = os.path.basename(str(raw))
    return re.sub(r"\.(tab|fna|fa|fasta|txt|csv|tsv)$", "", base)


def effect_label(r):
    ar = abs(r)
    if ar < 0.147:
        return "negligible"
    if ar < 0.33:
        return "small"
    if ar < 0.474:
        return "medium"
    return "large"


def mwu_with_effect(lab_vals, bac_vals, label):
    lab_vals = np.asarray(lab_vals, dtype=float)
    bac_vals = np.asarray(bac_vals, dtype=float)
    n1, n2 = len(lab_vals), len(bac_vals)
    u, p = stats.mannwhitneyu(lab_vals, bac_vals, alternative="two-sided")
    u_y = n1 * n2 - u
    r = 1 - (2 * u_y) / (n1 * n2)
    return {"metric": label, "n_LAB": n1, "n_Bacillus_group": n2,
            "mean_LAB": lab_vals.mean(), "mean_Bacillus_group": bac_vals.mean(),
            "p_value": p, "rank_biserial_r": r, "effect_size": effect_label(r),
            "direction": "LAB higher" if r > 0 else ("Bacillus_group higher" if r < 0 else "tied")}


def bh_fdr(pvals):
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order]
    adj = ranked * n / (np.arange(n) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    out = np.empty(n)
    out[order] = adj
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cog-ratio", required=True)
    ap.add_argument("--card-summary", required=True)
    ap.add_argument("--vfdb-summary", required=True)
    ap.add_argument("--antismash-summary", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    jeotgal_lab = set(JEOTGAL_LAB)
    jeotgal_bac = set(JEOTGAL_BAC)
    print(f"젓갈 LAB: {len(jeotgal_lab)}, Bacillus_group: {len(jeotgal_bac)}\n")

    all_results = []

    cog = pd.read_csv(args.cog_ratio, sep="\t").set_index("sample_id")
    cog_cols = [c for c in cog.columns if c not in
                {"genus_final", "species_final", "functional_group"} and len(c) <= 2]
    cog_rows = []
    for cat in cog_cols:
        lv = cog.loc[cog.index.isin(jeotgal_lab), cat].dropna().astype(float)
        bv = cog.loc[cog.index.isin(jeotgal_bac), cat].dropna().astype(float)
        if len(lv) < 3 or len(bv) < 3 or lv.std() == 0 or bv.std() == 0:
            continue
        cog_rows.append(mwu_with_effect(lv, bv, f"COG_{cat}"))
    cog_df = pd.DataFrame(cog_rows)
    cog_df["p_adj_BH"] = bh_fdr(cog_df["p_value"].values)
    cog_df["significant_FDR<0.05"] = cog_df["p_adj_BH"] < 0.05
    n_sig = cog_df["significant_FDR<0.05"].sum()
    print(f"=== COG (젓갈, 정정된 n={len(jeotgal_lab)}/{len(jeotgal_bac)}) ===")
    print(f"유의 카테고리: {n_sig}/{len(cog_df)}")
    print(cog_df[["metric", "rank_biserial_r", "effect_size", "direction", "p_adj_BH", "significant_FDR<0.05"]]
          .sort_values("p_adj_BH").to_string(index=False))
    cog_df.to_csv(os.path.join(args.outdir, "jeotgal_v2_COG.csv"), index=False)
    print()

    for label, path in [("CARD", args.card_summary), ("VFDB", args.vfdb_summary)]:
        df = pd.read_csv(path, sep="\t")
        df["sample_id"] = df[df.columns[0]].apply(derive_sample_id)
        df = df.set_index("sample_id")
        lv = df.loc[df.index.isin(jeotgal_lab), "NUM_FOUND"].dropna()
        bv = df.loc[df.index.isin(jeotgal_bac), "NUM_FOUND"].dropna()
        r = mwu_with_effect(lv, bv, f"{label}_burden")
        print(f"[{label}] n_LAB={r['n_LAB']}, n_Bacillus_group={r['n_Bacillus_group']}")
        print(f"  mean_LAB={r['mean_LAB']:.3f}, mean_Bacillus_group={r['mean_Bacillus_group']:.3f}, "
              f"p={r['p_value']:.4g}, r={r['rank_biserial_r']:.3f} ({r['effect_size']}, {r['direction']})")
        all_results.append(r)
    print()

    asum = pd.read_csv(args.antismash_summary, sep="\t")
    id_col = "sample_id" if "sample_id" in asum.columns else asum.columns[0]
    asum = asum.set_index(id_col)
    burden_col = "total_bgc_count" if "total_bgc_count" in asum.columns else asum.columns[0]
    lv = asum.loc[asum.index.isin(jeotgal_lab), burden_col].dropna()
    bv = asum.loc[asum.index.isin(jeotgal_bac), burden_col].dropna()
    r = mwu_with_effect(lv, bv, "antiSMASH_burden")
    print(f"[antiSMASH] n_LAB={r['n_LAB']}, n_Bacillus_group={r['n_Bacillus_group']}")
    print(f"  mean_LAB={r['mean_LAB']:.3f}, mean_Bacillus_group={r['mean_Bacillus_group']:.3f}, "
          f"p={r['p_value']:.4g}, r={r['rank_biserial_r']:.3f} ({r['effect_size']}, {r['direction']})")
    all_results.append(r)

    pd.DataFrame(all_results).to_csv(os.path.join(args.outdir, "jeotgal_v2_burden_results.csv"), index=False)
    print(f"\n완료. 결과 저장: {args.outdir}/")


if __name__ == "__main__":
    main()
