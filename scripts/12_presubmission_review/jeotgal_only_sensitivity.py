#!/usr/bin/env python3
"""
jeotgal_only_sensitivity.py

food source(김치/장류/젓갈)와 functional guild(LAB/Bacillus-group)가 완전히
분리되어 있지 않다는 점(Supplementary Table S3)을 활용해, "젓갈 유래 시료만"
따로 떼어내 그 안에서도 COG/CARD/antiSMASH의 LAB vs Bacillus-group 차이가
유지되는지 확인한다. 유지된다면 "이 차이가 김치 vs 장류라는 기질 차이 때문에
생긴 것"이라는 대안 설명을 직접 배제할 수 있다.

*** 중요 ***: functional_group은 반드시 공식 group_manifest.tsv
(grouped/by_functional_group/{LAB,Bacillus_group}/group_manifest.tsv)로만
판정한다. 폴더 스캔 방식은 Bacillus_group에서 5개가 누락되는 버그가 있었음
(투고 전 자체 검증 과정에서 발견, 이미 수정됨).

사용법:
    python3 jeotgal_only_sensitivity.py \
        --wgs-summaries WGS_Summaries.xlsx \
        --lab-manifest grouped/by_functional_group/LAB/group_manifest.tsv \
        --bacillus-manifest grouped/by_functional_group/Bacillus_group/group_manifest.tsv \
        --cog-ratio eggnog_cog_ratio_wide.tsv \
        --card-summary abricate_out/card_summary.tsv \
        --vfdb-summary abricate_out/vfdb_summary.tsv \
        --antismash-summary antismash_analysis/antismash_bgc_summary.tsv \
        --genome-qc genome_qc/genome_qc_metrics.tsv \
        --outdir jeotgal_only_result
"""
import argparse
import os
import re
import pandas as pd
import numpy as np
from scipy import stats


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
    if n1 < 3 or n2 < 3:
        return {"metric": label, "n_LAB": n1, "n_Bacillus_group": n2,
                "mean_LAB": np.nan, "mean_Bacillus_group": np.nan,
                "p_value": np.nan, "rank_biserial_r": np.nan, "effect_size": "n too small"}
    u, p = stats.mannwhitneyu(lab_vals, bac_vals, alternative="two-sided")
    u_y = n1 * n2 - u
    r = 1 - (2 * u_y) / (n1 * n2)
    return {
        "metric": label, "n_LAB": n1, "n_Bacillus_group": n2,
        "mean_LAB": lab_vals.mean(), "mean_Bacillus_group": bac_vals.mean(),
        "median_LAB": np.median(lab_vals), "median_Bacillus_group": np.median(bac_vals),
        "p_value": p, "rank_biserial_r": r, "effect_size": effect_label(r),
        "direction": "LAB higher" if r > 0 else ("Bacillus_group higher" if r < 0 else "tied"),
    }


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
    ap.add_argument("--wgs-summaries", required=True)
    ap.add_argument("--lab-manifest", required=True)
    ap.add_argument("--bacillus-manifest", required=True)
    ap.add_argument("--cog-ratio", required=True)
    ap.add_argument("--card-summary", required=True)
    ap.add_argument("--vfdb-summary", required=True)
    ap.add_argument("--antismash-summary", required=True)
    ap.add_argument("--genome-qc", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # ---- 1) 공식 manifest로 fg_map 구성 (폴더 스캔 금지) ----
    lab_ids = set(pd.read_csv(args.lab_manifest, sep="\t")["sample_id"])
    bac_ids = set(pd.read_csv(args.bacillus_manifest, sep="\t")["sample_id"])
    print(f"LAB(공식): {len(lab_ids)}, Bacillus_group(공식): {len(bac_ids)}")
    assert len(lab_ids) == 99, f"LAB이 99가 아님: {len(lab_ids)}"
    assert len(bac_ids) == 79, f"Bacillus_group이 79가 아님: {len(bac_ids)}"
    print("검증 통과\n")

    # ---- 2) food source 분류 ----
    food = pd.read_excel(args.wgs_summaries, sheet_name="Sheet1")
    food.columns = [c.strip() for c in food.columns]
    food["sample_id"] = food["Order Name"].astype(str) + "_" + food["Name"].astype(str)

    def broad(src):
        src = str(src)
        if "젓갈" in src:
            return "Jeotgal"
        if src in {"된장", "고추장", "청국장", "간장", "메주"}:
            return "Jang-type"
        if src == "김치":
            return "Kimchi"
        return "Other"

    food["food_broad"] = food["분리원"].apply(broad)
    jeotgal_ids = set(food.loc[food["food_broad"] == "Jeotgal", "sample_id"])
    print(f"젓갈(jeotgal) 유래 전체 시료: {len(jeotgal_ids)}개")

    jeotgal_lab = lab_ids & jeotgal_ids
    jeotgal_bac = bac_ids & jeotgal_ids
    print(f"  이 중 LAB: {len(jeotgal_lab)}개, Bacillus_group: {len(jeotgal_bac)}개\n")

    qc_full = pd.read_csv(args.genome_qc, sep="\t")
    qc = qc_full.set_index(qc_full.columns[0])
    genome_mb = qc["total_length_bp"] / 1_000_000 if "total_length_bp" in qc.columns else None

    all_results = []

    # ---- 3) COG (25개 카테고리) ----
    cog = pd.read_csv(args.cog_ratio, sep="\t").set_index("sample_id")
    cog_cols = [c for c in cog.columns if c not in
                {"genus_final", "species_final", "functional_group"} and len(c) <= 2]
    print("=== COG (젓갈 내부만) ===")
    cog_rows = []
    for cat in cog_cols:
        lab_vals = cog.loc[cog.index.isin(jeotgal_lab), cat].dropna()
        bac_vals = cog.loc[cog.index.isin(jeotgal_bac), cat].dropna()
        r = mwu_with_effect(lab_vals, bac_vals, f"COG_{cat}")
        cog_rows.append(r)
    cog_df = pd.DataFrame(cog_rows).dropna(subset=["p_value"])
    if len(cog_df) > 0:
        cog_df["p_adj_BH"] = bh_fdr(cog_df["p_value"].values)
        cog_df["significant_FDR<0.05"] = cog_df["p_adj_BH"] < 0.05
        n_sig = cog_df["significant_FDR<0.05"].sum()
        print(f"유의한 카테고리: {n_sig}/{len(cog_df)}")
        print(cog_df[["metric", "rank_biserial_r", "effect_size", "direction", "p_adj_BH", "significant_FDR<0.05"]]
              .sort_values("p_adj_BH").to_string(index=False))
        cog_df.to_csv(os.path.join(args.outdir, "jeotgal_only_COG.csv"), index=False)
    print()

    # ---- 4) CARD, VFDB ----
    for label, path in [("CARD", args.card_summary), ("VFDB", args.vfdb_summary)]:
        df = pd.read_csv(path, sep="\t")
        df["sample_id"] = df[df.columns[0]].apply(derive_sample_id)
        df = df.set_index("sample_id")
        lab_vals = df.loc[df.index.isin(jeotgal_lab), "NUM_FOUND"].dropna()
        bac_vals = df.loc[df.index.isin(jeotgal_bac), "NUM_FOUND"].dropna()
        r = mwu_with_effect(lab_vals, bac_vals, f"{label}_burden")
        print(f"[{label}, 젓갈 내부] n_LAB={r['n_LAB']}, n_Bacillus_group={r['n_Bacillus_group']}")
        if not np.isnan(r.get("p_value", np.nan)):
            print(f"  mean_LAB={r['mean_LAB']:.3f}, mean_Bacillus_group={r['mean_Bacillus_group']:.3f}, "
                  f"p={r['p_value']:.4g}, r={r['rank_biserial_r']:.3f} ({r['effect_size']}, {r['direction']})")
        all_results.append(r)
    print()

    # ---- 5) antiSMASH ----
    asum = pd.read_csv(args.antismash_summary, sep="\t")
    id_col = "sample_id" if "sample_id" in asum.columns else asum.columns[0]
    asum = asum.set_index(id_col)
    burden_col = "total_bgc_count" if "total_bgc_count" in asum.columns else asum.columns[0]
    lab_vals = asum.loc[asum.index.isin(jeotgal_lab), burden_col].dropna()
    bac_vals = asum.loc[asum.index.isin(jeotgal_bac), burden_col].dropna()
    r = mwu_with_effect(lab_vals, bac_vals, "antiSMASH_burden")
    print(f"[antiSMASH, 젓갈 내부] n_LAB={r['n_LAB']}, n_Bacillus_group={r['n_Bacillus_group']}")
    if not np.isnan(r.get("p_value", np.nan)):
        print(f"  mean_LAB={r['mean_LAB']:.3f}, mean_Bacillus_group={r['mean_Bacillus_group']:.3f}, "
              f"p={r['p_value']:.4g}, r={r['rank_biserial_r']:.3f} ({r['effect_size']}, {r['direction']})")
    all_results.append(r)

    pd.DataFrame(all_results).to_csv(os.path.join(args.outdir, "jeotgal_only_burden_results.csv"), index=False)
    print(f"\n완료. 결과 저장: {args.outdir}/")


if __name__ == "__main__":
    main()
