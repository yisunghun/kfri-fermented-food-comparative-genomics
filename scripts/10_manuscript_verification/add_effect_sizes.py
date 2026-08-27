#!/usr/bin/env python3
"""
add_effect_sizes.py

지금까지 COG(25개 카테고리), CARD, VFDB, antiSMASH 비교는 전부 p-value만
보고했다. p-value는 "차이가 유의한지"만 알려주고 "그 차이가 얼마나 큰지"는
알려주지 않으므로, Mann-Whitney U 검정의 표준 짝인 rank-biserial correlation
(Cliff's delta와 동일한 지표, 범위 -1~+1)을 모든 비교에 추가로 계산한다.

해석 기준(관례적):
    |r| < 0.1  : negligible
    0.1 <= |r| < 0.3 : small
    0.3 <= |r| < 0.5 : medium
    |r| >= 0.5 : large

사용법:
    python3 add_effect_sizes.py \
        --master master_table_qc.tsv \
        --functional-group-root grouped/by_functional_group \
        --cog-ratio eggnog_cog_ratio_wide.tsv \
        --card-summary abricate_out/card_summary.tsv \
        --vfdb-summary abricate_out/vfdb_summary.tsv \
        --antismash-summary antismash_analysis/antismash_bgc_summary.tsv \
        --outdir effect_sizes
"""
import argparse
import glob
import os
import re
import pandas as pd
import numpy as np
from scipy import stats


def load_functional_group_map(fg_root):
    mapping = {}
    for fg_name in os.listdir(fg_root):
        fg_dir = os.path.join(fg_root, fg_name)
        if not os.path.isdir(fg_dir):
            continue
        gff_dir = os.path.join(fg_dir, "genomes_gff")
        search_dir = gff_dir if os.path.isdir(gff_dir) else fg_dir
        for fpath in glob.glob(os.path.join(search_dir, "*.gff")):
            sid = os.path.splitext(os.path.basename(fpath))[0]
            mapping[sid] = fg_name
    return mapping


def derive_sample_id(raw):
    base = os.path.basename(str(raw))
    return re.sub(r"\.(tab|fna|fa|fasta|txt|csv|tsv)$", "", base)


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


def effect_label(r):
    ar = abs(r)
    if ar < 0.1:
        return "negligible"
    if ar < 0.3:
        return "small"
    if ar < 0.5:
        return "medium"
    return "large"


def mwu_with_effect(lab_vals, bac_vals, label):
    lab_vals = np.asarray(lab_vals, dtype=float)
    bac_vals = np.asarray(bac_vals, dtype=float)
    u, p = stats.mannwhitneyu(lab_vals, bac_vals, alternative="two-sided")
    n1, n2 = len(lab_vals), len(bac_vals)
    # rank-biserial correlation (LAB 기준; 양수 = LAB가 stochastically 더 큼)
    # 주의: scipy의 mannwhitneyu(x, y)가 반환하는 U는 U_x (x=lab_vals 기준)이며,
    # x가 stochastically 더 크면 U_x가 커진다. 표준 공식 r = 1 - 2*U_y/(n1*n2)
    # (U_y = n1*n2 - U_x)를 써야 "양수=LAB가 더 큼"이 정확히 성립한다.
    u_y = n1 * n2 - u
    r = 1 - (2 * u_y) / (n1 * n2)
    return {
        "metric": label, "n_LAB": n1, "n_Bacillus_group": n2,
        "mean_LAB": lab_vals.mean(), "mean_Bacillus_group": bac_vals.mean(),
        "median_LAB": np.median(lab_vals), "median_Bacillus_group": np.median(bac_vals),
        "U": u, "p_value": p, "rank_biserial_r": r, "effect_size": effect_label(r),
        "direction": "LAB higher" if r > 0 else ("Bacillus_group higher" if r < 0 else "tied"),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--master", required=True)
    ap.add_argument("--functional-group-root", required=True)
    ap.add_argument("--cog-ratio", required=True)
    ap.add_argument("--card-summary", required=True)
    ap.add_argument("--vfdb-summary", required=True)
    ap.add_argument("--antismash-summary", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    fg_map = load_functional_group_map(args.functional_group_root)
    print(f"functional_group 매핑: {len(fg_map)}개 시료\n")

    all_results = []

    # --- COG ---
    cog = pd.read_csv(args.cog_ratio, sep="\t").set_index("sample_id")
    cog["functional_group"] = cog.index.map(lambda s: fg_map.get(s))
    cog_cols = [c for c in cog.columns if c not in
                {"genus_final", "species_final", "functional_group"} and len(c) <= 2]

    print("=== COG (25개 카테고리) ===")
    cog_rows = []
    for cat in cog_cols:
        sub = cog[cog["functional_group"].isin(["LAB", "Bacillus_group"])]
        lab = sub.loc[sub["functional_group"] == "LAB", cat].dropna()
        bac = sub.loc[sub["functional_group"] == "Bacillus_group", cat].dropna()
        r = mwu_with_effect(lab, bac, f"COG_{cat}")
        cog_rows.append(r)
    cog_df = pd.DataFrame(cog_rows)
    cog_df["p_adj_BH"] = bh_fdr(cog_df["p_value"].values)
    cog_df["significant_FDR<0.05"] = cog_df["p_adj_BH"] < 0.05
    cog_df = cog_df.sort_values("p_adj_BH")
    cog_df.to_csv(os.path.join(args.outdir, "cog_effect_sizes.csv"), index=False)
    print(cog_df[["metric", "rank_biserial_r", "effect_size", "direction", "p_adj_BH", "significant_FDR<0.05"]]
          .round(3).to_string(index=False))
    all_results.append(("COG", cog_df))
    print()

    # --- CARD / VFDB ---
    for label, path in [("CARD", args.card_summary), ("VFDB", args.vfdb_summary)]:
        df = pd.read_csv(path, sep="\t")
        df["sample_id"] = df[df.columns[0]].apply(derive_sample_id)
        df = df.set_index("sample_id")
        df["functional_group"] = df.index.map(lambda s: fg_map.get(s))
        sub = df[df["functional_group"].isin(["LAB", "Bacillus_group"])]
        lab = sub.loc[sub["functional_group"] == "LAB", "NUM_FOUND"].dropna()
        bac = sub.loc[sub["functional_group"] == "Bacillus_group", "NUM_FOUND"].dropna()
        r = mwu_with_effect(lab, bac, f"{label}_burden")
        print(f"=== {label} burden ===")
        print(f"n_LAB={r['n_LAB']}, n_Bacillus_group={r['n_Bacillus_group']}")
        print(f"rank_biserial_r={r['rank_biserial_r']:.3f} ({r['effect_size']}, {r['direction']}), p={r['p_value']:.4g}")
        pd.DataFrame([r]).to_csv(os.path.join(args.outdir, f"{label.lower()}_effect_size.csv"), index=False)
        print()

    # --- antiSMASH ---
    asum = pd.read_csv(args.antismash_summary, sep="\t")
    id_col = "sample_id" if "sample_id" in asum.columns else asum.columns[0]
    asum = asum.set_index(id_col)
    burden_col = "total_bgc_count" if "total_bgc_count" in asum.columns else asum.columns[0]
    asum["functional_group"] = asum.index.map(lambda s: fg_map.get(s))
    sub = asum[asum["functional_group"].isin(["LAB", "Bacillus_group"])]
    lab = sub.loc[sub["functional_group"] == "LAB", burden_col].dropna()
    bac = sub.loc[sub["functional_group"] == "Bacillus_group", burden_col].dropna()
    r = mwu_with_effect(lab, bac, "antiSMASH_burden")
    print("=== antiSMASH burden ===")
    print(f"n_LAB={r['n_LAB']}, n_Bacillus_group={r['n_Bacillus_group']}")
    print(f"rank_biserial_r={r['rank_biserial_r']:.3f} ({r['effect_size']}, {r['direction']}), p={r['p_value']:.4g}")
    pd.DataFrame([r]).to_csv(os.path.join(args.outdir, "antismash_effect_size.csv"), index=False)

    print(f"\n완료. 결과 저장: {args.outdir}/")


if __name__ == "__main__":
    main()
