#!/usr/bin/env python3
"""
final_taxonomy_guild_sensitivity.py

기존 LAB(n=99)/Bacillus-group(n=79) guild 배정은 초기 BLAST 기반
genus_final(species 수준까지 못 간 경우 보수적으로 "unresolved" 처리)을
기준으로 확정된 것이다. 이후 GTDB-Tk(Section 3.3)가 이 "unresolved" 19개
중 다수를 포함해 220개 거의 전부에 확신 있는 genus/species를 부여했으므로,
최종 taxonomy 기준으로 guild를 다시 매겼을 때 LAB/Bacillus-group 비교
결론(COG/CARD/VFDB/antiSMASH)이 바뀌는지 확인한다.

사용법:
    python3 final_taxonomy_guild_sensitivity.py \
        --gtdbtk-summary gtdbtk_bac120_summary.tsv \
        --lab-manifest grouped/by_functional_group/LAB/group_manifest.tsv \
        --bacillus-manifest grouped/by_functional_group/Bacillus_group/group_manifest.tsv \
        --cog-ratio eggnog_cog_ratio_wide.tsv \
        --card-summary abricate_out/card_summary.tsv \
        --vfdb-summary abricate_out/vfdb_summary.tsv \
        --antismash-summary antismash_analysis/antismash_bgc_summary.tsv \
        --outdir final_taxonomy_guild_result
"""
import argparse
import os
import re
import pandas as pd
import numpy as np
from scipy import stats

LAB_GENERA = {
    "Lactiplantibacillus", "Levilactobacillus", "Latilactobacillus", "Lactobacillus",
    "Lactococcus", "Leuconostoc", "Weissella", "Pediococcus", "Enterococcus",
    "Tetragenococcus", "Lacticaseibacillus", "Limosilactobacillus", "Lentilactobacillus",
    "Loigolactobacillus", "Fructilactobacillus",
}
BACILLUS_GROUP_GENERA = {
    "Bacillus", "Paenibacillus", "Oceanobacillus", "Virgibacillus", "Priestia",
    "Rossellomorea", "Shouchella", "Halobacillus", "Alkalicoccobacillus", "Alkalihalobacillus",
}


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


def mwu_with_effect(lab_vals, bac_vals):
    lab_vals = np.asarray(lab_vals, dtype=float)
    bac_vals = np.asarray(bac_vals, dtype=float)
    n1, n2 = len(lab_vals), len(bac_vals)
    u, p = stats.mannwhitneyu(lab_vals, bac_vals, alternative="two-sided")
    u_y = n1 * n2 - u
    r = 1 - (2 * u_y) / (n1 * n2)
    return n1, n2, p, r, effect_label(r)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gtdbtk-summary", required=True)
    ap.add_argument("--lab-manifest", required=True)
    ap.add_argument("--bacillus-manifest", required=True)
    ap.add_argument("--cog-ratio", required=True)
    ap.add_argument("--card-summary", required=True)
    ap.add_argument("--vfdb-summary", required=True)
    ap.add_argument("--antismash-summary", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    gtdb = pd.read_csv(args.gtdbtk_summary, sep="\t")

    def parse_taxon(c, prefix):
        for part in str(c).split(";"):
            if part.startswith(prefix):
                return part[len(prefix):]
        return None
    gtdb["gtdb_genus_base"] = gtdb["classification"].apply(
        lambda c: re.sub(r"_[A-Z]+$", "", str(parse_taxon(c, "g__"))))
    gtdb["sample_id"] = gtdb["user_genome"]

    lab_ids = set(pd.read_csv(args.lab_manifest, sep="\t")["sample_id"])
    bac_ids = set(pd.read_csv(args.bacillus_manifest, sep="\t")["sample_id"])
    print(f"원래 LAB: {len(lab_ids)}, 원래 Bacillus_group: {len(bac_ids)}")

    outside = gtdb[~gtdb["sample_id"].isin(lab_ids | bac_ids)]
    add_to_lab = set(outside.loc[outside["gtdb_genus_base"].isin(LAB_GENERA), "sample_id"])
    add_to_bac = set(outside.loc[outside["gtdb_genus_base"].isin(BACILLUS_GROUP_GENERA), "sample_id"])
    print(f"GTDB-Tk 기준 LAB에 추가될 시료: {len(add_to_lab)}개 {sorted(add_to_lab)}")
    print(f"GTDB-Tk 기준 Bacillus_group에 추가될 시료: {len(add_to_bac)}개 {sorted(add_to_bac)}")

    lab_ext = lab_ids | add_to_lab
    bac_ext = bac_ids | add_to_bac
    print(f"확장판 LAB: {len(lab_ext)}, 확장판 Bacillus_group: {len(bac_ext)}\n")

    def compare(vals_df, col, label):
        orig = mwu_with_effect(
            vals_df.loc[vals_df.index.isin(lab_ids), col].dropna(),
            vals_df.loc[vals_df.index.isin(bac_ids), col].dropna())
        ext = mwu_with_effect(
            vals_df.loc[vals_df.index.isin(lab_ext), col].dropna(),
            vals_df.loc[vals_df.index.isin(bac_ext), col].dropna())
        print(f"[{label}]")
        print(f"  원래(n={orig[0]}/{orig[1]}):   p={orig[2]:.4g}, r={orig[3]:.3f} ({orig[4]})")
        print(f"  확장판(n={ext[0]}/{ext[1]}): p={ext[2]:.4g}, r={ext[3]:.3f} ({ext[4]})")
        print()

    for label, path in [("CARD", args.card_summary), ("VFDB", args.vfdb_summary)]:
        df = pd.read_csv(path, sep="\t")
        df["sample_id"] = df[df.columns[0]].apply(derive_sample_id)
        df = df.set_index("sample_id")
        compare(df, "NUM_FOUND", f"{label} burden")

    asum = pd.read_csv(args.antismash_summary, sep="\t")
    id_col = "sample_id" if "sample_id" in asum.columns else asum.columns[0]
    asum = asum.set_index(id_col)
    burden_col = "total_bgc_count" if "total_bgc_count" in asum.columns else asum.columns[0]
    compare(asum, burden_col, "antiSMASH burden")

    cog = pd.read_csv(args.cog_ratio, sep="\t").set_index("sample_id")
    cog_cols = [c for c in cog.columns if c not in
                {"genus_final", "species_final", "functional_group"} and len(c) <= 2]

    def bh_fdr(pvals):
        pvals = np.asarray(pvals, dtype=float)
        n = len(pvals)
        order = np.argsort(pvals)
        ranked = pvals[order]
        adj = ranked * n / (np.arange(n) + 1)
        adj = np.minimum.accumulate(adj[::-1])[::-1]
        return np.clip(adj, 0, 1)[np.argsort(order)]

    def cog_significant_count(lset, bset):
        ps, rs = [], []
        for cat in cog_cols:
            lv = cog.loc[cog.index.isin(lset), cat].dropna()
            bv = cog.loc[cog.index.isin(bset), cat].dropna()
            n1, n2, p, r, _ = mwu_with_effect(lv, bv)
            ps.append(p)
            rs.append(r)
        padj = bh_fdr(ps)
        return (padj < 0.05).sum(), list(zip(cog_cols, np.sign(rs)))

    n_sig_orig, dir_orig = cog_significant_count(lab_ids, bac_ids)
    n_sig_ext, dir_ext = cog_significant_count(lab_ext, bac_ext)
    print(f"[COG] 원래: {n_sig_orig}/25 유의, 확장판: {n_sig_ext}/25 유의")
    flips = [c for (c, s1), (_, s2) in zip(dir_orig, dir_ext) if s1 != s2 and s1 != 0 and s2 != 0]
    print(f"방향이 뒤집힌 카테고리: {flips if flips else '없음'}")


if __name__ == "__main__":
    main()
