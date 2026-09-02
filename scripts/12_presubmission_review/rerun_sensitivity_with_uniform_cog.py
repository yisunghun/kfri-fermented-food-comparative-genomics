#!/usr/bin/env python3
"""
rerun_sensitivity_with_uniform_cog.py

균일 eggNOG-mapper 재annotation(Section 2.6) 이후, dereplication(근연 시료
1개 대표 선정) 및 genus-stratified(속 단위 평균) 민감도 분석을 COG에 대해
새 데이터로 재실행한다. CARD/VFDB/antiSMASH는 COG 재annotation과 무관하므로
기존 검증된 수치를 그대로 사용하되, 동일한 시료 집합(dereplicated set,
genus-level)에 대해 함께 재확인한다.

사용법:
    python3 rerun_sensitivity_with_uniform_cog.py \
        --ani-matrix ani_analysis/ani_matrix.csv \
        --cog-ratio eggnog_reannotation/eggnog_cog_ratio_wide_UNIFORM.tsv \
        --card-summary abricate_out/card_summary.tsv \
        --vfdb-summary abricate_out/vfdb_summary.tsv \
        --antismash-summary antismash_analysis/antismash_bgc_summary.tsv \
        --lab-manifest grouped/by_functional_group/LAB/group_manifest.tsv \
        --bacillus-manifest grouped/by_functional_group/Bacillus_group/group_manifest.tsv \
        --outdir sensitivity_uniform_cog_result
"""
import argparse
import os
import re

import networkx as nx
import numpy as np
import pandas as pd
from scipy import stats

COG_COLS = list("ABCDEFGHIJKLMNOPQRSTUVWYZ")


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


def mwu_effect(lv, bv):
    lv = np.asarray(lv, dtype=float)
    bv = np.asarray(bv, dtype=float)
    n1, n2 = len(lv), len(bv)
    u, p = stats.mannwhitneyu(lv, bv, alternative="two-sided")
    r = 1 - (2 * (n1 * n2 - u)) / (n1 * n2)
    return n1, n2, p, r


def run_cog_comparison(cog_df, lab_ids, bac_ids, label):
    rows = []
    excluded = []
    for cat in COG_COLS:
        lv = cog_df.loc[cog_df.index.isin(lab_ids), cat].dropna().astype(float)
        bv = cog_df.loc[cog_df.index.isin(bac_ids), cat].dropna().astype(float)
        combined = pd.concat([lv, bv])
        if combined.std() == 0 or len(lv) < 3 or len(bv) < 3:
            excluded.append(cat)
            continue
        n1, n2, p, r = mwu_effect(lv, bv)
        rows.append({"cat": cat, "n1": n1, "n2": n2, "p": p, "r": r, "effect": effect_label(r)})
    df = pd.DataFrame(rows)
    df["p_adj"] = bh_fdr(df["p"].values)
    df["sig"] = df["p_adj"] < 0.05
    n_sig = df["sig"].sum()
    n_large_sig = ((df["effect"] == "large") & df["sig"]).sum()
    print(f"\n=== COG: {label} ===")
    print(f"제외(분산0/표본부족): {excluded}")
    print(f"유의: {n_sig}/{len(df)}, large+유의: {n_large_sig}")
    non_sig = df.loc[~df["sig"], "cat"].tolist()
    print(f"비유의 카테고리: {non_sig}")
    return df, excluded, n_sig, non_sig


def build_dereplicated_set(ani_matrix_path, threshold=99.9):
    ani = pd.read_csv(ani_matrix_path, index_col=0)
    ani.index = ani.index.astype(str)
    ani.columns = ani.columns.astype(str)
    G = nx.Graph()
    G.add_nodes_from(ani.index)
    for i in ani.index:
        for j in ani.columns:
            if i == j:
                continue
            val = ani.loc[i, j]
            if pd.notna(val) and val >= threshold:
                G.add_edge(i, j)
    groups = list(nx.connected_components(G))
    multi_groups = [g for g in groups if len(g) > 1]
    n_in_groups = sum(len(g) for g in multi_groups)
    print(f"근연(ANI>={threshold}%) 그룹: {len(multi_groups)}개, 포함 시료: {n_in_groups}개")

    representatives = set()
    for g in groups:
        rep = sorted(g)[0]
        representatives.add(rep)
    print(f"대표 시료 선정 후 남은 시료: {len(representatives)}개")
    return representatives, multi_groups


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ani-matrix", required=True)
    ap.add_argument("--cog-ratio", required=True)
    ap.add_argument("--card-summary", required=True)
    ap.add_argument("--vfdb-summary", required=True)
    ap.add_argument("--antismash-summary", required=True)
    ap.add_argument("--lab-manifest", required=True)
    ap.add_argument("--bacillus-manifest", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    lab_ids_full = set(pd.read_csv(args.lab_manifest, sep="\t")["sample_id"])
    bac_ids_full = set(pd.read_csv(args.bacillus_manifest, sep="\t")["sample_id"])
    print(f"원래 LAB: {len(lab_ids_full)}, Bacillus_group: {len(bac_ids_full)}")

    cog = pd.read_csv(args.cog_ratio, sep="\t").set_index("sample_id")

    # ---------------- 1) Dereplication ----------------
    representatives, multi_groups = build_dereplicated_set(args.ani_matrix)
    lab_derep = lab_ids_full & representatives
    bac_derep = bac_ids_full & representatives
    print(f"Dereplicated LAB: {len(lab_derep)}, Bacillus_group: {len(bac_derep)}")

    run_cog_comparison(cog, lab_derep, bac_derep, f"Dereplicated (n={len(lab_derep)}/{len(bac_derep)})")

    card = pd.read_csv(args.card_summary, sep="\t")
    card["sample_id"] = card[card.columns[0]].apply(derive_sample_id)
    card = card.set_index("sample_id")
    lv = card.loc[card.index.isin(lab_derep), "NUM_FOUND"].dropna()
    bv = card.loc[card.index.isin(bac_derep), "NUM_FOUND"].dropna()
    n1, n2, p, r = mwu_effect(lv, bv)
    print(f"[CARD, dereplicated] n={n1}/{n2}, mean={lv.mean():.2f}/{bv.mean():.2f}, p={p:.3g}, r={r:.3f}")

    vfdb = pd.read_csv(args.vfdb_summary, sep="\t")
    vfdb["sample_id"] = vfdb[vfdb.columns[0]].apply(derive_sample_id)
    vfdb = vfdb.set_index("sample_id")
    lv = vfdb.loc[vfdb.index.isin(lab_derep), "NUM_FOUND"].dropna()
    bv = vfdb.loc[vfdb.index.isin(bac_derep), "NUM_FOUND"].dropna()
    n1, n2, p, r = mwu_effect(lv, bv)
    print(f"[VFDB, dereplicated] n={n1}/{n2}, mean={lv.mean():.2f}/{bv.mean():.2f}, p={p:.3g}, r={r:.3f}")

    asum = pd.read_csv(args.antismash_summary, sep="\t")
    id_col = "sample_id" if "sample_id" in asum.columns else asum.columns[0]
    asum = asum.set_index(id_col)
    burden_col = "total_bgc_count" if "total_bgc_count" in asum.columns else asum.columns[0]
    lv = asum.loc[asum.index.isin(lab_derep), burden_col].dropna()
    bv = asum.loc[asum.index.isin(bac_derep), burden_col].dropna()
    n1, n2, p, r = mwu_effect(lv, bv)
    print(f"[antiSMASH, dereplicated] n={n1}/{n2}, mean={lv.mean():.2f}/{bv.mean():.2f}, p={p:.3g}, r={r:.3f}")

    # ---------------- 2) Genus-stratified (genus-level means) ----------------
    print("\n\n=== Genus-stratified 분석 준비 ===")
    lab_manifest_df = pd.read_csv(args.lab_manifest, sep="\t")
    bac_manifest_df = pd.read_csv(args.bacillus_manifest, sep="\t")
    genus_col = "genus_final" if "genus_final" in lab_manifest_df.columns else None
    if genus_col is None:
        print("!! group_manifest.tsv에 genus_final 컬럼이 없습니다 - 별도 마스터 파일과 조인이 필요할 수 있습니다.")
    else:
        lab_manifest_df = lab_manifest_df.set_index("sample_id")
        bac_manifest_df = bac_manifest_df.set_index("sample_id")

        def genus_means(cog_df, manifest_df, min_n=3):
            genus_series = manifest_df[genus_col].rename("_genus_for_grouping")
            merged = cog_df.join(genus_series, how="inner")
            counts = merged["_genus_for_grouping"].value_counts()
            keep_genera = counts[counts >= min_n].index
            sub = merged[merged["_genus_for_grouping"].isin(keep_genera)]
            return sub.groupby("_genus_for_grouping")[COG_COLS].mean(), len(keep_genera)

        lab_genus_means, n_lab_genera = genus_means(cog, lab_manifest_df)
        bac_genus_means, n_bac_genera = genus_means(cog, bac_manifest_df)
        print(f"LAB 속 개수(n>=3): {n_lab_genera}, Bacillus_group 속 개수(n>=3): {n_bac_genera}")

        rows = []
        excluded = []
        for cat in COG_COLS:
            lv = lab_genus_means[cat].dropna()
            bv = bac_genus_means[cat].dropna()
            combined = pd.concat([lv, bv])
            if combined.std() == 0 or len(lv) < 3 or len(bv) < 3:
                excluded.append(cat)
                continue
            n1, n2, p, r = mwu_effect(lv, bv)
            rows.append({"cat": cat, "n1": n1, "n2": n2, "p": p, "r": r, "effect": effect_label(r)})
        gdf = pd.DataFrame(rows)
        if len(gdf) > 0:
            gdf["p_adj"] = bh_fdr(gdf["p"].values)
            gdf["sig"] = gdf["p_adj"] < 0.05
            print(f"\n=== COG: Genus-stratified (n_LAB genera={n_lab_genera}, n_Bacillus genera={n_bac_genera}) ===")
            print(f"제외: {excluded}")
            print(f"유의: {gdf['sig'].sum()}/{len(gdf)}")
            print(gdf.sort_values('p_adj').to_string(index=False))
            gdf.to_csv(os.path.join(args.outdir, "genus_stratified_COG_uniform.csv"), index=False)

    print(f"\n완료. 결과는 {args.outdir}/ 에 저장됨")


if __name__ == "__main__":
    main()