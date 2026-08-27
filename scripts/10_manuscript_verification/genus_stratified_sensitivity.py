#!/usr/bin/env python3
"""
genus_stratified_sensitivity.py

LAB vs Bacillus-group 비교(COG, CARD, VFDB, antiSMASH)가 특정 속(예: Bacillus
n=63, Enterococcus n=23)의 시료 수 편중에 의해 좌우되는 것은 아닌지 검증한다.

방법: isolate 단위가 아니라 "genus 단위 평균"을 계산해, 각 속이 시료 수와
무관하게 동일한 가중치(1속=1값)를 갖도록 한 뒤 LAB-genera vs
Bacillus-group-genera를 Mann-Whitney U로 재검정한다.
(n>=3 시료를 가진 genus만 포함 - 나머지 파이프라인과 기준 일치)

사용법:
    python3 genus_stratified_sensitivity.py \
        --master master_table_qc.tsv \
        --cog-ratio eggnog_cog_ratio_wide.tsv \
        --card-summary abricate_out/card_summary.tsv \
        --vfdb-summary abricate_out/vfdb_summary.tsv \
        --antismash-summary antismash_analysis/antismash_bgc_summary.tsv \
        --outdir genus_stratified_result \
        --min-n 3
"""
import argparse
import glob
import os
import re
import pandas as pd
import numpy as np
from scipy import stats


def load_functional_group_map(fg_root):
    """grouped/by_functional_group/<LAB|Bacillus_group|...>/ 폴더 구조를 근거로
    sample_id -> functional_group 매핑을 만든다 (master_table_qc.tsv에는 이
    컬럼이 없는 것으로 확인되어, 짐작 대신 이미 확정된 분류를 그대로 사용)."""
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
        manifest_path = os.path.join(fg_dir, "group_manifest.tsv")
        if os.path.isfile(manifest_path):
            try:
                mdf = pd.read_csv(manifest_path, sep="\t")
                mdf.columns = [re.sub(r"^\ufeff", "", c) for c in mdf.columns]
                if "sample_id" in mdf.columns:
                    for sid in mdf["sample_id"]:
                        mapping[sid] = fg_name
            except Exception:
                pass
    return mapping


def derive_sample_id(raw_value):
    base = os.path.basename(str(raw_value))
    return re.sub(r"\.(tab|fna|fa|fasta|txt|csv|tsv)$", "", base)


def load_master(path):
    master = pd.read_csv(path, sep="\t")
    master.columns = [re.sub(r"^\ufeff", "", c) for c in master.columns]  # BOM 제거
    return master.set_index("sample_id")


def genus_level_means(value_series_by_sample, master, fg_map, min_n, value_name):
    """sample_id -> 값 매핑을 genus_final 기준으로 평균낸 DataFrame 반환"""
    df = value_series_by_sample.rename(value_name).to_frame()
    df["genus_final"] = df.index.map(lambda s: master.loc[s, "genus_final"] if s in master.index else np.nan)
    df["functional_group"] = df.index.map(lambda s: fg_map.get(s, np.nan))
    df = df.dropna(subset=["genus_final", "functional_group"])
    df = df[df["functional_group"].isin(["LAB", "Bacillus_group"])]

    genus_n = df.groupby("genus_final").size()
    valid_genera = genus_n[genus_n >= min_n].index
    df = df[df["genus_final"].isin(valid_genera)]

    genus_means = df.groupby(["genus_final", "functional_group"])[value_name].mean().reset_index()
    return genus_means


def compare_genus_means(genus_means, value_name, label):
    lab = genus_means[genus_means["functional_group"] == "LAB"][value_name]
    bac = genus_means[genus_means["functional_group"] == "Bacillus_group"][value_name]
    n_lab_genera = len(lab)
    n_bac_genera = len(bac)
    if n_lab_genera < 2 or n_bac_genera < 2:
        print(f"[{label}] 속 개수가 너무 적어(LAB {n_lab_genera}속, Bacillus_group {n_bac_genera}속) 검정 생략")
        return None
    u, p = stats.mannwhitneyu(lab, bac, alternative="two-sided")
    print(f"[{label}] LAB {n_lab_genera}개 속(평균 {lab.mean():.3f}) vs "
          f"Bacillus_group {n_bac_genera}개 속(평균 {bac.mean():.3f}) "
          f"-> Mann-Whitney p = {p:.4g}")
    return {"metric": label, "n_lab_genera": n_lab_genera, "n_bacillus_genera": n_bac_genera,
            "mean_lab": lab.mean(), "mean_bacillus_group": bac.mean(), "p_value": p}


def bh_fdr(pvals):
    """의존성 없는 수동 Benjamini-Hochberg FDR 보정"""
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order]
    adj = ranked * n / (np.arange(n) + 1)
    # 단조 감소 방향으로 누적 최소값 적용 (표준 BH 절차)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    out = np.empty(n)
    out[order] = adj
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--master", required=True)
    ap.add_argument("--functional-group-root", required=True,
                     help="grouped/by_functional_group 폴더 경로 (LAB/Bacillus_group 하위폴더 포함)")
    ap.add_argument("--cog-ratio", required=True)
    ap.add_argument("--card-summary", required=True)
    ap.add_argument("--vfdb-summary", required=True)
    ap.add_argument("--antismash-summary", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--min-n", type=int, default=3)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    master = load_master(args.master)
    fg_map = load_functional_group_map(args.functional_group_root)
    print(f"functional_group 매핑 로드됨: {len(fg_map)}개 시료")

    results = []

    print("=== 어느 속들이 genus 단위 비교에 포함되는지 ===")
    tmp = pd.Series(fg_map, name="functional_group").to_frame()
    tmp["genus_final"] = tmp.index.map(lambda s: master.loc[s, "genus_final"] if s in master.index else np.nan)
    genus_n = tmp.groupby(["functional_group", "genus_final"]).size()
    print(genus_n[genus_n >= args.min_n].to_string())
    print()

    # --- COG: 카테고리별로 genus-stratified 비교 ---
    cog = pd.read_csv(args.cog_ratio, sep="\t").set_index("sample_id")
    cog_cols = [c for c in cog.columns if c not in
                {"genus_final", "species_final", "functional_group"} and len(c) <= 2]
    print("=== COG 카테고리별 genus-stratified 비교 (FDR 보정 전 p-value) ===")
    cog_rows = []
    for cat in cog_cols:
        gm = genus_level_means(cog[cat], master, fg_map, args.min_n, cat)
        r = compare_genus_means(gm, cat, f"COG_{cat}")
        if r:
            cog_rows.append(r)
    if cog_rows:
        cog_df = pd.DataFrame(cog_rows)
        cog_df["p_adj_BH"] = bh_fdr(cog_df["p_value"].values)
        cog_df["significant_FDR<0.05"] = cog_df["p_adj_BH"] < 0.05
        cog_df.to_csv(os.path.join(args.outdir, "cog_genus_stratified.csv"), index=False)
        n_sig = cog_df["significant_FDR<0.05"].sum()
        print(f"\nCOG genus-stratified: {n_sig}/{len(cog_df)} 카테고리 유의 (FDR<0.05)")
        results.append(("COG", n_sig, len(cog_df)))

    # --- CARD/VFDB burden: genus-stratified ---
    for label, path in [("CARD", args.card_summary), ("VFDB", args.vfdb_summary)]:
        df = pd.read_csv(path, sep="\t")
        id_col = df.columns[0]
        df["sample_id"] = df[id_col].apply(derive_sample_id)
        df = df.set_index("sample_id")
        burden_col = "NUM_FOUND" if "NUM_FOUND" in df.columns else None
        if burden_col is None:
            print(f"[{label}] NUM_FOUND 컬럼을 못 찾음, 건너뜀")
            continue
        gm = genus_level_means(df[burden_col], master, fg_map, args.min_n, "burden")
        r = compare_genus_means(gm, "burden", f"{label}_burden")
        if r:
            pd.DataFrame([r]).to_csv(os.path.join(args.outdir, f"{label.lower()}_genus_stratified.csv"), index=False)

    # --- antiSMASH burden: genus-stratified ---
    asum = pd.read_csv(args.antismash_summary, sep="\t")
    id_col = "sample_id" if "sample_id" in asum.columns else asum.columns[0]
    asum = asum.set_index(id_col)
    burden_col = "total_bgc_count" if "total_bgc_count" in asum.columns else asum.columns[0]
    gm = genus_level_means(asum[burden_col], master, fg_map, args.min_n, "bgc_burden")
    r = compare_genus_means(gm, "bgc_burden", "antiSMASH_burden")
    if r:
        pd.DataFrame([r]).to_csv(os.path.join(args.outdir, "antismash_genus_stratified.csv"), index=False)

    print(f"\n완료. 결과 파일들: {args.outdir}/")


if __name__ == "__main__":
    main()
