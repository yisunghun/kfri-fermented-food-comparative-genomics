#!/usr/bin/env python3
"""
normalize_by_genome_size.py

CARD/VFDB 유전자 개수, antiSMASH BGC 개수는 지금까지 "genome당 절대 개수"로
비교했다. 그런데 Bacillus-group genome이 LAB보다 평균 1.7배 크므로(Section 3.1;
Supplementary Table S1), 이 차이만으로도 절대 개수 차이가 일부 설명될 수 있다.
이 스크립트는 각 지표를 "genome 크기(Mb)당 개수"로 정규화한 뒤 LAB vs
Bacillus-group을 재검정하여, genome-size 효과를 제거해도 결과가 유지되는지
확인한다.

사용법:
    python3 normalize_by_genome_size.py \
        --genome-qc genome_qc_metrics.tsv \
        --functional-group-root grouped/by_functional_group \
        --card-summary abricate_out/card_summary.tsv \
        --vfdb-summary abricate_out/vfdb_summary.tsv \
        --antismash-summary antismash_analysis/antismash_bgc_summary.tsv \
        --outdir genome_size_normalized
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


def compare(df, value_col, label, outdir):
    lab = df.loc[df["functional_group"] == "LAB", value_col].dropna()
    bac = df.loc[df["functional_group"] == "Bacillus_group", value_col].dropna()
    u, p = stats.mannwhitneyu(lab, bac, alternative="two-sided")
    print(f"[{label}] LAB n={len(lab)} (mean={lab.mean():.4f}, median={lab.median():.4f}) vs "
          f"Bacillus_group n={len(bac)} (mean={bac.mean():.4f}, median={bac.median():.4f}) "
          f"-> Mann-Whitney p = {p:.4g}")
    return {"metric": label, "n_LAB": len(lab), "n_Bacillus_group": len(bac),
            "mean_LAB": lab.mean(), "mean_Bacillus_group": bac.mean(),
            "median_LAB": lab.median(), "median_Bacillus_group": bac.median(),
            "p_value": p}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome-qc", required=True)
    ap.add_argument("--functional-group-root", required=True)
    ap.add_argument("--card-summary", required=True)
    ap.add_argument("--vfdb-summary", required=True)
    ap.add_argument("--antismash-summary", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    qc = pd.read_csv(args.genome_qc, sep="\t").set_index("sample_id")
    genome_mb = qc["total_length_bp"] / 1_000_000

    fg_map = load_functional_group_map(args.functional_group_root)
    print(f"functional_group 매핑: {len(fg_map)}개 시료\n")

    results = []

    # --- CARD ---
    card = pd.read_csv(args.card_summary, sep="\t")
    card["sample_id"] = card[card.columns[0]].apply(derive_sample_id)
    card = card.set_index("sample_id")
    card["burden_raw"] = card["NUM_FOUND"]
    card["genome_mb"] = genome_mb
    card["burden_per_mb"] = card["burden_raw"] / card["genome_mb"]
    card["functional_group"] = card.index.map(lambda s: fg_map.get(s, np.nan))
    card = card[card["functional_group"].isin(["LAB", "Bacillus_group"])]

    print("=== CARD ===")
    results.append(compare(card, "burden_raw", "CARD_raw_count", args.outdir))
    results.append(compare(card, "burden_per_mb", "CARD_per_Mb", args.outdir))
    print()

    # --- VFDB ---
    vfdb = pd.read_csv(args.vfdb_summary, sep="\t")
    vfdb["sample_id"] = vfdb[vfdb.columns[0]].apply(derive_sample_id)
    vfdb = vfdb.set_index("sample_id")
    vfdb["burden_raw"] = vfdb["NUM_FOUND"]
    vfdb["genome_mb"] = genome_mb
    vfdb["burden_per_mb"] = vfdb["burden_raw"] / vfdb["genome_mb"]
    vfdb["functional_group"] = vfdb.index.map(lambda s: fg_map.get(s, np.nan))
    vfdb = vfdb[vfdb["functional_group"].isin(["LAB", "Bacillus_group"])]

    print("=== VFDB ===")
    results.append(compare(vfdb, "burden_raw", "VFDB_raw_count", args.outdir))
    results.append(compare(vfdb, "burden_per_mb", "VFDB_per_Mb", args.outdir))
    print()

    # --- antiSMASH ---
    asum = pd.read_csv(args.antismash_summary, sep="\t")
    id_col = "sample_id" if "sample_id" in asum.columns else asum.columns[0]
    asum = asum.set_index(id_col)
    burden_col = "total_bgc_count" if "total_bgc_count" in asum.columns else asum.columns[0]
    asum["burden_raw"] = asum[burden_col]
    asum["genome_mb"] = genome_mb
    asum["burden_per_mb"] = asum["burden_raw"] / asum["genome_mb"]
    asum["functional_group"] = asum.index.map(lambda s: fg_map.get(s, np.nan))
    asum = asum[asum["functional_group"].isin(["LAB", "Bacillus_group"])]

    print("=== antiSMASH ===")
    results.append(compare(asum, "burden_raw", "antiSMASH_raw_count", args.outdir))
    results.append(compare(asum, "burden_per_mb", "antiSMASH_per_Mb", args.outdir))

    out_df = pd.DataFrame(results)
    out_df.to_csv(os.path.join(args.outdir, "genome_size_normalized_comparison.csv"), index=False)
    print(f"\n완료. 결과: {os.path.join(args.outdir, 'genome_size_normalized_comparison.csv')}")


if __name__ == "__main__":
    main()
