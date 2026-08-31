#!/usr/bin/env python3
"""
recompute_with_correct_fg_map.py

이전 세션에서 만든 여러 스크립트(add_effect_sizes.py, normalize_by_genome_size.py,
genus_stratified_sensitivity.py, check_nheABC_denominator.py)가 전부 "폴더 스캔"
방식으로 fg_map을 만들었는데, 이 방식이 Bacillus_group에서 79개 중 5개를
누락시켜(74개로 인식) 일부 결과에 영향을 준 것으로 확인됐다. 이 스크립트는
공식 group_manifest.tsv(각 functional group 폴더 안에 있는, faa_ok/gff_ok로
검증된 진짜 명단)를 기준으로 fg_map을 정확하게 재구성하고, 영향받은 핵심
분석들을 n=79(Bacillus_group), n=99(LAB) 기준으로 재실행한다.

사용법:
    python3 recompute_with_correct_fg_map.py \
        --functional-group-root grouped/by_functional_group \
        --master master_table_qc.tsv \
        --card-summary abricate_out/card_summary.tsv \
        --vfdb-summary abricate_out/vfdb_summary.tsv \
        --antismash-summary antismash_analysis/antismash_bgc_summary.tsv \
        --genome-qc genome_qc/genome_qc_metrics.tsv \
        --outdir corrected_reanalysis
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


def load_official_fg_map(fg_root):
    fg_map = {}
    for fg_name in os.listdir(fg_root):
        manifest_path = os.path.join(fg_root, fg_name, "group_manifest.tsv")
        if not os.path.isfile(manifest_path):
            continue
        mdf = pd.read_csv(manifest_path, sep="\t")
        mdf.columns = [re.sub(r"^\ufeff", "", c) for c in mdf.columns]
        for sid in mdf["sample_id"]:
            fg_map[sid] = fg_name
    return fg_map


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
    ap.add_argument("--functional-group-root", required=True)
    ap.add_argument("--master", required=True)
    ap.add_argument("--card-summary", required=True)
    ap.add_argument("--vfdb-summary", required=True)
    ap.add_argument("--antismash-summary", required=True)
    ap.add_argument("--genome-qc", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    fg_map = load_official_fg_map(args.functional_group_root)
    n_by_group = pd.Series(fg_map).value_counts()
    print("=== 공식 manifest 기준 functional_group 개수 (검증용) ===")
    print(n_by_group.to_string())
    assert n_by_group.get("LAB") == 99, "LAB이 99가 아님 - manifest 확인 필요"
    assert n_by_group.get("Bacillus_group") == 79, "Bacillus_group이 79가 아님 - manifest 확인 필요"
    print("검증 통과: LAB=99, Bacillus_group=79\n")

    qc = pd.read_csv(args.genome_qc, sep="\t")
    qc = qc.set_index(qc.columns[0])
    genome_mb = qc["total_length_bp"] / 1_000_000 if "total_length_bp" in qc.columns else None

    all_results = []

    for label, path, count_col in [
        ("CARD", args.card_summary, "NUM_FOUND"),
        ("VFDB", args.vfdb_summary, "NUM_FOUND"),
    ]:
        df = pd.read_csv(path, sep="\t")
        df["sample_id"] = df[df.columns[0]].apply(derive_sample_id)
        df = df.set_index("sample_id")
        df["functional_group"] = df.index.map(lambda s: fg_map.get(s))
        sub = df[df["functional_group"].isin(["LAB", "Bacillus_group"])]
        lab = sub.loc[sub["functional_group"] == "LAB", count_col].dropna()
        bac = sub.loc[sub["functional_group"] == "Bacillus_group", count_col].dropna()
        r = mwu_with_effect(lab, bac, f"{label}_burden")
        print(f"[{label}] n_LAB={r['n_LAB']}, n_Bacillus_group={r['n_Bacillus_group']}")
        print(f"  mean_LAB={r['mean_LAB']:.3f}, mean_Bacillus_group={r['mean_Bacillus_group']:.3f}, "
              f"p={r['p_value']:.4g}, r={r['rank_biserial_r']:.3f} ({r['effect_size']}, {r['direction']})")
        all_results.append(r)

        if genome_mb is not None:
            sub2 = sub.copy()
            sub2["genome_mb"] = genome_mb
            sub2["per_mb"] = sub2[count_col] / sub2["genome_mb"]
            lab_mb = sub2.loc[sub2["functional_group"] == "LAB", "per_mb"].dropna()
            bac_mb = sub2.loc[sub2["functional_group"] == "Bacillus_group", "per_mb"].dropna()
            r_mb = mwu_with_effect(lab_mb, bac_mb, f"{label}_per_Mb")
            print(f"  [{label}_per_Mb] mean_LAB={r_mb['mean_LAB']:.3f}, "
                  f"mean_Bacillus_group={r_mb['mean_Bacillus_group']:.3f}, p={r_mb['p_value']:.4g}")
            all_results.append(r_mb)
        print()

    asum = pd.read_csv(args.antismash_summary, sep="\t")
    id_col = "sample_id" if "sample_id" in asum.columns else asum.columns[0]
    asum = asum.set_index(id_col)
    burden_col = "total_bgc_count" if "total_bgc_count" in asum.columns else asum.columns[0]
    asum["functional_group"] = asum.index.map(lambda s: fg_map.get(s))
    sub = asum[asum["functional_group"].isin(["LAB", "Bacillus_group"])]
    lab = sub.loc[sub["functional_group"] == "LAB", burden_col].dropna()
    bac = sub.loc[sub["functional_group"] == "Bacillus_group", burden_col].dropna()
    r = mwu_with_effect(lab, bac, "antiSMASH_burden")
    print(f"[antiSMASH] n_LAB={r['n_LAB']}, n_Bacillus_group={r['n_Bacillus_group']}")
    print(f"  mean_LAB={r['mean_LAB']:.3f}, mean_Bacillus_group={r['mean_Bacillus_group']:.3f}, "
          f"p={r['p_value']:.4g}, r={r['rank_biserial_r']:.3f} ({r['effect_size']}, {r['direction']})\n")
    all_results.append(r)

    out_df = pd.DataFrame(all_results)
    out_df.to_csv(os.path.join(args.outdir, "corrected_burden_effect_sizes.csv"), index=False)

    vfdb_raw = pd.read_csv(args.vfdb_summary, sep="\t")
    vfdb_raw["sample_id"] = vfdb_raw[vfdb_raw.columns[0]].apply(derive_sample_id)
    nhe_cols = [c for c in vfdb_raw.columns if re.match(r"^nhe[ABC]$", c, re.IGNORECASE)]
    if len(nhe_cols) == 3:
        vfdb_raw = vfdb_raw.set_index("sample_id")

        def has_gene(v):
            return str(v).strip() not in {".", "nan", ""}
        vfdb_raw["nheABC_positive"] = vfdb_raw[nhe_cols].apply(
            lambda row: all(has_gene(row[c]) for c in nhe_cols), axis=1)
        vfdb_raw["functional_group"] = vfdb_raw.index.map(lambda s: fg_map.get(s))
        bac_all = vfdb_raw[vfdb_raw["functional_group"] == "Bacillus_group"]
        n_pos = bac_all["nheABC_positive"].sum()
        n_total = len(bac_all)
        print(f"=== nheABC 재확인 (정정된 n=79 기준) ===")
        print(f"nheABC 양성: {n_pos} / {n_total} ({100*n_pos/n_total:.1f}%)")

    print(f"\n완료. 결과 저장: {args.outdir}/corrected_burden_effect_sizes.csv")


if __name__ == "__main__":
    main()#!/usr/bin/env python3
"""
recompute_with_correct_fg_map.py

이전 세션에서 만든 여러 스크립트(add_effect_sizes.py, normalize_by_genome_size.py,
genus_stratified_sensitivity.py, check_nheABC_denominator.py)가 전부 "폴더 스캔"
방식으로 fg_map을 만들었는데, 이 방식이 Bacillus_group에서 79개 중 5개를
누락시켜(74개로 인식) 일부 결과에 영향을 준 것으로 확인됐다. 이 스크립트는
공식 group_manifest.tsv(각 functional group 폴더 안에 있는, faa_ok/gff_ok로
검증된 진짜 명단)를 기준으로 fg_map을 정확하게 재구성하고, 영향받은 핵심
분석들을 n=79(Bacillus_group), n=99(LAB) 기준으로 재실행한다.

사용법:
    python3 recompute_with_correct_fg_map.py \
        --functional-group-root grouped/by_functional_group \
        --master master_table_qc.tsv \
        --card-summary abricate_out/card_summary.tsv \
        --vfdb-summary abricate_out/vfdb_summary.tsv \
        --antismash-summary antismash_analysis/antismash_bgc_summary.tsv \
        --genome-qc genome_qc/genome_qc_metrics.tsv \
        --outdir corrected_reanalysis
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


def load_official_fg_map(fg_root):
    fg_map = {}
    for fg_name in os.listdir(fg_root):
        manifest_path = os.path.join(fg_root, fg_name, "group_manifest.tsv")
        if not os.path.isfile(manifest_path):
            continue
        mdf = pd.read_csv(manifest_path, sep="\t")
        mdf.columns = [re.sub(r"^\ufeff", "", c) for c in mdf.columns]
        for sid in mdf["sample_id"]:
            fg_map[sid] = fg_name
    return fg_map


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
    ap.add_argument("--functional-group-root", required=True)
    ap.add_argument("--master", required=True)
    ap.add_argument("--card-summary", required=True)
    ap.add_argument("--vfdb-summary", required=True)
    ap.add_argument("--antismash-summary", required=True)
    ap.add_argument("--genome-qc", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    fg_map = load_official_fg_map(args.functional_group_root)
    n_by_group = pd.Series(fg_map).value_counts()
    print("=== 공식 manifest 기준 functional_group 개수 (검증용) ===")
    print(n_by_group.to_string())
    assert n_by_group.get("LAB") == 99, "LAB이 99가 아님 - manifest 확인 필요"
    assert n_by_group.get("Bacillus_group") == 79, "Bacillus_group이 79가 아님 - manifest 확인 필요"
    print("검증 통과: LAB=99, Bacillus_group=79\n")

    qc = pd.read_csv(args.genome_qc, sep="\t")
    qc = qc.set_index(qc.columns[0])
    genome_mb = qc["total_length_bp"] / 1_000_000 if "total_length_bp" in qc.columns else None

    all_results = []

    for label, path, count_col in [
        ("CARD", args.card_summary, "NUM_FOUND"),
        ("VFDB", args.vfdb_summary, "NUM_FOUND"),
    ]:
        df = pd.read_csv(path, sep="\t")
        df["sample_id"] = df[df.columns[0]].apply(derive_sample_id)
        df = df.set_index("sample_id")
        df["functional_group"] = df.index.map(lambda s: fg_map.get(s))
        sub = df[df["functional_group"].isin(["LAB", "Bacillus_group"])]
        lab = sub.loc[sub["functional_group"] == "LAB", count_col].dropna()
        bac = sub.loc[sub["functional_group"] == "Bacillus_group", count_col].dropna()
        r = mwu_with_effect(lab, bac, f"{label}_burden")
        print(f"[{label}] n_LAB={r['n_LAB']}, n_Bacillus_group={r['n_Bacillus_group']}")
        print(f"  mean_LAB={r['mean_LAB']:.3f}, mean_Bacillus_group={r['mean_Bacillus_group']:.3f}, "
              f"p={r['p_value']:.4g}, r={r['rank_biserial_r']:.3f} ({r['effect_size']}, {r['direction']})")
        all_results.append(r)

        if genome_mb is not None:
            sub2 = sub.copy()
            sub2["genome_mb"] = genome_mb
            sub2["per_mb"] = sub2[count_col] / sub2["genome_mb"]
            lab_mb = sub2.loc[sub2["functional_group"] == "LAB", "per_mb"].dropna()
            bac_mb = sub2.loc[sub2["functional_group"] == "Bacillus_group", "per_mb"].dropna()
            r_mb = mwu_with_effect(lab_mb, bac_mb, f"{label}_per_Mb")
            print(f"  [{label}_per_Mb] mean_LAB={r_mb['mean_LAB']:.3f}, "
                  f"mean_Bacillus_group={r_mb['mean_Bacillus_group']:.3f}, p={r_mb['p_value']:.4g}")
            all_results.append(r_mb)
        print()

    asum = pd.read_csv(args.antismash_summary, sep="\t")
    id_col = "sample_id" if "sample_id" in asum.columns else asum.columns[0]
    asum = asum.set_index(id_col)
    burden_col = "total_bgc_count" if "total_bgc_count" in asum.columns else asum.columns[0]
    asum["functional_group"] = asum.index.map(lambda s: fg_map.get(s))
    sub = asum[asum["functional_group"].isin(["LAB", "Bacillus_group"])]
    lab = sub.loc[sub["functional_group"] == "LAB", burden_col].dropna()
    bac = sub.loc[sub["functional_group"] == "Bacillus_group", burden_col].dropna()
    r = mwu_with_effect(lab, bac, "antiSMASH_burden")
    print(f"[antiSMASH] n_LAB={r['n_LAB']}, n_Bacillus_group={r['n_Bacillus_group']}")
    print(f"  mean_LAB={r['mean_LAB']:.3f}, mean_Bacillus_group={r['mean_Bacillus_group']:.3f}, "
          f"p={r['p_value']:.4g}, r={r['rank_biserial_r']:.3f} ({r['effect_size']}, {r['direction']})\n")
    all_results.append(r)

    out_df = pd.DataFrame(all_results)
    out_df.to_csv(os.path.join(args.outdir, "corrected_burden_effect_sizes.csv"), index=False)

    vfdb_raw = pd.read_csv(args.vfdb_summary, sep="\t")
    vfdb_raw["sample_id"] = vfdb_raw[vfdb_raw.columns[0]].apply(derive_sample_id)
    nhe_cols = [c for c in vfdb_raw.columns if re.match(r"^nhe[ABC]$", c, re.IGNORECASE)]
    if len(nhe_cols) == 3:
        vfdb_raw = vfdb_raw.set_index("sample_id")

        def has_gene(v):
            return str(v).strip() not in {".", "nan", ""}
        vfdb_raw["nheABC_positive"] = vfdb_raw[nhe_cols].apply(
            lambda row: all(has_gene(row[c]) for c in nhe_cols), axis=1)
        vfdb_raw["functional_group"] = vfdb_raw.index.map(lambda s: fg_map.get(s))
        bac_all = vfdb_raw[vfdb_raw["functional_group"] == "Bacillus_group"]
        n_pos = bac_all["nheABC_positive"].sum()
        n_total = len(bac_all)
        print(f"=== nheABC 재확인 (정정된 n=79 기준) ===")
        print(f"nheABC 양성: {n_pos} / {n_total} ({100*n_pos/n_total:.1f}%)")

    print(f"\n완료. 결과 저장: {args.outdir}/corrected_burden_effect_sizes.csv")


if __name__ == "__main__":
    main()
