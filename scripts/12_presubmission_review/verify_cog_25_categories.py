#!/usr/bin/env python3
"""
verify_cog_25_categories.py

COG 25개 카테고리 전체에 대해 "원래 guild(n=99/79)"와
"GTDB-Tk 기준 확장판 guild(n=100/87)" 두 조건에서 Mann-Whitney U + BH-FDR을
계산하고 나란히 비교한다.

사용법:
    python3 verify_cog_25_categories.py \
        --cog-ratio eggnog_cog_ratio_wide.tsv \
        --lab-manifest grouped/by_functional_group/LAB/group_manifest.tsv \
        --bacillus-manifest grouped/by_functional_group/Bacillus_group/group_manifest.tsv
"""
import argparse
import pandas as pd
import numpy as np
from scipy import stats

ADD_TO_LAB = {"HN00200749_MK2-46"}
ADD_TO_BAC = {
    "HN00171167_F3034", "HN00175796_F3369", "HN00175796_F3370", "HN00222446_FA0508",
    "HN00251139_BMX25007", "HN00251139_KCUT25010", "HN00280011_BSO25040", "HN00280011_L11-30",
}

COG_COLS = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M",
            "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "Y", "Z"]


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


def run(cog, lset, bset, label):
    ps, rs = [], []
    for cat in COG_COLS:
        lv = cog.loc[cog.index.isin(lset), cat].dropna().astype(float)
        bv = cog.loc[cog.index.isin(bset), cat].dropna().astype(float)
        n1, n2 = len(lv), len(bv)
        u, p = stats.mannwhitneyu(lv, bv, alternative="two-sided")
        r = 1 - (2 * (n1 * n2 - u)) / (n1 * n2)
        ps.append(p)
        rs.append(r)
    padj = bh_fdr(ps)
    n_sig = int((padj < 0.05).sum())
    print(f"[{label}] 유의 카테고리: {n_sig}/25")
    for cat, p, r, pa in zip(COG_COLS, ps, rs, padj):
        sig = "YES" if pa < 0.05 else "no"
        print(f"  {cat}: r={r:+.3f}, p={p:.3g}, padj={pa:.3g}, sig={sig}")
    return dict(zip(COG_COLS, np.sign(rs))), n_sig


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cog-ratio", required=True)
    ap.add_argument("--lab-manifest", required=True)
    ap.add_argument("--bacillus-manifest", required=True)
    args = ap.parse_args()

    cog = pd.read_csv(args.cog_ratio, sep="\t").set_index("sample_id")
    lab_ids = set(pd.read_csv(args.lab_manifest, sep="\t")["sample_id"])
    bac_ids = set(pd.read_csv(args.bacillus_manifest, sep="\t")["sample_id"])
    lab_ext = lab_ids | ADD_TO_LAB
    bac_ext = bac_ids | ADD_TO_BAC

    print(f"원래: LAB={len(lab_ids)}, Bacillus_group={len(bac_ids)}")
    print(f"확장판: LAB={len(lab_ext)}, Bacillus_group={len(bac_ext)}\n")

    d1, n1 = run(cog, lab_ids, bac_ids, "원래 (n=99/79)")
    print()
    d2, n2 = run(cog, lab_ext, bac_ext, "확장판 (n=100/87)")
    print()

    flips = [c for c in COG_COLS if d1[c] != 0 and d2[c] != 0 and d1[c] != d2[c]]
    print(f"요약: 원래 {n1}/25 유의 -> 확장판 {n2}/25 유의")
    print(f"방향이 뒤집힌 카테고리: {flips if flips else '없음'}")


if __name__ == "__main__":
    main()
