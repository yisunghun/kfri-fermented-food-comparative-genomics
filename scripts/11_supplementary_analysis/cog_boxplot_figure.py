#!/usr/bin/env python3
"""
cog_boxplot_figure.py

기존 Fig. 1(COG 카테고리별 "평균" 막대그래프)은 Mann-Whitney(비모수 검정)를
쓰면서 정작 그림은 평균만 보여줘서 분포/이상치 정보가 없다는 지적을 반영,
effect size(rank-biserial |r|) 상위 카테고리들을 box plot(중앙값/IQR/이상치)
으로 다시 그린다. 25개 전부를 box plot으로 그리면 지나치게 복잡해지므로,
Supplementary Table S2에서 이미 확인한 effect size 기준 상위 10개 카테고리만
선택하고, 전체 25개 요약은 기존 막대그래프(Supplementary Fig.)로 남긴다.

사용법:
    python3 cog_boxplot_figure.py \
        --cog-ratio eggnog_cog_ratio_wide.tsv \
        --functional-group-root grouped/by_functional_group \
        --outdir cog_boxplot_result
"""
import argparse
import glob
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# effect size 기준 상위 10개 카테고리 (Supplementary Table S2 기준, |r| 내림차순)
TOP10_CATEGORIES = ["F", "N", "J", "Q", "T", "U", "D", "S", "B", "L"]
COG_LABELS = {
    "F": "F: Nucleotide\nmetabolism", "N": "N: Cell\nmotility",
    "T": "T: Signal\ntransduction", "J": "J: Translation/\nribosome",
    "Q": "Q: Secondary\nmetabolites", "L": "L: Replication/\nrepair",
    "B": "B: Chromatin\nstructure", "C": "C: Energy\nproduction",
    "V": "V: Defense\nmechanisms", "G": "G: Carbohydrate\nmetabolism",
    "U": "U: Intracellular\ntrafficking/secretion", "D": "D: Cell cycle/\ndivision",
    "S": "S: Function\nunknown",
}


def load_functional_group_map(fg_root):
    # 공식 group_manifest.tsv 기준 (폴더 스캔 방식은 Bacillus_group에서 5개
    # 누락되는 버그가 있었음 - 투고 전 자체 검증 과정에서 발견, 수정됨)
    import pandas as pd
    mapping = {}
    for fg_name in os.listdir(fg_root):
        manifest_path = os.path.join(fg_root, fg_name, "group_manifest.tsv")
        if not os.path.isfile(manifest_path):
            continue
        mdf = pd.read_csv(manifest_path, sep="\t")
        for sid in mdf["sample_id"]:
            mapping[sid] = fg_name
    return mapping


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cog-ratio", required=True)
    ap.add_argument("--functional-group-root", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    df = pd.read_csv(args.cog_ratio, sep="\t").set_index("sample_id")
    fg_map = load_functional_group_map(args.functional_group_root)
    df["functional_group"] = df.index.map(lambda s: fg_map.get(s))
    df = df[df["functional_group"].isin(["LAB", "Bacillus_group"])]

    fig, axes = plt.subplots(2, 5, figsize=(18, 8))
    axes = axes.flatten()

    for i, cat in enumerate(TOP10_CATEGORIES):
        ax = axes[i]
        lab_vals = df.loc[df["functional_group"] == "LAB", cat].dropna()
        bac_vals = df.loc[df["functional_group"] == "Bacillus_group", cat].dropna()

        bp = ax.boxplot(
            [lab_vals, bac_vals], tick_labels=["LAB", "Bacillus\n-group"],
            patch_artist=True, widths=0.5, showfliers=True,
            flierprops=dict(marker="o", markersize=3, alpha=0.5),
        )
        colors = ["#E69F00", "#0072B2"]  # 앞서 쓴 색상 팔레트와 통일
        for patch, color in zip(bp["boxes"], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.6)

        ax.set_title(COG_LABELS.get(cat, cat), fontsize=10)
        ax.set_ylabel("Relative abundance (%)", fontsize=8)
        ax.tick_params(axis="both", labelsize=8)

    plt.tight_layout()
    out_path = os.path.join(args.outdir, "cog_top10_boxplot.pdf")
    plt.savefig(out_path, dpi=200, bbox_inches="tight")
    print(f"저장됨: {out_path}")


if __name__ == "__main__":
    main()
