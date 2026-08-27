#!/usr/bin/env python3
"""
cog_pca_clr.py

COG 카테고리 비율(%)은 시료당 합이 100%로 고정되는 compositional data이다.
Raw percentage에 그냥 PCA/유클리드 기반 분석을 적용하면 이 "닫힘 제약
(closure constraint)" 때문에 카테고리 간에 인위적인(spurious) 음의 상관관계가
생길 수 있다(Aitchison, 1986). 표준적인 해결책은 PCA 전에 CLR(centered
log-ratio) 변환을 적용하는 것이다. 이 스크립트는 기존 raw-percentage PCA와
CLR-변환 PCA를 나란히 비교하여, LAB vs Bacillus-group 분리 패턴이 변환 후에도
유지되는지 확인한다.

사용법:
    python3 cog_pca_clr.py \
        --cog-ratio eggnog_cog_ratio_wide.tsv \
        --outdir cog_pca_clr_result
"""
import argparse
import os
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def clr_transform(df, pseudocount=None):
    """centered log-ratio 변환. 0값은 관례대로 작은 pseudocount로 대체."""
    X = df.values.astype(float)
    if pseudocount is None:
        # 0이 아닌 최솟값의 절반을 pseudocount로 사용 (표준적 관행)
        nonzero_min = X[X > 0].min()
        pseudocount = nonzero_min / 2
    X = np.where(X == 0, pseudocount, X)
    log_X = np.log(X)
    clr = log_X - log_X.mean(axis=1, keepdims=True)
    return clr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cog-ratio", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    df = pd.read_csv(args.cog_ratio, sep="\t").set_index("sample_id")
    cog_cols = [c for c in df.columns if c not in
                {"genus_final", "species_final", "functional_group"} and len(c) <= 2]

    cog_data = df[cog_cols].fillna(0)
    fg = df["functional_group"] if "functional_group" in df.columns else None

    # 1) 기존 방식: raw percentage PCA
    pca_raw = PCA(n_components=2)
    coords_raw = pca_raw.fit_transform(cog_data.values)
    var_raw = pca_raw.explained_variance_ratio_ * 100

    # 2) CLR 변환 PCA
    clr_data = clr_transform(cog_data)
    pca_clr = PCA(n_components=2)
    coords_clr = pca_clr.fit_transform(clr_data)
    var_clr = pca_clr.explained_variance_ratio_ * 100

    print(f"Raw percentage PCA: PC1 {var_raw[0]:.1f}%, PC2 {var_raw[1]:.1f}% 설명")
    print(f"CLR-transformed PCA: PC1 {var_clr[0]:.1f}%, PC2 {var_clr[1]:.1f}% 설명")

    # 그림 비교
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    colors = {"LAB": "tab:orange", "Bacillus_group": "tab:blue",
              "Other_Environmental": "tab:green", "Unresolved": "tab:red"}
    for ax, coords, var, title in [
        (axes[0], coords_raw, var_raw, "Raw percentage PCA (original)"),
        (axes[1], coords_clr, var_clr, "CLR-transformed PCA (compositional-data-appropriate)"),
    ]:
        for g, c in colors.items():
            mask = (fg == g).values if fg is not None else np.ones(len(coords), dtype=bool)
            ax.scatter(coords[mask, 0], coords[mask, 1], label=g, alpha=0.7, color=c, s=25)
        ax.set_xlabel(f"PC1 ({var[0]:.1f}%)")
        ax.set_ylabel(f"PC2 ({var[1]:.1f}%)")
        ax.set_title(title)
        ax.legend(fontsize=8)
    plt.tight_layout()
    out_fig = os.path.join(args.outdir, "cog_pca_raw_vs_clr.pdf")
    plt.savefig(out_fig, dpi=150)
    print(f"비교 그림 저장: {out_fig}")

    # LAB vs Bacillus_group PC1 분리가 유지되는지 정량 확인 (Mann-Whitney)
    if fg is not None:
        from scipy import stats
        lab_mask = (fg == "LAB").values
        bac_mask = (fg == "Bacillus_group").values
        for coords, label in [(coords_raw, "Raw"), (coords_clr, "CLR")]:
            u, p = stats.mannwhitneyu(coords[lab_mask, 0], coords[bac_mask, 0], alternative="two-sided")
            print(f"{label} PC1: LAB vs Bacillus_group Mann-Whitney p = {p:.4g}")

    out_csv = os.path.join(args.outdir, "pca_comparison_summary.csv")
    pd.DataFrame({
        "method": ["raw_percentage", "CLR_transformed"],
        "PC1_variance_pct": [var_raw[0], var_clr[0]],
        "PC2_variance_pct": [var_raw[1], var_clr[1]],
    }).to_csv(out_csv, index=False)
    print(f"요약 저장: {out_csv}")


if __name__ == "__main__":
    main()
