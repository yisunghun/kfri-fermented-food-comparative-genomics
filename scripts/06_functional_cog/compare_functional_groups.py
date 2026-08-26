#!/usr/bin/env python3
"""
compare_functional_groups.py

eggnog_cog_ratio_wide.tsv를 이용해 기능군(functional_group, 기본 LAB vs
Bacillus_group)간 COG 카테고리 구성을 비교합니다.

산출물:
  1) cog_comparison_barplot.pdf - 카테고리별 그룹 평균 비율 막대그래프
  2) cog_group_stats.csv - 카테고리별 Mann-Whitney U 검정 + BH-FDR 보정 p-value
  3) cog_pca.pdf - 시료들을 COG 프로파일 기준 PCA로 투영, 그룹별 색상

사용법:
  python3 compare_functional_groups.py \
      --ratio-tsv eggnog_cog_ratio_wide.tsv \
      --outdir ./functional_comparison \
      --group-a LAB --group-b Bacillus_group
"""
import argparse
import os

import numpy as np
import pandas as pd
from scipy import stats

NON_CATEGORY_COLS = {"Total", "-", "---------------------------------------------------------"}
META_COLS = {"sample_id", "genus_final", "species_final", "functional_group"}

COG_DESC = {
    "J": "Translation/ribosome", "A": "RNA processing", "K": "Transcription",
    "L": "Replication/repair", "B": "Chromatin structure", "D": "Cell cycle/division",
    "Y": "Nuclear structure", "V": "Defense mechanisms", "T": "Signal transduction",
    "M": "Cell wall/membrane", "N": "Cell motility", "Z": "Cytoskeleton",
    "W": "Extracellular structures", "U": "Intracellular trafficking/secretion",
    "O": "Posttranslational modification/chaperones", "C": "Energy production",
    "G": "Carbohydrate metabolism", "E": "Amino acid metabolism",
    "F": "Nucleotide metabolism", "H": "Coenzyme metabolism", "I": "Lipid metabolism",
    "P": "Inorganic ion metabolism", "Q": "Secondary metabolites",
    "R": "General function prediction", "S": "Function unknown",
}


def setup_korean_font():
    import matplotlib.font_manager as fm
    candidates = ["NanumGothic", "Nanum Gothic", "Malgun Gothic", "AppleGothic",
                  "Noto Sans CJK KR", "Noto Sans KR"]
    available = {f.name for f in fm.fontManager.ttflist}
    for c in candidates:
        if c in available:
            plt.rcParams["font.family"] = c
            plt.rcParams["axes.unicode_minus"] = False
            return c
    return None


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR 보정 (statsmodels 없이 직접 구현).
    NaN p-value(분산 0인 카테고리 등)는 검정 불가로 보고 1.0(비유의)로 처리."""
    pvals = np.where(np.isnan(pvals), 1.0, pvals)
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order] * n / (np.arange(n) + 1)
    ranked_min = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked_min, 0, 1)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ratio-tsv", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--group-a", default="LAB")
    ap.add_argument("--group-b", default="Bacillus_group")
    args = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    globals()["plt"] = plt
    setup_korean_font()

    os.makedirs(args.outdir, exist_ok=True)

    df = pd.read_csv(args.ratio_tsv, sep="\t", encoding="utf-8-sig")
    cog_cols = [c for c in df.columns if c not in META_COLS and c not in NON_CATEGORY_COLS]
    print(f"COG 카테고리 컬럼: {cog_cols}")

    sub = df[df["functional_group"].isin([args.group_a, args.group_b])].copy()
    print(f"{args.group_a}: {(sub['functional_group']==args.group_a).sum()}개, "
          f"{args.group_b}: {(sub['functional_group']==args.group_b).sum()}개")

    # ---------- 1) 막대그래프 ----------
    means = sub.groupby("functional_group")[cog_cols].mean().T
    means = means.sort_values(by=args.group_a, ascending=False)

    fig, ax = plt.subplots(figsize=(14, 6))
    x = np.arange(len(means))
    width = 0.35
    ax.bar(x - width / 2, means[args.group_a], width, label=args.group_a)
    ax.bar(x + width / 2, means[args.group_b], width, label=args.group_b)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{c}\n({COG_DESC.get(c, '')})" for c in means.index],
                        rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Mean relative abundance (%)")
    ax.set_title(f"{args.group_a} vs {args.group_b} COG Category Composition")
    ax.legend()
    fig.tight_layout()
    bar_path = os.path.join(args.outdir, "cog_comparison_barplot.pdf")
    fig.savefig(bar_path)
    print(f"막대그래프 저장: {bar_path}")

    # ---------- 2) 통계검정 (Mann-Whitney U + BH-FDR) ----------
    rows = []
    a_vals = sub[sub["functional_group"] == args.group_a]
    b_vals = sub[sub["functional_group"] == args.group_b]
    for c in cog_cols:
        try:
            u, p = stats.mannwhitneyu(a_vals[c], b_vals[c], alternative="two-sided")
        except ValueError:
            # 두 그룹 모두 분산이 0(완전 동일값)인 경우 등 검정 자체가 성립하지 않음
            u, p = np.nan, np.nan
        rows.append({
            "cog_category": c, "description": COG_DESC.get(c, ""),
            f"mean_{args.group_a}": round(a_vals[c].mean(), 3),
            f"mean_{args.group_b}": round(b_vals[c].mean(), 3),
            "diff": round(a_vals[c].mean() - b_vals[c].mean(), 3),
            "U_stat": u, "p_value": p,
        })
    stats_df = pd.DataFrame(rows)
    stats_df["p_adj_BH"] = bh_fdr(stats_df["p_value"].to_numpy())
    stats_df["significant_FDR<0.05"] = stats_df["p_adj_BH"] < 0.05
    stats_df = stats_df.sort_values("p_adj_BH")
    stats_csv = os.path.join(args.outdir, "cog_group_stats.csv")
    stats_df.to_csv(stats_csv, index=False, encoding="utf-8-sig")
    print(f"통계검정 결과 저장: {stats_csv}")
    print(f"\n[유의미하게 차이나는 카테고리 (FDR<0.05)]")
    sig = stats_df[stats_df["significant_FDR<0.05"]]
    print(sig.to_string(index=False) if not sig.empty else "  없음")

    # ---------- 3) PCA ----------
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler

    X = df[cog_cols].to_numpy()
    X_scaled = StandardScaler().fit_transform(X)
    pca = PCA(n_components=2)
    coords = pca.fit_transform(X_scaled)

    fig2, ax2 = plt.subplots(figsize=(9, 7))
    groups = df["functional_group"].unique()
    cmap = plt.colormaps["tab10"]
    for i, g in enumerate(sorted(groups)):
        mask = df["functional_group"] == g
        ax2.scatter(coords[mask, 0], coords[mask, 1], label=g, alpha=0.7, color=cmap(i))
    ax2.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
    ax2.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
    ax2.set_title("PCA of COG Functional Profiles (all groups)")
    ax2.legend()
    fig2.tight_layout()
    pca_path = os.path.join(args.outdir, "cog_pca.pdf")
    fig2.savefig(pca_path)
    print(f"PCA 플롯 저장: {pca_path}")


if __name__ == "__main__":
    main()
