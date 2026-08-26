#!/usr/bin/env python3
"""
genus_cog_heatmap.py

eggnog_cog_ratio_wide.tsv를 이용해 속(genus_final) 단위로 COG 카테고리
평균 프로파일 히트맵을 그리고, 카테고리별 Kruskal-Wallis 검정으로
"어떤 기능이 속 사이에서 가장 크게 갈리는지"를 정리합니다.

산출물:
  1) genus_cog_heatmap.pdf - 속 x COG카테고리 히트맵 (계층적 클러스터링으로 속 재정렬)
  2) genus_cog_means.csv - 속별 COG카테고리 평균 비율 표
  3) genus_cog_kruskal.csv - 카테고리별 Kruskal-Wallis H-검정 + BH-FDR

사용법:
  python3 genus_cog_heatmap.py \
      --ratio-tsv eggnog_cog_ratio_wide.tsv \
      --outdir ./genus_cog_analysis \
      --min-n 3
"""
import argparse
import os

import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import dendrogram, leaves_list, linkage

NON_CATEGORY_COLS = {"Total", "-", "---------------------------------------------------------"}
META_COLS = {"sample_id", "genus_final", "species_final", "functional_group"}

COG_ORDER = list("JAKLBDYVT MNZWUOCGEFHIPQRS".replace(" ", ""))  # 기능적으로 대략 묶어 정렬


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
    ap.add_argument("--min-n", type=int, default=3)
    args = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    globals()["plt"] = plt
    setup_korean_font()

    os.makedirs(args.outdir, exist_ok=True)

    df = pd.read_csv(args.ratio_tsv, sep="\t", encoding="utf-8-sig")
    cog_cols_all = [c for c in df.columns if c not in META_COLS and c not in NON_CATEGORY_COLS]
    cog_cols = [c for c in COG_ORDER if c in cog_cols_all] + \
               [c for c in cog_cols_all if c not in COG_ORDER]

    genus_counts = df["genus_final"].value_counts()
    valid_genera = genus_counts[genus_counts >= args.min_n].index.tolist()
    valid_genera = [g for g in valid_genera if g != "unresolved"]
    sub = df[df["genus_final"].isin(valid_genera)].copy()
    print(f"분석 대상 속: {len(valid_genera)}개 (n>={args.min_n}), 총 {len(sub)}개 시료")

    # ---------- 1) 속별 평균 + 히트맵 ----------
    means = sub.groupby("genus_final")[cog_cols].mean()
    means_csv = os.path.join(args.outdir, "genus_cog_means.csv")
    means.to_csv(means_csv, encoding="utf-8-sig")
    print(f"속별 평균 저장: {means_csv}")

    # 행(속)을 유사도 기준으로 재정렬 (z-score 후 클러스터링)
    z = (means - means.mean()) / means.std(ddof=0).replace(0, 1)
    row_linkage = linkage(z.fillna(0), method="average", metric="euclidean")
    row_order = leaves_list(row_linkage)
    means_ordered = means.iloc[row_order]

    fig, ax = plt.subplots(figsize=(max(10, len(cog_cols) * 0.5), max(6, len(means_ordered) * 0.4)))
    im = ax.imshow(means_ordered.to_numpy(), aspect="auto", cmap="YlOrRd")
    ax.set_xticks(range(len(cog_cols)))
    ax.set_xticklabels(cog_cols, fontsize=9)
    ax.set_yticks(range(len(means_ordered)))
    ax.set_yticklabels(means_ordered.index, fontsize=9)
    ax.set_title("Mean COG Category Relative Abundance (%) by Genus")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Mean relative abundance (%)")

    # 셀에 값 표기
    for i in range(means_ordered.shape[0]):
        for j in range(means_ordered.shape[1]):
            val = means_ordered.iat[i, j]
            ax.text(j, i, f"{val:.1f}", ha="center", va="center", fontsize=6,
                     color="black" if val < means_ordered.to_numpy().max() * 0.6 else "white")

    fig.tight_layout()
    heatmap_path = os.path.join(args.outdir, "genus_cog_heatmap.pdf")
    fig.savefig(heatmap_path)
    print(f"히트맵 저장: {heatmap_path}")

    # ---------- 2) 카테고리별 Kruskal-Wallis (여러 속 간 전체 비교) ----------
    rows = []
    for c in cog_cols:
        groups = [sub.loc[sub["genus_final"] == g, c].to_numpy() for g in valid_genera]
        try:
            h, p = stats.kruskal(*groups)
        except ValueError:
            h, p = np.nan, np.nan
        rows.append({"cog_category": c, "H_stat": h, "p_value": p})
    kruskal_df = pd.DataFrame(rows)
    kruskal_df["p_adj_BH"] = bh_fdr(kruskal_df["p_value"].to_numpy())
    kruskal_df["significant_FDR<0.05"] = kruskal_df["p_adj_BH"] < 0.05
    kruskal_df = kruskal_df.sort_values("p_adj_BH")
    kruskal_csv = os.path.join(args.outdir, "genus_cog_kruskal.csv")
    kruskal_df.to_csv(kruskal_csv, index=False, encoding="utf-8-sig")
    print(f"Kruskal-Wallis 결과 저장: {kruskal_csv}")
    print("\n[속들 간 가장 크게 갈리는 카테고리 상위 10개]")
    print(kruskal_df.head(10).to_string(index=False))


if __name__ == "__main__":
    main()
