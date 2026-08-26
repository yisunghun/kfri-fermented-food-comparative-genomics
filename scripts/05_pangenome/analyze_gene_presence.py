#!/usr/bin/env python3
"""
analyze_gene_presence.py

각 그룹의 gene_presence_absence.Rtab을 읽어 유전자 패밀리별 '존재 비율'
(해당 유전자가 전체 시료 중 몇 %에 있는지)을 직접 계산합니다.
Panaroo의 고정 99% core 기준 대신, 여러 임계값에서 몇 개 유전자가
core-like 한지 유연하게 보여주고, 그룹별 분포를 히스토그램으로 시각화합니다.

사용법:
  python3 analyze_gene_presence.py \
      --pangenome-root /mnt/f/WGS_Consolidated/pangenome \
      --outdir /mnt/f/WGS_Consolidated/pangenome/presence_analysis
"""
import argparse
import os

import numpy as np
import pandas as pd

THRESHOLDS = [100, 99, 95, 90, 80, 50]


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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pangenome-root", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    globals()["plt"] = plt
    setup_korean_font()

    os.makedirs(args.outdir, exist_ok=True)

    group_dirs = []
    for root, dirs, files in os.walk(args.pangenome_root):
        if "gene_presence_absence.Rtab" in files:
            rel = os.path.relpath(root, args.pangenome_root)
            parts = rel.split(os.sep)
            if len(parts) >= 2:
                group_dirs.append((parts[0], parts[1], root))

    summary_rows = []
    n_groups = len(group_dirs)
    ncols = 4
    nrows = int(np.ceil(n_groups / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 3.5 * nrows))
    axes = np.array(axes).reshape(-1)

    for i, (category, group, root) in enumerate(sorted(group_dirs)):
        rtab_path = os.path.join(root, "gene_presence_absence.Rtab")
        df = pd.read_csv(rtab_path, sep="\t", index_col=0)
        n_genomes = df.shape[1]
        pct_present = (df.sum(axis=1) / n_genomes * 100)

        row = {"category": category, "group": group, "n_genomes": n_genomes,
               "n_gene_families": len(pct_present)}
        for th in THRESHOLDS:
            row[f"n_genes_present_>={th}pct"] = int((pct_present >= th).sum())
        summary_rows.append(row)

        ax = axes[i]
        ax.hist(pct_present, bins=20, range=(0, 100), color="steelblue", edgecolor="white")
        ax.set_title(f"{category}/{group} (n={n_genomes})", fontsize=9)
        ax.set_xlabel("% of genomes containing gene", fontsize=8)
        ax.set_ylabel("Number of gene families", fontsize=8)
        ax.tick_params(labelsize=7)

    for j in range(len(group_dirs), len(axes)):
        axes[j].axis("off")

    fig.suptitle("Gene presence-frequency distribution by group "
                 "(U-shaped = distinct core/accessory; flat/left-skewed = open pangenome)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig_path = os.path.join(args.outdir, "presence_histograms.pdf")
    fig.savefig(fig_path)
    print(f"히스토그램 저장: {fig_path}")

    summary_df = pd.DataFrame(summary_rows).sort_values(["category", "n_genomes"], ascending=[True, False])
    summary_csv = os.path.join(args.outdir, "presence_threshold_summary.csv")
    summary_df.to_csv(summary_csv, index=False, encoding="utf-8-sig")
    print(f"임계값별 요약 저장: {summary_csv}\n")
    print(summary_df.to_string(index=False))


if __name__ == "__main__":
    main()
