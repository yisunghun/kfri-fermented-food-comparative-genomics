#!/usr/bin/env python3
"""
global_tree_heatmap.py (Fig. 5 단순화 개정판)

투고 전 자체 검토에서, 본문 Fig. 5에 220개 isolate label이 전부 들어가서
본문 크기에서는 사실상 해석 불가능하다는 지적을 받아 개조했다.

- 기본(라벨 숨김) 모드: guild 색상 tick만 남기고 텍스트 라벨은 전부 제거,
  --highlight-sample로 지정한 시료(K. pneumoniae 이상치)만 굵게 표시.
  -> 본문 Fig. 5용.
- --show-all-labels 플래그를 주면 기존처럼 220개 라벨을 전부 표시.
  -> Supplementary Fig. S1(4-panel 분할)용, 기존 동작과 동일.

Usage (본문 Fig. 5, 단순화):
  python3 global_tree_heatmap.py \
      --ani-matrix /mnt/f/WGS_Consolidated/ani_analysis/ani_matrix.csv \
      --master /mnt/f/WGS_Consolidated/master_table_qc.tsv \
      --antismash-summary /mnt/f/WGS_Consolidated/antismash_analysis/antismash_bgc_summary.tsv \
      --card-summary /mnt/f/WGS_Consolidated/abricate_out/card_summary.tsv \
      --vfdb-summary /mnt/f/WGS_Consolidated/abricate_out/vfdb_summary.tsv \
      --highlight-sample HN00179262_F4055 \
      --outdir /mnt/f/WGS_Consolidated/global_summary

Usage (Supplementary Fig. S1, 전체 라벨):
  (위 명령에 --show-all-labels 추가, --highlight-sample 생략 가능)
"""
import argparse
import os

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform

LAB_GENERA = {
    "Lactiplantibacillus", "Levilactobacillus", "Latilactobacillus", "Lactobacillus",
    "Lactococcus", "Leuconostoc", "Weissella", "Pediococcus", "Enterococcus",
    "Tetragenococcus", "Lacticaseibacillus", "Limosilactobacillus", "Lentilactobacillus",
    "Loigolactobacillus", "Fructilactobacillus",
}
BACILLUS_GROUP_GENERA = {
    "Bacillus", "Paenibacillus", "Oceanobacillus", "Virgibacillus", "Priestia",
    "Rossellomorea", "Shouchella", "Halobacillus",
}
GROUP_COLORS = {"LAB": "#E69F00", "Bacillus_group": "#0072B2",
                "Other_Environmental": "#009E73", "Unresolved": "#D55E00"}


def functional_group(genus: str) -> str:
    if genus == "unresolved":
        return "Unresolved"
    if genus in LAB_GENERA:
        return "LAB"
    if genus in BACILLUS_GROUP_GENERA:
        return "Bacillus_group"
    return "Other_Environmental"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ani-matrix", required=True)
    ap.add_argument("--master", required=True)
    ap.add_argument("--antismash-summary", required=False)
    ap.add_argument("--card-summary", required=False)
    ap.add_argument("--vfdb-summary", required=False)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--show-all-labels", action="store_true",
                     help="220개 시료 라벨을 전부 표시 (Supplementary Fig. S1용). "
                          "기본값(False)이면 본문 Fig. 5용으로 라벨을 숨긴다.")
    ap.add_argument("--highlight-sample", default=None,
                     help="라벨을 숨기는 모드에서 예외적으로 굵게 강조 표시할 "
                          "단일 시료 (예: K. pneumoniae 이상치 sample_id)")
    args = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    os.makedirs(args.outdir, exist_ok=True)

    ani = pd.read_csv(args.ani_matrix, index_col=0)
    master = pd.read_csv(args.master, sep="\t").set_index("sample_id")
    master["functional_group"] = master["genus_final"].apply(functional_group)

    dist_arr = (100.0 - ani.fillna(70.0)).to_numpy(copy=True)
    dist_arr = np.clip(dist_arr, 0.0, None)
    np.fill_diagonal(dist_arr, 0.0)
    condensed = squareform(dist_arr, checks=False)
    Z = linkage(condensed, method="average")

    n = len(ani.index)
    fig_height = max(8, n * 0.16) if args.show_all_labels else 11.0
    fig = plt.figure(figsize=(9, fig_height), constrained_layout=True)
    gs = fig.add_gridspec(1, 2, width_ratios=[2.0, 1.0])
    ax_tree = fig.add_subplot(gs[0, 0])
    ax_heat = fig.add_subplot(gs[0, 1])

    dendro = dendrogram(Z, labels=ani.index.tolist(), orientation="left", ax=ax_tree,
                         leaf_font_size=7, no_labels=not args.show_all_labels)
    leaf_order = dendro["ivl"]
    ax_tree.set_title("ANI-based hierarchical clustering (all 220 genomes;\nnot a phylogenetic tree)")
    ax_tree.set_xlabel("Distance (100 - ANI)")

    fig.get_layout_engine().set(wspace=0.03, w_pad=0.01)

    if args.show_all_labels:
        old_labels = [t.get_text() for t in ax_tree.get_ymajorticklabels()]
        short_labels = []
        for sid in old_labels:
            if sid in master.index and "short_id" in master.columns and pd.notna(master.loc[sid, "short_id"]):
                short_labels.append(str(master.loc[sid, "short_id"]))
            else:
                short_labels.append(sid)
        ax_tree.set_yticklabels(short_labels, fontsize=7)
        for old_sid, tick_label in zip(old_labels, ax_tree.get_ymajorticklabels()):
            fg = master.loc[old_sid, "functional_group"] if old_sid in master.index else "Unresolved"
            tick_label.set_color(GROUP_COLORS.get(fg, "black"))
    else:
        ax_tree.set_yticks(np.arange(5, n * 10 + 5, 10))
        tick_colors = [GROUP_COLORS.get(
            master.loc[s, "functional_group"] if s in master.index else "Unresolved", "black"
        ) for s in leaf_order]
        ax_tree.set_yticklabels([""] * n)
        for tick, color in zip(ax_tree.yaxis.get_major_ticks(), tick_colors):
            tick.tick1line.set_markersize(4)
            tick.tick1line.set_markeredgewidth(1.5)
            tick.tick1line.set_color(color)

        if args.highlight_sample and args.highlight_sample in leaf_order:
            idx = leaf_order.index(args.highlight_sample)
            y_pos = idx * 10 + 5
            display_name = (master.loc[args.highlight_sample, "short_id"]
                             if args.highlight_sample in master.index
                             and "short_id" in master.columns
                             and pd.notna(master.loc[args.highlight_sample, "short_id"])
                             else args.highlight_sample)
            ax_tree.annotate(
                f"{display_name}\n(K. pneumoniae outlier)",
                xy=(0, y_pos), xytext=(15, y_pos),
                textcoords="data", fontsize=8, fontweight="bold", color="firebrick",
                va="center", ha="left", clip_on=False, annotation_clip=False,
                arrowprops=dict(arrowstyle="-", color="firebrick", lw=1.2),
            )

    tracks = pd.DataFrame(index=leaf_order)
    tracks["functional_group_code"] = [
        list(GROUP_COLORS.keys()).index(master.loc[s, "functional_group"])
        if s in master.index else np.nan for s in leaf_order
    ]

    if args.antismash_summary and os.path.isfile(args.antismash_summary):
        asum = pd.read_csv(args.antismash_summary, sep="\t", encoding="utf-8-sig").set_index("sample_id")
        tracks["antiSMASH_BGC_count"] = asum.reindex(leaf_order)["total_bgc_count"]

    for label, path in [("CARD_gene_count", args.card_summary), ("VFDB_gene_count", args.vfdb_summary)]:
        if path and os.path.isfile(path):
            s = pd.read_csv(path, sep="\t")
            s["sample_id"] = s["#FILE"].apply(lambda p: os.path.basename(str(p)).replace(".tab", ""))
            s = s.set_index("sample_id")
            tracks[label] = s.reindex(leaf_order)["NUM_FOUND"]

    numeric_cols = [c for c in tracks.columns if c != "functional_group_code"]
    if numeric_cols:
        norm_data = tracks[numeric_cols].apply(lambda col: (col - col.min()) / (col.max() - col.min() + 1e-9))
        im = ax_heat.imshow(norm_data.to_numpy(dtype=float), aspect="auto", cmap="YlOrRd",
                             extent=[1, len(numeric_cols) + 1, n + 0.5, 0.5])
        ax_heat.set_xticks(np.arange(len(numeric_cols)) + 1.5)
        ax_heat.set_xticklabels(numeric_cols, rotation=45, ha="right", fontsize=8)
        fig.colorbar(im, ax=ax_heat, fraction=0.03, pad=0.02, label="Normalized value (0-1 per column)")

        if not args.show_all_labels and args.highlight_sample and args.highlight_sample in leaf_order:
            idx = leaf_order.index(args.highlight_sample)
            ax_heat.annotate("", xy=(len(numeric_cols) + 1.3, idx + 1), xytext=(len(numeric_cols) + 2.3, idx + 1),
                              annotation_clip=False,
                              arrowprops=dict(arrowstyle="->", color="firebrick", lw=1.5))
    ax_heat.set_yticks([])
    ax_heat.set_title("Screening summary\n(CARD/VFDB gene counts, antiSMASH BGC count)")

    legend_elems = [Patch(facecolor=c, label=g) for g, c in GROUP_COLORS.items()]
    fig.legend(handles=legend_elems, loc="outside upper right", title="Functional group", fontsize=9)

    suffix = "_full_labels" if args.show_all_labels else "_simplified"
    out_path = os.path.join(args.outdir, f"global_tree_screening_heatmap{suffix}.pdf")
    fig.savefig(out_path)
    print(f"Combined figure saved: {out_path}")

    tracks_path = os.path.join(args.outdir, f"global_tracks_table{suffix}.csv")
    tracks.to_csv(tracks_path, encoding="utf-8-sig")
    print(f"Underlying track data saved: {tracks_path}")


if __name__ == "__main__":
    main()