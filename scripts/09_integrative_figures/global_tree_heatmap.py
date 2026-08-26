#!/usr/bin/env python3
"""
global_tree_heatmap.py

Builds the final integrative figure: ANI-based hierarchical clustering tree
for all genomes, with an adjacent multi-track heatmap showing:
  - functional_group (categorical: LAB / Bacillus_group / Other_Environmental / Unresolved)
  - antiSMASH total BGC count
  - CARD gene count (antibiotic resistance burden)
  - VFDB gene count (virulence factor burden)

Usage:
  python3 global_tree_heatmap.py \
      --ani-matrix /mnt/f/WGS_Consolidated/ani_analysis/ani_matrix.csv \
      --master /mnt/f/WGS_Consolidated/master_table_qc.tsv \
      --antismash-summary /mnt/f/WGS_Consolidated/antismash_analysis/antismash_bgc_summary.tsv \
      --card-summary /mnt/f/WGS_Consolidated/abricate_out/card_summary.tsv \
      --vfdb-summary /mnt/f/WGS_Consolidated/abricate_out/vfdb_summary.tsv \
      --outdir /mnt/f/WGS_Consolidated/global_summary
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
    args = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    os.makedirs(args.outdir, exist_ok=True)

    ani = pd.read_csv(args.ani_matrix, index_col=0)
    master = pd.read_csv(args.master, sep="\t").set_index("sample_id")
    master["functional_group"] = master["genus_final"].apply(functional_group)

    # ---------- build dendrogram (reuse ANI-distance approach) ----------
    dist_arr = (100.0 - ani.fillna(70.0)).to_numpy(copy=True)
    dist_arr = np.clip(dist_arr, 0.0, None)
    np.fill_diagonal(dist_arr, 0.0)
    condensed = squareform(dist_arr, checks=False)
    Z = linkage(condensed, method="average")

    n = len(ani.index)
    fig = plt.figure(figsize=(18, max(10, n * 0.13)))
    gs = fig.add_gridspec(1, 2, width_ratios=[2.5, 1.5], wspace=0.02)
    ax_tree = fig.add_subplot(gs[0, 0])
    ax_heat = fig.add_subplot(gs[0, 1])

    dendro = dendrogram(Z, labels=ani.index.tolist(), orientation="left", ax=ax_tree,
                         leaf_font_size=5, no_labels=False)
    leaf_order = dendro["ivl"]
    ax_tree.set_title("ANI-based hierarchical clustering (all genomes)")
    ax_tree.set_xlabel("Distance (100 - ANI)")

    # colour tip labels by functional group
    for tick_label in ax_tree.get_ymajorticklabels():
        sid = tick_label.get_text()
        fg = master.loc[sid, "functional_group"] if sid in master.index else "Unresolved"
        tick_label.set_color(GROUP_COLORS.get(fg, "black"))

    # ---------- assemble numeric tracks ----------
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

    # ---------- plot: categorical strip for functional_group + numeric heatmap for the rest ----------
    numeric_cols = [c for c in tracks.columns if c != "functional_group_code"]
    if numeric_cols:
        norm_data = tracks[numeric_cols].apply(lambda col: (col - col.min()) / (col.max() - col.min() + 1e-9))
        im = ax_heat.imshow(norm_data.to_numpy(dtype=float), aspect="auto", cmap="YlOrRd",
                             extent=[1, len(numeric_cols) + 1, n + 0.5, 0.5])
        ax_heat.set_xticks(np.arange(len(numeric_cols)) + 1.5)
        ax_heat.set_xticklabels(numeric_cols, rotation=45, ha="right", fontsize=8)
        fig.colorbar(im, ax=ax_heat, fraction=0.03, pad=0.02, label="Normalized value (0-1 per column)")
    ax_heat.set_yticks([])
    ax_heat.set_title("Screening summary")

    legend_elems = [Patch(facecolor=c, label=g) for g, c in GROUP_COLORS.items()]
    fig.legend(handles=legend_elems, loc="upper right", title="Functional group", fontsize=9)

    fig.tight_layout()
    out_path = os.path.join(args.outdir, "global_tree_screening_heatmap.pdf")
    fig.savefig(out_path)
    print(f"Combined figure saved: {out_path}")

    tracks_path = os.path.join(args.outdir, "global_tracks_table.csv")
    tracks.to_csv(tracks_path, encoding="utf-8-sig")
    print(f"Underlying track data saved: {tracks_path}")


if __name__ == "__main__":
    main()
