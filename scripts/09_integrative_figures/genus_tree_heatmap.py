#!/usr/bin/env python3
"""
genus_tree_heatmap.py

For a given genus, builds a phylogenetic tree from Panaroo's
core_gene_alignment.aln (via FastTree) and plots it side-by-side with a
COG functional-profile heatmap for the same genomes.

Requires: fasttree (binary in PATH), biopython

Usage:
  python3 genus_tree_heatmap.py \
      --genus Bacillus \
      --pangenome-root /mnt/f/WGS_Consolidated/pangenome \
      --ratio-tsv /mnt/f/WGS_Consolidated/eggnog_cog_ratio_wide.tsv \
      --master /mnt/f/WGS_Consolidated/master_table_qc.tsv \
      --outdir /mnt/f/WGS_Consolidated/tree_heatmap
"""
import argparse
import os
import subprocess

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform

NON_CATEGORY_COLS = {"Total", "-", "---------------------------------------------------------"}
META_COLS = {"sample_id", "genus_final", "species_final", "functional_group"}
COG_ORDER = list("JAKLBDYVTMNZWUOCGEFHIPQRS")


def build_ani_fallback_tree(ani_matrix_path, sample_ids, outdir, genus):
    """core_gene_alignment.aln이 없을 때(core gene 0개) ANI 매트릭스 기반 계층적
    클러스터링으로 대체 트리를 만든다. 반환: (leaf_order, Z_linkage) - Phylo 트리 대신
    scipy dendrogram으로 그림."""
    ani = pd.read_csv(ani_matrix_path, index_col=0)
    sub = ani.reindex(index=sample_ids, columns=sample_ids)
    dist_arr = (100.0 - sub.fillna(70.0)).to_numpy(copy=True)
    dist_arr = np.clip(dist_arr, 0.0, None)
    np.fill_diagonal(dist_arr, 0.0)
    condensed = squareform(dist_arr, checks=False)
    Z = linkage(condensed, method="average")
    return Z, sub.index.tolist()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genus", required=True)
    ap.add_argument("--pangenome-root", required=True)
    ap.add_argument("--ratio-tsv", required=True)
    ap.add_argument("--master", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--fasttree-bin", default="fasttree")
    ap.add_argument("--ani-matrix", default=None,
                     help="core_gene_alignment.aln이 없는 속(core gene 0개)일 때 "
                          "대신 사용할 ANI 매트릭스 csv 경로 (analyze_ani.py 산출물)")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    aln_path = os.path.join(args.pangenome_root, "by_genus", args.genus, "core_gene_alignment.aln")
    use_fallback = not (os.path.isfile(aln_path) and os.path.getsize(aln_path) > 0)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    master = pd.read_csv(args.master, sep="\t").set_index("sample_id")

    if use_fallback:
        if not args.ani_matrix:
            raise SystemExit(
                f"ERROR: core_gene_alignment.aln missing/empty for genus '{args.genus}' "
                f"(likely 0 core genes). Provide --ani-matrix to use an ANI-based tree instead."
            )
        print(f"[INFO] No core-genome alignment for '{args.genus}' (0 core genes) - "
              f"falling back to ANI-based clustering.")
        genus_samples = master[master["genus_final"] == args.genus].index.tolist()
        Z, leaf_order = build_ani_fallback_tree(args.ani_matrix, genus_samples, args.outdir, args.genus)
        n = len(leaf_order)

        fig = plt.figure(figsize=(16, max(6, n * 0.28)))
        gs = fig.add_gridspec(1, 2, width_ratios=[2.2, 1], wspace=0.05)
        ax_tree = fig.add_subplot(gs[0, 0])
        ax_heat = fig.add_subplot(gs[0, 1])

        def leaf_label(sid):
            sp = master.loc[sid, "species_final"] if sid in master.index else "?"
            return f"{sid} ({sp})"

        dendro = dendrogram(Z, labels=[leaf_label(s) for s in leaf_order], orientation="left",
                             ax=ax_tree, leaf_font_size=7)
        # dendrogram의 leaf 순서(ivl)는 원래 sample_id가 아니라 leaf_label 문자열이므로,
        # 히트맵 정렬을 위해 dendro['leaves']의 인덱스로 원래 순서를 복원
        leaf_order = [leaf_order[i] for i in dendro["leaves"]]
        ax_tree.set_title(f"ANI-based clustering: {args.genus} (no core genes for alignment)")
        ax_tree.set_xlabel("Distance (100 - ANI)")
    else:
        from Bio import Phylo

        nwk_path = os.path.join(args.outdir, f"{args.genus}.nwk")
        print(f"Building tree with FastTree from {aln_path} ...")
        with open(nwk_path, "w") as out_f:
            subprocess.run([args.fasttree_bin, "-nt", aln_path], stdout=out_f, check=True,
                            stderr=subprocess.PIPE)
        print(f"Tree saved: {nwk_path}")

        tree = Phylo.read(nwk_path, "newick")
        tree.ladderize()
        leaf_order = [t.name for t in tree.get_terminals()]
        n = len(leaf_order)
        print(f"Number of leaves: {n}")

        def label_func(clade):
            if clade.name and clade.name in master.index:
                sp = master.loc[clade.name, "species_final"]
                return f"{clade.name} ({sp})"
            return clade.name

        fig = plt.figure(figsize=(16, max(6, n * 0.28)))
        gs = fig.add_gridspec(1, 2, width_ratios=[2.2, 1], wspace=0.05)
        ax_tree = fig.add_subplot(gs[0, 0])
        ax_heat = fig.add_subplot(gs[0, 1])

        Phylo.draw(tree, do_show=False, axes=ax_tree, label_func=label_func, show_confidence=False)
        ax_tree.set_title(f"Core-genome phylogeny: {args.genus}")
        for side in ["top", "right", "left"]:
            ax_tree.spines[side].set_visible(False)

    # ---------- load COG ratio data for these samples ----------
    df = pd.read_csv(args.ratio_tsv, sep="\t", encoding="utf-8-sig")
    cog_cols_all = [c for c in df.columns if c not in META_COLS and c not in NON_CATEGORY_COLS]
    cog_cols = [c for c in COG_ORDER if c in cog_cols_all] + [c for c in cog_cols_all if c not in COG_ORDER]
    df = df.set_index("sample_id")

    missing = [s for s in leaf_order if s not in df.index]
    if missing:
        print(f"[WARN] {len(missing)} tree leaves missing from COG data: {missing[:5]}...")
    heat_data = df.reindex(leaf_order)[cog_cols].to_numpy(dtype=float)

    im = ax_heat.imshow(heat_data, aspect="auto", cmap="YlOrRd",
                         extent=[0, len(cog_cols), n + 0.5, 0.5])
    ax_heat.set_xticks(np.arange(len(cog_cols)) + 0.5)
    ax_heat.set_xticklabels(cog_cols, fontsize=7, rotation=90)
    ax_heat.set_yticks([])
    ax_heat.set_title("COG profile (%)")
    fig.colorbar(im, ax=ax_heat, fraction=0.03, pad=0.02, label="Relative abundance (%)")

    fig.tight_layout()
    out_path = os.path.join(args.outdir, f"{args.genus}_tree_cog_heatmap.pdf")
    fig.savefig(out_path)
    print(f"Combined figure saved: {out_path}")


if __name__ == "__main__":
    main()
