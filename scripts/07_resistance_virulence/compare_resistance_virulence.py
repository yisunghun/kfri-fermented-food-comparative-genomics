#!/usr/bin/env python3
"""
compare_resistance_virulence.py

abricate --summary output (card_summary.tsv or vfdb_summary.tsv) compared
between two functional groups (default: LAB vs Bacillus_group).

For each gene: presence/absence Fisher's exact test + BH-FDR correction.
Also compares total gene burden (NUM_FOUND) between groups (Mann-Whitney U).

Outputs:
  1) <db_label>_gene_stats.csv        - per-gene Fisher's exact test results
  2) <db_label>_top_genes_barplot.pdf - top significant genes, presence % by group
  3) <db_label>_burden_boxplot.pdf    - total gene count per genome, by group

Usage:
  python3 compare_resistance_virulence.py \
      --summary-tsv card_summary.tsv \
      --master master_table_qc.tsv \
      --db-label CARD \
      --outdir ./resistance_comparison \
      --group-a LAB --group-b Bacillus_group
"""
import argparse
import os
import re

import numpy as np
import pandas as pd
from scipy import stats

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


def functional_group(genus: str) -> str:
    if genus == "unresolved":
        return "Unresolved"
    if genus in LAB_GENERA:
        return "LAB"
    if genus in BACILLUS_GROUP_GENERA:
        return "Bacillus_group"
    return "Other_Environmental"


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    pvals = np.where(np.isnan(pvals), 1.0, pvals)
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order] * n / (np.arange(n) + 1)
    ranked_min = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked_min, 0, 1)
    return out


def sample_id_from_file(path: str) -> str:
    base = os.path.basename(str(path))
    return re.sub(r"\.tab$", "", base)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary-tsv", required=True)
    ap.add_argument("--master", required=True)
    ap.add_argument("--db-label", required=True, help="e.g. CARD or VFDB (used in filenames/titles)")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--group-a", default="LAB")
    ap.add_argument("--group-b", default="Bacillus_group")
    ap.add_argument("--top-n", type=int, default=20)
    args = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(args.outdir, exist_ok=True)

    summary = pd.read_csv(args.summary_tsv, sep="\t")
    summary["sample_id"] = summary["#FILE"].apply(sample_id_from_file)
    gene_cols = [c for c in summary.columns if c not in {"#FILE", "NUM_FOUND", "sample_id"}]
    presence = summary[gene_cols].map(lambda v: 0 if str(v).strip() in {".", "nan", ""} else 1)
    presence["sample_id"] = summary["sample_id"]
    presence["NUM_FOUND"] = summary["NUM_FOUND"]

    master = pd.read_csv(args.master, sep="\t")
    master["functional_group"] = master["genus_final"].apply(functional_group)

    merged = presence.merge(master[["sample_id", "genus_final", "functional_group"]], on="sample_id", how="inner")
    sub = merged[merged["functional_group"].isin([args.group_a, args.group_b])].copy()
    n_a = (sub["functional_group"] == args.group_a).sum()
    n_b = (sub["functional_group"] == args.group_b).sum()
    print(f"{args.db_label}: {args.group_a} n={n_a}, {args.group_b} n={n_b}")

    # ---------- 1) per-gene Fisher's exact test ----------
    rows = []
    a_mask = sub["functional_group"] == args.group_a
    b_mask = sub["functional_group"] == args.group_b
    for gene in gene_cols:
        a_pos = int(sub.loc[a_mask, gene].sum())
        b_pos = int(sub.loc[b_mask, gene].sum())
        if a_pos == 0 and b_pos == 0:
            continue  # skip genes absent in both groups entirely
        table = [[a_pos, n_a - a_pos], [b_pos, n_b - b_pos]]
        odds_ratio, p = stats.fisher_exact(table)
        rows.append({
            "gene": gene,
            f"n_positive_{args.group_a}": a_pos, f"pct_{args.group_a}": round(a_pos / n_a * 100, 1),
            f"n_positive_{args.group_b}": b_pos, f"pct_{args.group_b}": round(b_pos / n_b * 100, 1),
            "odds_ratio": odds_ratio, "p_value": p,
        })
    stats_df = pd.DataFrame(rows)
    if not stats_df.empty:
        stats_df["p_adj_BH"] = bh_fdr(stats_df["p_value"].to_numpy())
        stats_df["significant_FDR<0.05"] = stats_df["p_adj_BH"] < 0.05
        stats_df = stats_df.sort_values("p_adj_BH")
    stats_csv = os.path.join(args.outdir, f"{args.db_label}_gene_stats.csv")
    stats_df.to_csv(stats_csv, index=False, encoding="utf-8-sig")
    print(f"Gene-level stats saved: {stats_csv} ({len(stats_df)} genes tested)")
    sig = stats_df[stats_df.get("significant_FDR<0.05", False) == True] if not stats_df.empty else stats_df
    print(f"Significant genes (FDR<0.05): {len(sig)}")
    if not sig.empty:
        print(sig.head(20).to_string(index=False))

    # ---------- 2) top genes bar plot ----------
    if not stats_df.empty:
        top = stats_df.head(args.top_n)
        fig, ax = plt.subplots(figsize=(10, max(4, len(top) * 0.35)))
        y = np.arange(len(top))
        ax.barh(y - 0.2, top[f"pct_{args.group_a}"], height=0.4, label=args.group_a)
        ax.barh(y + 0.2, top[f"pct_{args.group_b}"], height=0.4, label=args.group_b)
        ax.set_yticks(y)
        ax.set_yticklabels(top["gene"], fontsize=8)
        ax.invert_yaxis()
        ax.set_xlabel("% of genomes with gene detected")
        ax.set_title(f"{args.db_label}: Top {len(top)} genes by group difference "
                     f"({args.group_a} vs {args.group_b})")
        ax.legend()
        fig.tight_layout()
        bar_path = os.path.join(args.outdir, f"{args.db_label}_top_genes_barplot.pdf")
        fig.savefig(bar_path)
        print(f"Bar plot saved: {bar_path}")

    # ---------- 3) total gene burden comparison ----------
    fig2, ax2 = plt.subplots(figsize=(6, 5))
    data_a = sub.loc[a_mask, "NUM_FOUND"]
    data_b = sub.loc[b_mask, "NUM_FOUND"]
    ax2.boxplot([data_a, data_b], tick_labels=[args.group_a, args.group_b])
    u, p_burden = stats.mannwhitneyu(data_a, data_b, alternative="two-sided")
    ax2.set_ylabel(f"Number of {args.db_label} genes detected per genome")
    ax2.set_title(f"{args.db_label} gene burden: {args.group_a} vs {args.group_b}\n"
                  f"Mann-Whitney U p={p_burden:.2e}")
    fig2.tight_layout()
    box_path = os.path.join(args.outdir, f"{args.db_label}_burden_boxplot.pdf")
    fig2.savefig(box_path)
    print(f"Burden box plot saved: {box_path}")
    print(f"Mean {args.db_label} gene count - {args.group_a}: {data_a.mean():.2f}, "
          f"{args.group_b}: {data_b.mean():.2f} (Mann-Whitney p={p_burden:.2e})")


if __name__ == "__main__":
    main()
