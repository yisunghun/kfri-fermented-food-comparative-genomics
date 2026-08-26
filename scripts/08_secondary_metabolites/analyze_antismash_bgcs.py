#!/usr/bin/env python3
"""
analyze_antismash_bgcs.py

Parses antiSMASH JSON output (antismash_out/<sample_id>/<sample_id>.json) to
build a sample x BGC-product-type count matrix, then compares two functional
groups (default LAB vs Bacillus_group): total BGC burden (Mann-Whitney U) and
per-product-type presence (Fisher's exact test + BH-FDR).

JSON schema (antiSMASH 7+):
  records: [ { areas: [ { start, end, products: [str, ...], ... }, ... ] }, ... ]

Outputs:
  1) antismash_bgc_summary.tsv     - sample x product-type count matrix
  2) antismash_burden_boxplot.pdf  - total BGC count per genome, by group
  3) antismash_product_barplot.pdf - mean BGC count per product type, by group
  4) antismash_product_stats.csv   - per product-type Fisher's exact test (presence)

Usage:
  python3 analyze_antismash_bgcs.py \
      --antismash-root /mnt/f/WGS_Consolidated/antismash_out \
      --master /mnt/f/WGS_Consolidated/master_table_qc.tsv \
      --outdir /mnt/f/WGS_Consolidated/antismash_analysis \
      --group-a LAB --group-b Bacillus_group
"""
import argparse
import json
import os

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


def parse_sample_bgcs(json_path: str) -> tuple[int, list[str]]:
    """Returns (total_bgc_count, flattened list of product-type strings)."""
    with open(json_path) as f:
        d = json.load(f)
    total = 0
    products = []
    for record in d.get("records", []):
        for area in record.get("areas", []):
            total += 1
            for p in area.get("products", []):
                products.append(p)
    return total, products


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--antismash-root", required=True)
    ap.add_argument("--master", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--group-a", default="LAB")
    ap.add_argument("--group-b", default="Bacillus_group")
    args = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(args.outdir, exist_ok=True)

    master = pd.read_csv(args.master, sep="\t")
    master["functional_group"] = master["genus_final"].apply(functional_group)

    rows = []
    n_missing = 0
    for sample_id in master["sample_id"]:
        json_path = os.path.join(args.antismash_root, sample_id, f"{sample_id}.json")
        if not os.path.isfile(json_path):
            n_missing += 1
            continue
        total, products = parse_sample_bgcs(json_path)
        row = {"sample_id": sample_id, "total_bgc_count": total}
        for p in products:
            row[p] = row.get(p, 0) + 1
        rows.append(row)

    if n_missing:
        print(f"[WARN] antiSMASH json not found for {n_missing} samples (skipped)")

    wide = pd.DataFrame(rows).fillna(0)
    product_cols = [c for c in wide.columns if c not in {"sample_id", "total_bgc_count"}]
    wide[product_cols] = wide[product_cols].astype(int)

    merged = wide.merge(master[["sample_id", "genus_final", "functional_group"]], on="sample_id", how="left")
    summary_path = os.path.join(args.outdir, "antismash_bgc_summary.tsv")
    merged.to_csv(summary_path, sep="\t", index=False, encoding="utf-8-sig")
    print(f"Sample-level BGC summary saved: {summary_path} ({len(merged)} samples, "
          f"{len(product_cols)} product types)")

    sub = merged[merged["functional_group"].isin([args.group_a, args.group_b])].copy()
    a_mask = sub["functional_group"] == args.group_a
    b_mask = sub["functional_group"] == args.group_b
    n_a, n_b = a_mask.sum(), b_mask.sum()
    print(f"{args.group_a} n={n_a}, {args.group_b} n={n_b}")

    # ---------- 1) total BGC burden boxplot ----------
    data_a = sub.loc[a_mask, "total_bgc_count"]
    data_b = sub.loc[b_mask, "total_bgc_count"]
    u, p_burden = stats.mannwhitneyu(data_a, data_b, alternative="two-sided")

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.boxplot([data_a, data_b], tick_labels=[args.group_a, args.group_b])
    ax.set_ylabel("Total BGC count per genome")
    ax.set_title(f"antiSMASH BGC burden: {args.group_a} vs {args.group_b}\n"
                 f"Mann-Whitney U p={p_burden:.2e}")
    fig.tight_layout()
    burden_path = os.path.join(args.outdir, "antismash_burden_boxplot.pdf")
    fig.savefig(burden_path)
    print(f"Burden box plot saved: {burden_path}")
    print(f"Mean BGC count - {args.group_a}: {data_a.mean():.2f}, {args.group_b}: {data_b.mean():.2f} "
          f"(Mann-Whitney p={p_burden:.2e})")

    # ---------- 2) per product-type Fisher's exact test (presence) ----------
    stat_rows = []
    for p in product_cols:
        a_pos = int((sub.loc[a_mask, p] > 0).sum())
        b_pos = int((sub.loc[b_mask, p] > 0).sum())
        if a_pos == 0 and b_pos == 0:
            continue
        table = [[a_pos, n_a - a_pos], [b_pos, n_b - b_pos]]
        odds_ratio, p_val = stats.fisher_exact(table)
        stat_rows.append({
            "bgc_product_type": p,
            f"n_positive_{args.group_a}": a_pos, f"pct_{args.group_a}": round(a_pos / n_a * 100, 1),
            f"n_positive_{args.group_b}": b_pos, f"pct_{args.group_b}": round(b_pos / n_b * 100, 1),
            "odds_ratio": odds_ratio, "p_value": p_val,
        })
    stats_df = pd.DataFrame(stat_rows)
    if not stats_df.empty:
        stats_df["p_adj_BH"] = bh_fdr(stats_df["p_value"].to_numpy())
        stats_df["significant_FDR<0.05"] = stats_df["p_adj_BH"] < 0.05
        stats_df = stats_df.sort_values("p_adj_BH")
    stats_csv = os.path.join(args.outdir, "antismash_product_stats.csv")
    stats_df.to_csv(stats_csv, index=False, encoding="utf-8-sig")
    print(f"Product-type stats saved: {stats_csv}")
    if not stats_df.empty:
        sig = stats_df[stats_df["significant_FDR<0.05"]]
        print(f"Significant BGC product types (FDR<0.05): {len(sig)}")
        print(sig.to_string(index=False) if not sig.empty else "  none")

    # ---------- 3) mean BGC count per product type bar plot ----------
    if product_cols:
        means = sub.groupby("functional_group")[product_cols].mean().T
        means = means.reindex(means.sum(axis=1).sort_values(ascending=False).index)
        fig2, ax2 = plt.subplots(figsize=(max(8, len(means) * 0.5), 6))
        x = np.arange(len(means))
        width = 0.35
        ax2.bar(x - width / 2, means[args.group_a], width, label=args.group_a)
        ax2.bar(x + width / 2, means[args.group_b], width, label=args.group_b)
        ax2.set_xticks(x)
        ax2.set_xticklabels(means.index, rotation=45, ha="right", fontsize=8)
        ax2.set_ylabel("Mean BGC count per genome")
        ax2.set_title(f"antiSMASH BGC product types: {args.group_a} vs {args.group_b}")
        ax2.legend()
        fig2.tight_layout()
        bar_path = os.path.join(args.outdir, "antismash_product_barplot.pdf")
        fig2.savefig(bar_path)
        print(f"Product-type bar plot saved: {bar_path}")


if __name__ == "__main__":
    main()
