#!/usr/bin/env python3
"""
replot_antismash_burden.py

analyze_antismash_bgcs.py가 원본 JSON 경로 문제로 실패해서, 이미 만들어진
antismash_bgc_summary.tsv를 직접 읽어 깔끔한 라벨(matplotlib 변수명 노출 없음,
p-value 제목에서 제거)로 antiSMASH burden boxplot을 다시 그린다.

사용법:
    python3 replot_antismash_burden.py \
        --bgc-summary antismash_analysis/antismash_bgc_summary.tsv \
        --master master_table_qc_with_fg.tsv \
        --outdir figures_fixed
"""
import argparse
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

DISP = {"LAB": "LAB", "Bacillus_group": "Bacillus-group"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bgc-summary", required=True)
    ap.add_argument("--master", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    bgc = pd.read_csv(args.bgc_summary, sep="\t")
    id_col = "sample_id" if "sample_id" in bgc.columns else bgc.columns[0]
    count_col = "total_bgc_count" if "total_bgc_count" in bgc.columns else bgc.columns[1]
    bgc = bgc[[id_col, count_col]].rename(columns={id_col: "sample_id", count_col: "bgc_count"})

    master = pd.read_csv(args.master, sep="\t")
    master.columns = [c.replace("\ufeff", "") for c in master.columns]
    merged = bgc.merge(master[["sample_id", "functional_group"]], on="sample_id", how="left")

    data_a = merged.loc[merged["functional_group"] == "LAB", "bgc_count"].dropna()
    data_b = merged.loc[merged["functional_group"] == "Bacillus_group", "bgc_count"].dropna()
    u, p = stats.mannwhitneyu(data_a, data_b, alternative="two-sided")
    print(f"n_LAB={len(data_a)}, n_Bacillus_group={len(data_b)}, p={p:.3g}")

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.boxplot([data_a, data_b], tick_labels=[DISP["LAB"], DISP["Bacillus_group"]])
    ax.set_ylabel("Total BGC count per genome")
    ax.set_title("antiSMASH BGC burden by functional group")
    fig.tight_layout()
    out_path = f"{args.outdir}/antismash_burden_boxplot_clean.pdf"
    fig.savefig(out_path)
    print(f"저장 완료: {out_path}")


if __name__ == "__main__":
    main()
