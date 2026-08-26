#!/usr/bin/env python3
"""
summarize_pangenomes.py

pangenome/<category>/<group>/summary_statistics.txt 들을 모두 모아
그룹별 core/soft-core/shell/cloud 유전자 개수와 pangenome 크기를 한 표로 정리합니다.

Panaroo의 summary_statistics.txt 포맷 예시:
    Core genes  (99% <= strains <= 100%)   2134
    Soft core genes (95% <= strains < 99%) 45
    Shell genes (15% <= strains < 95%)     612
    Cloud genes (0% <= strains < 15%)      1890
    Total genes (0% <= strains <= 100%)    4681

사용법:
  python3 summarize_pangenomes.py --pangenome-root /mnt/f/WGS_Consolidated/pangenome \
      --out /mnt/f/WGS_Consolidated/pangenome/pangenome_summary.csv
"""
import argparse
import os
import re

import pandas as pd

CATEGORY_MAP = {
    "Core genes": "core_genes",
    "Soft core genes": "soft_core_genes",
    "Shell genes": "shell_genes",
    "Cloud genes": "cloud_genes",
    "Total genes": "total_genes",
}


def parse_summary(path: str) -> dict:
    result = {}
    with open(path, "r") as f:
        for line in f:
            for label, key in CATEGORY_MAP.items():
                if line.strip().startswith(label):
                    m = re.search(r"(\d+)\s*$", line.strip())
                    if m:
                        result[key] = int(m.group(1))
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pangenome-root", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    rows = []
    for root, dirs, files in os.walk(args.pangenome_root):
        if "summary_statistics.txt" in files:
            stats_path = os.path.join(root, "summary_statistics.txt")
            rel = os.path.relpath(root, args.pangenome_root)
            parts = rel.split(os.sep)
            category = parts[0] if len(parts) > 0 else "?"
            group = parts[1] if len(parts) > 1 else "?"

            gff_dir_guess = None  # n(genome) 개수는 group_manifest.tsv에서 가져옴
            manifest_path = None
            for cand in [os.path.join(args.pangenome_root, "..", "grouped", category, group, "group_manifest.tsv")]:
                if os.path.isfile(cand):
                    manifest_path = cand
                    break

            stats = parse_summary(stats_path)
            n_genomes = None
            if manifest_path:
                try:
                    n_genomes = len(pd.read_csv(manifest_path, sep="\t"))
                except Exception:
                    pass

            row = {"category": category, "group": group, "n_genomes": n_genomes}
            row.update(stats)
            rows.append(row)

    df = pd.DataFrame(rows)
    for col in ["core_genes", "soft_core_genes", "shell_genes", "cloud_genes", "total_genes"]:
        if col not in df.columns:
            df[col] = pd.NA

    # 보조 지표: core genome 비율(core / total), accessory 비율
    df["core_pct_of_total"] = (df["core_genes"] / df["total_genes"] * 100).round(1)
    df["accessory_genes"] = df["shell_genes"].fillna(0) + df["cloud_genes"].fillna(0)

    df = df.sort_values(["category", "n_genomes"], ascending=[True, False])
    df.to_csv(args.out, index=False, encoding="utf-8-sig")

    print(f"저장 완료: {args.out}\n")
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
