#!/usr/bin/env python3
"""
diagnose_pangenome_genomes.py

각 그룹(genus/functional_group)의 Panaroo 결과에서, 시료별로
  - 원본 gff에 실제 어노테이션된 CDS 개수 (원본 데이터 품질)
  - Panaroo gene_presence_absence 매트릭스에서 그 시료가 '보유'로 잡힌 유전자 패밀리 수
를 비교합니다. 두 값의 비율이 비정상적으로 낮은 시료는 어노테이션이 부실하거나
(consolidate 단계에서 contig 누락 등) panaroo 처리 중 문제가 생겼을 가능성이 높습니다.
core_genes=0으로 나오는 그룹들의 원인 진단용입니다.

사용법:
  python3 diagnose_pangenome_genomes.py \
      --grouped-root /mnt/f/WGS_Consolidated/grouped \
      --pangenome-root /mnt/f/WGS_Consolidated/pangenome \
      --out /mnt/f/WGS_Consolidated/pangenome/genome_qc_report.csv \
      --ratio-threshold 0.5
"""
import argparse
import os

import pandas as pd


def count_cds(gff_path: str) -> int:
    if not os.path.isfile(gff_path):
        return -1
    n = 0
    with open(gff_path, "r", errors="ignore") as f:
        for line in f:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3 and fields[2] == "CDS":
                n += 1
    return n


def load_presence_counts(pangenome_group_dir: str) -> dict:
    """gene_presence_absence.Rtab (선호) 또는 .csv에서 시료별 보유 유전자 패밀리 수를 반환."""
    rtab_path = os.path.join(pangenome_group_dir, "gene_presence_absence.Rtab")
    if os.path.isfile(rtab_path):
        df = pd.read_csv(rtab_path, sep="\t", index_col=0)
        return df.sum(axis=0).to_dict()

    csv_path = os.path.join(pangenome_group_dir, "gene_presence_absence.csv")
    if os.path.isfile(csv_path):
        df = pd.read_csv(csv_path, low_memory=False)
        # panaroo csv: 앞쪽 메타컬럼들 이후 각 시료 컬럼에 유전자ID(있음) 또는 빈값(없음)
        meta_cols = {"Gene", "Non-unique Gene name", "Annotation"}
        sample_cols = [c for c in df.columns if c not in meta_cols]
        return {c: df[c].notna().sum() for c in sample_cols}

    return {}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grouped-root", required=True)
    ap.add_argument("--pangenome-root", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--ratio-threshold", type=float, default=0.5,
                     help="panaroo_gene_count / original_cds_count 가 이 값 미만이면 의심 시료로 표시")
    args = ap.parse_args()

    rows = []
    for root, dirs, files in os.walk(args.pangenome_root):
        if "summary_statistics.txt" not in files:
            continue
        rel = os.path.relpath(root, args.pangenome_root)
        parts = rel.split(os.sep)
        if len(parts) < 2:
            continue
        category, group = parts[0], parts[1]

        presence_counts = load_presence_counts(root)
        if not presence_counts:
            print(f"[경고] {category}/{group}: gene_presence_absence 파일을 못 찾음")
            continue

        gff_dir = os.path.join(args.grouped_root, category, group, "genomes_gff")

        for sample_id, panaroo_n in presence_counts.items():
            gff_path = os.path.join(gff_dir, f"{sample_id}.gff")
            orig_cds = count_cds(gff_path)
            ratio = (panaroo_n / orig_cds) if orig_cds and orig_cds > 0 else None
            rows.append({
                "category": category, "group": group, "sample_id": sample_id,
                "original_cds_count": orig_cds,
                "panaroo_gene_count": panaroo_n,
                "ratio": round(ratio, 3) if ratio is not None else None,
            })

    df = pd.DataFrame(rows)
    df.to_csv(args.out, index=False, encoding="utf-8-sig")
    print(f"저장 완료: {args.out} (총 {len(df)}건)")

    suspect = df[(df["ratio"].notna()) & (df["ratio"] < args.ratio_threshold)]
    suspect = suspect.sort_values("ratio")
    print(f"\n[의심 시료: ratio < {args.ratio_threshold}] {len(suspect)}건")
    if not suspect.empty:
        print(suspect.to_string(index=False))
    else:
        print("  없음 - 시료별 데이터 품질은 정상으로 보입니다. (core=0 원인은 다른 곳에 있을 수 있음)")

    print(f"\n[그룹별 original_cds_count / panaroo_gene_count 중앙값 비교]")
    med = df.groupby(["category", "group"])[["original_cds_count", "panaroo_gene_count"]].median()
    print(med.to_string())


if __name__ == "__main__":
    main()
