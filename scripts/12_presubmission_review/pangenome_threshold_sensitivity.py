#!/usr/bin/env python3
"""
pangenome_threshold_sensitivity.py

3.4절의 species-level pangenome(B. velezensis, L. plantarum)이 Panaroo의 기본
core-gene 임계값(99%)에 지나치게 민감한 결과는 아닌지 확인하기 위해,
gene_presence_absence.csv를 이용해 여러 임계값(100/99/95/90/85/80%)에서
core genome 비율을 재계산한다.

사용법:
    python3 pangenome_threshold_sensitivity.py \
        --panaroo-csv species_pangenome/Bacillus_velezensis/panaroo_out/gene_presence_absence.csv \
        --label "Bacillus velezensis" \
        --outdir threshold_sensitivity_result
    (L. plantarum도 동일하게 별도 실행)
"""
import argparse
import os
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panaroo-csv", required=True,
                     help="Panaroo의 gene_presence_absence.csv 경로")
    ap.add_argument("--label", required=True, help="종 이름 (결과 표시용)")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    df = pd.read_csv(args.panaroo_csv, low_memory=False)
    # Panaroo 표준 메타데이터 컬럼 (genome이 아닌 것들) 제외하고 genome 컬럼만 추출
    meta_cols = {
        "Gene", "Non-unique Gene name", "Annotation", "No. isolates",
        "No. sequences", "Avg sequences per isolate", "Genome Fragment",
        "Order within Fragment", "Accessory Fragment", "Accessory Order with Fragment",
        "QC", "Min group size nuc", "Max group size nuc", "Avg group size nuc",
    }
    genome_cols = [c for c in df.columns if c not in meta_cols]
    n_genomes = len(genome_cols)
    print(f"=== {args.label}: genome 수 = {n_genomes}, gene family 수 = {len(df)} ===")

    # 각 gene family가 몇 개 genome에 존재하는지 (빈 문자열/NaN이 아니면 존재)
    presence = df[genome_cols].notna() & (df[genome_cols] != "")
    n_present = presence.sum(axis=1)
    pct_present = n_present / n_genomes * 100

    thresholds = [100, 99, 95, 90, 85, 80]
    rows = []
    for t in thresholds:
        n_core = (pct_present >= t).sum()
        pct_core = n_core / len(df) * 100
        rows.append({"species": args.label, "threshold_pct": t,
                      "n_core_genes": n_core, "pangenome_size": len(df),
                      "core_pct_of_pangenome": round(pct_core, 1)})
        print(f"  threshold >= {t}%: core genes = {n_core} ({pct_core:.1f}% of pangenome)")

    out_df = pd.DataFrame(rows)
    safe_label = args.label.replace(" ", "_").replace(".", "")
    out_path = os.path.join(args.outdir, f"{safe_label}_threshold_sensitivity.csv")
    out_df.to_csv(out_path, index=False)
    print(f"저장됨: {out_path}")


if __name__ == "__main__":
    main()
