#!/usr/bin/env python3
"""
assess_genome_qc.py

220개 isolate의 consolidated genome fasta(genomes_fna/<sample_id>.fna)로부터
기본 assembly QC 지표를 계산한다:
  - 총 assembly 길이
  - contig 개수
  - N50, L50
  - 최대 contig 길이
  - GC%
  - 20kb 이상 contig 개수/비율

리뷰 지적사항(assembly quality 정보가 논문에 전혀 없음) 대응용 Supplementary
Table 생성 스크립트.

사용법:
    python3 assess_genome_qc.py \
        --genomes-dir genomes_fna \
        --master master_table_qc.tsv \
        --outdir genome_qc
"""
import argparse
import glob
import os
import pandas as pd
import numpy as np


def parse_fasta_lengths(path):
    """fasta 파일에서 각 레코드(contig)의 길이와 GC 개수를 계산 (Biopython 없이 순수 파싱)"""
    lengths = []
    gc_count = 0
    at_count = 0
    cur_len = 0
    cur_gc = 0
    cur_at = 0
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_len > 0:
                    lengths.append(cur_len)
                    gc_count += cur_gc
                    at_count += cur_at
                cur_len = 0
                cur_gc = 0
                cur_at = 0
            else:
                seq = line.upper()
                cur_len += len(seq)
                cur_gc += seq.count("G") + seq.count("C")
                cur_at += seq.count("A") + seq.count("T")
        if cur_len > 0:
            lengths.append(cur_len)
            gc_count += cur_gc
            at_count += cur_at
    return lengths, gc_count, at_count


def calc_n50_l50(lengths):
    lengths_sorted = sorted(lengths, reverse=True)
    total = sum(lengths_sorted)
    half = total / 2
    cum = 0
    for i, L in enumerate(lengths_sorted, 1):
        cum += L
        if cum >= half:
            return L, i
    return 0, 0


def load_functional_group_map(fg_root):
    """grouped/by_functional_group/<LAB|Bacillus_group|Other_Environmental|Unresolved>/
    폴더 구조 자체를 근거로 sample_id -> functional_group 매핑을 만든다.
    (genus 목록을 추측해서 다시 만들지 않고, 파이프라인이 이미 확정한 분류를 그대로 사용)"""
    mapping = {}
    if not fg_root or not os.path.isdir(fg_root):
        return mapping
    for fg_name in os.listdir(fg_root):
        fg_dir = os.path.join(fg_root, fg_name)
        if not os.path.isdir(fg_dir):
            continue
        # genomes_gff 하위의 파일명에서 sample_id 추출 (그룹 폴더 구조: <fg>/genomes_gff/<sample_id>.gff)
        gff_dir = os.path.join(fg_dir, "genomes_gff")
        search_dir = gff_dir if os.path.isdir(gff_dir) else fg_dir
        for fpath in glob.glob(os.path.join(search_dir, "*.gff")):
            sid = os.path.splitext(os.path.basename(fpath))[0]
            mapping[sid] = fg_name
        # manifest가 있으면 그것도 보조로 사용 (더 정확할 수 있음)
        manifest_path = os.path.join(fg_dir, "group_manifest.tsv")
        if os.path.isfile(manifest_path):
            try:
                mdf = pd.read_csv(manifest_path, sep="\t")
                if "sample_id" in mdf.columns:
                    for sid in mdf["sample_id"]:
                        mapping[sid] = fg_name
            except Exception:
                pass
    return mapping


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genomes-dir", required=True)
    ap.add_argument("--master", required=True)
    ap.add_argument("--functional-group-root", default=None,
                     help="grouped/by_functional_group 폴더 경로 (LAB/Bacillus_group/Other_Environmental/Unresolved 하위폴더가 있는 곳). "
                          "지정하면 이 폴더 구조를 근거로 functional_group을 정확히 매핑함.")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    master = pd.read_csv(args.master, sep="\t")
    master = master.set_index("sample_id")

    fg_map = load_functional_group_map(args.functional_group_root)
    if fg_map:
        print(f"functional_group 매핑 로드됨: {len(fg_map)}개 시료 ({args.functional_group_root} 기준)")
    else:
        print("[참고] --functional-group-root 미지정 또는 폴더 없음 -> functional_group 컬럼은 비어있게 됩니다.")

    fasta_files = sorted(glob.glob(os.path.join(args.genomes_dir, "*.fna")))
    print(f"발견된 genome fasta: {len(fasta_files)}개")

    rows = []
    for fpath in fasta_files:
        sample_id = os.path.splitext(os.path.basename(fpath))[0]
        lengths, gc, at = parse_fasta_lengths(fpath)
        if not lengths:
            print(f"[경고] {sample_id}: 서열을 읽지 못함")
            continue
        total_len = sum(lengths)
        n_contigs = len(lengths)
        n50, l50 = calc_n50_l50(lengths)
        largest = max(lengths)
        gc_pct = 100 * gc / (gc + at) if (gc + at) > 0 else np.nan
        n_ge20kb = sum(1 for L in lengths if L >= 20000)

        row = {
            "sample_id": sample_id,
            "total_length_bp": total_len,
            "n_contigs": n_contigs,
            "N50": n50,
            "L50": l50,
            "largest_contig_bp": largest,
            "GC_pct": round(gc_pct, 2),
            "n_contigs_ge20kb": n_ge20kb,
            "pct_length_in_contigs_ge20kb": round(
                100 * sum(L for L in lengths if L >= 20000) / total_len, 2
            ) if total_len > 0 else np.nan,
        }
        if sample_id in master.index:
            row["genus_final"] = master.loc[sample_id].get("genus_final", np.nan)
            row["species_final"] = master.loc[sample_id].get("species_final", np.nan)
        row["functional_group"] = fg_map.get(sample_id, np.nan)
        rows.append(row)

    df = pd.DataFrame(rows)
    out_path = os.path.join(args.outdir, "genome_qc_metrics.tsv")
    df.to_csv(out_path, sep="\t", index=False)

    print(f"\n=== 전체 220개 요약 통계 ===")
    print(df[["total_length_bp", "n_contigs", "N50", "GC_pct"]].describe().to_string())

    print(f"\n=== 이상치 후보 (참고용 기준) ===")
    # 지나치게 많은 contig 수 (>200), 지나치게 작은 N50 (<20kb), 비정상적 총 길이
    flagged = df[(df["n_contigs"] > 200) | (df["N50"] < 20000)]
    if len(flagged) > 0:
        print(f"contig>200 또는 N50<20kb 인 시료: {len(flagged)}개")
        print(flagged[["sample_id", "genus_final", "total_length_bp", "n_contigs", "N50"]].to_string(index=False))
    else:
        print("특별한 이상치 없음")

    print(f"\n저장됨: {out_path}")

    # 기능군별 요약
    if "functional_group" in df.columns:
        print(f"\n=== 기능군별 QC 요약 (평균) ===")
        summary = df.groupby("functional_group")[["total_length_bp", "n_contigs", "N50", "GC_pct"]].mean().round(1)
        print(summary.to_string())
        summary.to_csv(os.path.join(args.outdir, "genome_qc_by_functional_group.tsv"), sep="\t")


if __name__ == "__main__":
    main()
