#!/usr/bin/env python3
"""
filter_to_dereplicated.py

대표 시료 목록(representative_samples.txt)만 남기도록 기존 데이터 파일들을
필터링한다. 이렇게 만든 _dereplicated 파일들을 기존 분석 스크립트
(compare_functional_groups.py, compare_resistance_virulence.py 등)에 그대로
다시 입력하면, 통계 로직을 새로 짜지 않고도 민감도 분석(클론 제거 후 재검정)이
가능하다.

사용법:
    python3 filter_to_dereplicated.py \
        --representatives dereplication/representative_samples.txt \
        --master master_table_qc.tsv \
        --extra-tsv eggnog_cog_ratio_wide.tsv \
        --extra-tsv card_summary.tsv \
        --extra-tsv vfdb_summary.tsv \
        --outdir dereplication/filtered
"""
import argparse
import os
import pandas as pd


import re


def derive_sample_id(raw_value):
    """#FILE 같은 컬럼 값(예: HN00157641_AMT60212.tab, 혹은 경로 포함)에서
    확장자/경로를 제거해 순수 sample_id로 정규화한다."""
    base = os.path.basename(str(raw_value))
    base = re.sub(r"\.(tab|fna|fa|fasta|txt|csv|tsv)$", "", base)
    return base


def filter_tsv(path, keep_samples, outdir):
    df = pd.read_csv(path, sep="\t")
    # sample_id 컬럼이 있으면 그걸 그대로, #FILE류 컬럼이면 확장자/경로를 제거해 비교
    if "sample_id" in df.columns:
        id_col = "sample_id"
        compare_series = df[id_col].astype(str)
    else:
        id_col = df.columns[0]
        compare_series = df[id_col].apply(derive_sample_id)

    before = len(df)
    filtered = df[compare_series.isin(keep_samples)].copy()
    after = len(filtered)

    base = os.path.splitext(os.path.basename(path))[0]
    out_path = os.path.join(outdir, f"{base}_dereplicated.tsv")
    filtered.to_csv(out_path, sep="\t", index=False)
    print(f"{path}: {before} -> {after}행 (식별자 컬럼: {id_col}) -> {out_path}")
    if after == 0:
        print(f"  [경고] 0행으로 필터링됨! 컬럼 '{id_col}' 값 예시: {df[id_col].head(3).tolist()}")
    return out_path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--representatives", required=True, help="representative_samples.txt")
    ap.add_argument("--master", required=True)
    ap.add_argument("--extra-tsv", action="append", default=[],
                     help="추가로 필터링할 TSV (여러 번 지정 가능): eggnog_cog_ratio_wide.tsv, card_summary.tsv, vfdb_summary.tsv 등")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    with open(args.representatives) as f:
        keep_samples = set(line.strip() for line in f if line.strip())
    print(f"대표 시료 수: {len(keep_samples)}개\n")

    filter_tsv(args.master, keep_samples, args.outdir)
    for tsv in args.extra_tsv:
        filter_tsv(tsv, keep_samples, args.outdir)

    print("\n완료. 위 _dereplicated.tsv 파일들을 기존 분석 스크립트에 다시 입력하면 됩니다. 예:")
    print("  python3 compare_functional_groups.py --ratio-tsv <outdir>/eggnog_cog_ratio_wide_dereplicated.tsv --outdir dereplication/cog_result --group-a LAB --group-b Bacillus_group")
    print("  python3 compare_resistance_virulence.py --summary-tsv <outdir>/card_summary_dereplicated.tsv --master <outdir>/master_table_qc_dereplicated.tsv --db-label CARD --outdir dereplication/card_result --group-a LAB --group-b Bacillus_group")


if __name__ == "__main__":
    main()
