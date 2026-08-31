#!/usr/bin/env python3
"""
diagnose_missing_5_genomes.py

이번 세션에서 만든 여러 스크립트(add_effect_sizes.py, normalize_by_genome_size.py,
check_nheABC_denominator.py 등)가 공통으로 쓴 "폴더 기반 fg_map"이
Bacillus-group을 79개가 아니라 74개로 인식하는 문제의 정확한 원인을 진단한다.

1) master_table_qc.tsv 기준 진짜 Bacillus-group 79개 목록 확보
2) grouped/by_functional_group/Bacillus_group/genomes_gff/ 폴더에 있는 실제 파일 목록 확보
3) 두 목록을 대조해서 "빠진 5개"를 정확히 특정
4) 그 5개가 abricate_out/vfdb_summary.tsv, card_summary.tsv에 실제로 존재하는지 확인
   -> 존재하면: 폴더 매핑 스크립트의 버그 (데이터 자체는 있음, 재분석 필요)
   -> 존재 안 하면: 진짜로 이 5개는 VFDB/CARD 스크리닝이 안 된 것 (Methods에 사유 명시 필요)

사용법:
    python3 diagnose_missing_5_genomes.py \
        --master master_table_qc.tsv \
        --functional-group-root grouped/by_functional_group \
        --vfdb-summary abricate_out/vfdb_summary.tsv \
        --card-summary abricate_out/card_summary.tsv
"""
import argparse
import glob
import os
import re
import pandas as pd


def derive_sample_id(raw):
    base = os.path.basename(str(raw))
    return re.sub(r"\.(tab|fna|fa|fasta|txt|csv|tsv)$", "", base)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--master", required=True)
    ap.add_argument("--functional-group-root", required=True)
    ap.add_argument("--vfdb-summary", required=True)
    ap.add_argument("--card-summary", required=True)
    args = ap.parse_args()

    master = pd.read_csv(args.master, sep="\t")
    master.columns = [re.sub(r"^\ufeff", "", c) for c in master.columns]

    # 1) master_table_qc.tsv 기준 진짜 Bacillus-group 목록
    true_bacillus_group = set(master.loc[master["genus_final"].isin(
        ["Bacillus", "Paenibacillus", "Oceanobacillus"]), "sample_id"])
    print(f"=== master_table_qc.tsv 기준 Bacillus-group: {len(true_bacillus_group)}개 ===")

    # 2) 폴더 기반 fg_map이 실제로 찾는 목록
    fg_dir = os.path.join(args.functional_group_root, "Bacillus_group", "genomes_gff")
    if not os.path.isdir(fg_dir):
        candidates = glob.glob(os.path.join(args.functional_group_root, "*acillus*group*"))
        print(f"[경고] {fg_dir} 없음. 후보: {candidates}")
        return
    folder_found = set()
    for fpath in glob.glob(os.path.join(fg_dir, "*.gff")):
        sid = os.path.splitext(os.path.basename(fpath))[0]
        folder_found.add(sid)
    print(f"=== 폴더(genomes_gff) 기준 Bacillus-group: {len(folder_found)}개 ===")

    # 3) 빠진 5개 특정
    missing = true_bacillus_group - folder_found
    print(f"\n=== 폴더 매핑에서 빠진 시료: {len(missing)}개 ===")
    for s in sorted(missing):
        print(f"  {s}")

    if not missing:
        print("빠진 시료 없음 - 문제 없음")
        return

    # 4) 이 5개가 실제 VFDB/CARD summary에 존재하는지 확인
    vfdb = pd.read_csv(args.vfdb_summary, sep="\t")
    vfdb["sample_id"] = vfdb[vfdb.columns[0]].apply(derive_sample_id)
    card = pd.read_csv(args.card_summary, sep="\t")
    card["sample_id"] = card[card.columns[0]].apply(derive_sample_id)

    print(f"\n=== 빠진 시료들의 실제 VFDB/CARD 데이터 존재 여부 ===")
    for s in sorted(missing):
        in_vfdb = s in set(vfdb["sample_id"])
        in_card = s in set(card["sample_id"])
        print(f"  {s}: VFDB 데이터 존재={in_vfdb}, CARD 데이터 존재={in_card}")

    n_have_data = sum(1 for s in missing if s in set(vfdb["sample_id"]))
    print(f"\n결론: 빠진 {len(missing)}개 중 {n_have_data}개가 실제로는 VFDB 데이터를 갖고 있음.")
    if n_have_data == len(missing):
        print(">>> 이건 폴더 매핑 스크립트(fg_map)의 버그입니다. 데이터 자체는 정상이며,")
        print(">>> 이번 세션에서 만든 effect size/genome-size normalization/VFDB 방향 검증")
        print(">>> 스크립트들을 n=79 기준으로 재실행해야 합니다.")
    elif n_have_data == 0:
        print(">>> 이 5개는 실제로 VFDB/CARD 스크리닝 대상에서 빠진 것입니다.")
        print(">>> Methods에 이유를 명시해야 합니다 (예: 5개 QC-제외 genome과 일치하는지 확인 필요).")
    else:
        print(">>> 일부만 데이터가 있음 - 개별 확인 필요.")


if __name__ == "__main__":
    main()
