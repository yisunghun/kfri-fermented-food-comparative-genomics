#!/usr/bin/env python3
"""
extract_eggnog_summary.py

각 시료의 annotation_EggNOG.xlsx 안 'Eggnog_Count' 시트
(컬럼: Eggnog, Description, Count, Ratio (%))를 모아
시료 x COG카테고리 매트릭스(Count 기준, Ratio% 기준 둘 다)로 만듭니다.
genus_final / functional_group 라벨도 같이 붙여서 그룹별 비교에 바로 쓸 수 있게 합니다.

사용법:
  python3 extract_eggnog_summary.py \
      --master master_table_qc.tsv \
      --out-count eggnog_cog_count_wide.tsv \
      --out-ratio eggnog_cog_ratio_wide.tsv
"""
import argparse
import os

import pandas as pd

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


def normalize_sheet_name(name: str) -> str:
    return name.lower().replace(" ", "").replace("_", "")


def find_matching_sheet(xls: pd.ExcelFile, keyword: str) -> str:
    """정규화(공백/언더스코어 제거, 소문자화) 후 keyword를 포함하는 시트명을 찾음."""
    target = normalize_sheet_name(keyword)
    for name in xls.sheet_names:
        if normalize_sheet_name(name) == target:
            return name
    # 완전일치 실패 시 부분포함으로 재시도
    for name in xls.sheet_names:
        if target in normalize_sheet_name(name):
            return name
    return None


def eggnog_path_for_sample(row) -> str:
    consensus_src = row.get("consensus_fasta_src")
    short_id = row.get("short_id")
    if not isinstance(consensus_src, str) or not isinstance(short_id, str):
        return ""
    assembly_dir = os.path.dirname(consensus_src)
    sample_top_dir = os.path.dirname(assembly_dir)
    return os.path.join(sample_top_dir, f"{short_id}_FunctionalAnnotation",
                         "FunctionalAnnotation", "annotation_EggNOG.xlsx")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--master", required=True)
    ap.add_argument("--out-count", required=True)
    ap.add_argument("--out-ratio", required=True)
    ap.add_argument("--sheet-name", default="Eggnog_Count")
    args = ap.parse_args()

    master = pd.read_csv(args.master, sep="\t")
    master["eggnog_path"] = master.apply(eggnog_path_for_sample, axis=1)
    master["functional_group"] = master["genus_final"].apply(functional_group)

    count_rows = []
    ratio_rows = []
    n_missing = 0
    n_error = 0

    for _, row in master.iterrows():
        sid = row["sample_id"]
        p = row["eggnog_path"]
        if not p or not os.path.isfile(p):
            n_missing += 1
            continue
        try:
            xls = pd.ExcelFile(p)
            sheet = find_matching_sheet(xls, args.sheet_name)
            if sheet is None:
                print(f"[경고] {sid}: '{args.sheet_name}' 계열 시트를 못 찾음 "
                      f"(실제 시트: {xls.sheet_names}), 스킵")
                n_error += 1
                continue
            df = pd.read_excel(xls, sheet_name=sheet)
            # 같은 카테고리 코드가 여러 행으로 중복될 수 있어 groupby sum으로 합산
            count_series = df.groupby("Eggnog")["Count"].sum()
            ratio_series = df.groupby("Eggnog")["Ratio (%)"].sum()
            count_rows.append(count_series.rename(sid))
            ratio_rows.append(ratio_series.rename(sid))
        except Exception as e:
            n_error += 1
            print(f"[오류] {sid}: {e}")

    if n_missing:
        print(f"[경고] eggNOG xlsx 파일 누락: {n_missing}건")
    if n_error:
        print(f"[경고] 파싱 오류: {n_error}건")

    meta_cols = ["sample_id", "genus_final", "species_final", "functional_group"]

    for rows, out_path, label in [(count_rows, args.out_count, "Count"),
                                   (ratio_rows, args.out_ratio, "Ratio(%)")]:
        wide = pd.DataFrame(rows).fillna(0)
        wide.index.name = "sample_id"
        wide = wide.reset_index()
        merged = master[meta_cols].merge(wide, on="sample_id", how="inner")
        merged.to_csv(out_path, sep="\t", index=False, encoding="utf-8-sig")
        print(f"저장 완료 ({label} 기준): {out_path} ({len(merged)}개 시료)")


if __name__ == "__main__":
    main()
