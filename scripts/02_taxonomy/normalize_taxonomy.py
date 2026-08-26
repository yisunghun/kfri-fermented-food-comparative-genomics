#!/usr/bin/env python3
"""
normalize_taxonomy.py

master_table.tsv의 species_guess에서 구(舊) Lactobacillus 속명을
2020년 재분류(Zheng et al., 2020, "A taxonomic note on the genus Lactobacillus")
기준 신규 속명으로 정규화합니다.

배경: BLAST 참조 DB(NCBI nt 등)에 신/구 명명법 레코드가 혼재되어 있어서,
완전히 동일한 종(strain-level로 거의 동일한 genome, ANI 99%+)인데도
'Lactobacillus plantarum' vs 'Lactiplantibacillus plantarum'처럼
서로 다른 속명으로 판정되는 경우가 많습니다. ani_species_mismatch.tsv에서
나온 불일치 대다수가 이 케이스에 해당합니다.

주의: 딕셔너리는 이번 데이터셋에서 관측된 종 위주로 구성했습니다.
     새로운 종이 추가되면 SPECIES_TO_NEW_GENUS에 계속 보강하세요.

사용법:
  python3 normalize_taxonomy.py --master master_table.tsv --out master_table_normalized.tsv
"""
import argparse

import pandas as pd

# 종소명(epithet) -> 신규 속명
SPECIES_TO_NEW_GENUS = {
    "plantarum": "Lactiplantibacillus",
    "paraplantarum": "Lactiplantibacillus",
    "pentosus": "Lactiplantibacillus",
    "brevis": "Levilactobacillus",
    "curvatus": "Latilactobacillus",
    "sakei": "Latilactobacillus",
    "casei": "Lacticaseibacillus",
    "paracasei": "Lacticaseibacillus",
    "rhamnosus": "Lacticaseibacillus",
    "fermentum": "Limosilactobacillus",
    "reuteri": "Limosilactobacillus",
    "buchneri": "Lentilactobacillus",
    "coryniformis": "Loigolactobacillus",
    "acidophilus": "Lactobacillus",  # delbrueckii group -> 속명 유지 (구분용으로 명시)
}


def normalize_species(species_guess: str) -> str:
    if not isinstance(species_guess, str) or " " not in species_guess:
        return species_guess
    genus, epithet = species_guess.split(" ", 1)
    epithet_key = epithet.split()[0].rstrip(".")
    if genus == "Lactobacillus" and epithet_key in SPECIES_TO_NEW_GENUS:
        new_genus = SPECIES_TO_NEW_GENUS[epithet_key]
        return f"{new_genus} {epithet}"
    return species_guess


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--master", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    df = pd.read_csv(args.master, sep="\t")
    df["species_normalized"] = df["species_guess"].apply(normalize_species)
    changed_mask = df["species_normalized"] != df["species_guess"]
    df.to_csv(args.out, sep="\t", index=False, encoding="utf-8-sig")

    print(f"정규화 완료: {changed_mask.sum()}건 속명 변경 -> {args.out}")
    if changed_mask.any():
        print("\n[변경 예시 (원본 종명 기준 중복제거)]")
        diff = df.loc[changed_mask, ["species_guess", "species_normalized"]].drop_duplicates()
        print(diff.to_string(index=False))

    print("\n[정규화 안 된 'Lactobacillus ...' 잔여 항목] (SPECIES_TO_NEW_GENUS 보강 필요할 수 있음)")
    remaining = df[df["species_normalized"].str.startswith("Lactobacillus ", na=False)]
    if remaining.empty:
        print("  없음")
    else:
        print(remaining["species_normalized"].value_counts().to_string())
