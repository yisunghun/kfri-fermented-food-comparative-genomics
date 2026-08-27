#!/usr/bin/env python3
"""
identify_species_pangenome_candidates.py

master_table_qc.tsv의 species_final 컬럼 기준으로 종별 시료 수를 집계해서,
종 수준(species-level) pangenome 분석이 의미 있을 만큼 시료가 충분한 종
(기본 임계값 n>=5)을 추천한다.

사용법:
    python3 identify_species_pangenome_candidates.py \
        --master master_table_qc.tsv \
        --min-n 5
"""
import argparse
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--master", required=True)
    ap.add_argument("--min-n", type=int, default=5)
    args = ap.parse_args()

    master = pd.read_csv(args.master, sep="\t")

    # sp.로 끝나는(속만 판정된) 것은 종 수준 분석 대상에서 제외
    resolved = master[~master["species_final"].astype(str).str.endswith("sp.")]

    counts = resolved.groupby(["genus_final", "species_final"]).size().reset_index(name="n")
    counts = counts.sort_values("n", ascending=False)

    print("=== 종(species_final)별 시료 수 (많은 순) ===")
    print(counts.to_string(index=False))

    candidates = counts[counts["n"] >= args.min_n]
    print(f"\n=== n >= {args.min_n} 인 종 수준 pangenome 분석 후보 ({len(candidates)}개) ===")
    print(candidates.to_string(index=False))

    if len(candidates) > 0:
        print("\n추천: 위 후보 중 LAB 1개 + Bacillus-group 1개를 선택하면 두 기능군을 대표하는")
        print("종 수준 pangenome 예시를 보여줄 수 있습니다.")


if __name__ == "__main__":
    main()
