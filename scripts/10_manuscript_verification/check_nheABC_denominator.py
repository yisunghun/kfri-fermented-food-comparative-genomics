#!/usr/bin/env python3
"""
check_nheABC_denominator.py

nheA/nheB/nheC(Bacillus cereus non-hemolytic enterotoxin)를 모두 가진 시료가
1) 전체 Bacillus-group 중 몇 %인지 (지금까지 원고에 쓴 방식)
2) 진짜 의미 있는 분모인 "B. cereus sensu lato(종복합체)" 내에서는 몇 %인지
를 정확히 교차확인한다. 또한 nheABC 양성이면서 B. cereus 종복합체 밖에 있는
시료가 있는지도 확인한다 (있다면 흥미로운 별도 소견).

사용법:
    python3 check_nheABC_denominator.py \
        --vfdb-summary abricate_out/vfdb_summary.tsv \
        --master master_table_qc.tsv \
        --functional-group-root grouped/by_functional_group
"""
import argparse
import glob
import os
import re
import pandas as pd

# B. cereus sensu lato (종복합체) 소속 종 목록 (Liu et al. 2015; Gupta et al. 2020 기준
# 통용되는 species complex 구성원)
B_CEREUS_SENSU_LATO = {
    "cereus", "thuringiensis", "anthracis", "mycoides", "pseudomycoides",
    "weihenstephanensis", "toyonensis", "cytotoxicus", "paranthracis",
    "albus", "mobilis", "proteolyticus", "tropicus", "wiedmannii",
    "paramycoides", "luti", "nitratireducens", "pacificus",
}


def derive_sample_id(raw):
    base = os.path.basename(str(raw))
    return re.sub(r"\.(tab|fna|fa|fasta|txt|csv|tsv)$", "", base)


def load_functional_group_map(fg_root):
    mapping = {}
    for fg_name in os.listdir(fg_root):
        fg_dir = os.path.join(fg_root, fg_name)
        if not os.path.isdir(fg_dir):
            continue
        gff_dir = os.path.join(fg_dir, "genomes_gff")
        search_dir = gff_dir if os.path.isdir(gff_dir) else fg_dir
        for fpath in glob.glob(os.path.join(search_dir, "*.gff")):
            sid = os.path.splitext(os.path.basename(fpath))[0]
            mapping[sid] = fg_name
    return mapping


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vfdb-summary", required=True)
    ap.add_argument("--master", required=True)
    ap.add_argument("--functional-group-root", required=True)
    args = ap.parse_args()

    vfdb = pd.read_csv(args.vfdb_summary, sep="\t")
    vfdb["sample_id"] = vfdb[vfdb.columns[0]].apply(derive_sample_id)
    vfdb = vfdb.set_index("sample_id")

    # nheA/nheB/nheC 관련 컬럼 자동 탐색 (대소문자/표기 편차 대응)
    nhe_cols = [c for c in vfdb.columns if re.match(r"^nhe[ABC]$", c, re.IGNORECASE)]
    print(f"발견된 nhe 관련 컬럼: {nhe_cols}")
    if len(nhe_cols) < 3:
        print("[경고] nheA/nheB/nheC 3개를 모두 못 찾음. 컬럼명을 확인하세요:")
        print([c for c in vfdb.columns if "nhe" in c.lower()])
        return

    def has_gene(val):
        return str(val).strip() not in {".", "nan", ""}

    vfdb["nheABC_positive"] = vfdb[nhe_cols].apply(
        lambda row: all(has_gene(row[c]) for c in nhe_cols), axis=1
    )

    master = pd.read_csv(args.master, sep="\t")
    master.columns = [re.sub(r"^\ufeff", "", c) for c in master.columns]
    master = master.set_index("sample_id")

    fg_map = load_functional_group_map(args.functional_group_root)

    vfdb["functional_group"] = vfdb.index.map(lambda s: fg_map.get(s))
    vfdb["species_final"] = vfdb.index.map(lambda s: master.loc[s, "species_final"] if s in master.index else None)
    vfdb["species_epithet"] = vfdb["species_final"].apply(
        lambda s: str(s).split()[1] if isinstance(s, str) and len(str(s).split()) > 1 else None
    )
    vfdb["is_b_cereus_sl"] = vfdb["species_epithet"].apply(lambda e: e in B_CEREUS_SENSU_LATO)

    bacillus_group = vfdb[vfdb["functional_group"] == "Bacillus_group"]
    n_bacillus_group = len(bacillus_group)
    n_nhe_pos_in_bacillus_group = bacillus_group["nheABC_positive"].sum()

    b_cereus_sl = vfdb[vfdb["is_b_cereus_sl"]]
    n_b_cereus_sl = len(b_cereus_sl)
    n_nhe_pos_in_b_cereus_sl = b_cereus_sl["nheABC_positive"].sum()

    print(f"\n=== 분모 1: 전체 Bacillus-group ===")
    print(f"nheABC 양성: {n_nhe_pos_in_bacillus_group} / {n_bacillus_group} "
          f"({100*n_nhe_pos_in_bacillus_group/n_bacillus_group:.1f}%)")

    print(f"\n=== 분모 2: B. cereus sensu lato (종복합체)만 ===")
    print(f"B. cereus 종복합체 시료 목록:")
    print(b_cereus_sl[["species_final"]].to_string())
    print(f"\nnheABC 양성: {n_nhe_pos_in_b_cereus_sl} / {n_b_cereus_sl} "
          f"({100*n_nhe_pos_in_b_cereus_sl/n_b_cereus_sl:.1f}%)" if n_b_cereus_sl > 0 else "B. cereus 종복합체 시료 없음")

    print(f"\n=== nheABC 양성인데 B. cereus 종복합체 밖인 시료 (있으면 특이 소견) ===")
    outside = vfdb[vfdb["nheABC_positive"] & ~vfdb["is_b_cereus_sl"]]
    if len(outside) > 0:
        print(outside[["species_final", "functional_group"]].to_string())
    else:
        print("없음 (nheABC 양성 시료는 전부 B. cereus 종복합체 소속)")


if __name__ == "__main__":
    main()
