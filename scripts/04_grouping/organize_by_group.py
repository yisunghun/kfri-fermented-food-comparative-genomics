#!/usr/bin/env python3
"""
organize_by_group.py

master_table_final.tsv (species_final/genus_final 포함)를 기준으로
genomes_faa/*.faa, genomes_gff/*.gff 파일들을 두 가지 구조로 재편합니다.

  1) by_genus/<Genus>/genomes_faa/*.faa, genomes_gff/*.gff
     - 표본 수가 --min-genus-n 미만인 속은 by_genus/_minor_genera/ 로 몰아넣음
       (Roary/Panaroo는 그룹 내 시료가 너무 적으면 의미있는 core/accessory 구분이 어려움)

  2) by_functional_group/<Group>/genomes_faa/*.faa, genomes_gff/*.gff
     - LAB(유산균류): Lactiplantibacillus, Levilactobacillus, Latilactobacillus,
       Lactobacillus, Lactococcus, Leuconostoc, Weissella, Pediococcus,
       Enterococcus, Tetragenococcus 등
     - Bacillus_group(포자형성 발효균, 장류 관련): Bacillus, Paenibacillus,
       Oceanobacillus, Virgibacillus, Priestia, Rossellomorea, Shouchella, Halobacillus
     - Other_Environmental: 그 외 전부
     - Unresolved: genus_final이 'unresolved'인 시료

각 그룹 폴더에는 genomes_faa/*.faa, genomes_gff/*.gff와 함께
group_manifest.tsv(그 그룹에 속한 시료 목록)도 같이 저장됩니다.

사용법:
  python3 organize_by_group.py \
      --master master_table_final.tsv \
      --outdir /mnt/f/WGS_Consolidated/grouped \
      --min-genus-n 3 \
      --mode copy      # 또는 symlink (WSL/NTFS 심볼릭 링크 지원 시 용량 절약)
"""
import argparse
import os
import shutil

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


def place_file(src: str, dst: str, mode: str):
    if not src or not isinstance(src, str) or not os.path.isfile(src):
        return False
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    if os.path.exists(dst):
        return True
    if mode == "symlink":
        os.symlink(os.path.abspath(src), dst)
    else:
        shutil.copy2(src, dst)
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--master", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--min-genus-n", type=int, default=3)
    ap.add_argument("--mode", choices=["copy", "symlink"], default="copy")
    args = ap.parse_args()

    master = pd.read_csv(args.master, sep="\t")
    genus_counts = master["genus_final"].value_counts()

    genus_root = os.path.join(args.outdir, "by_genus")
    func_root = os.path.join(args.outdir, "by_functional_group")

    genus_manifest_rows = []
    func_manifest_rows = []

    n_missing_faa = 0
    n_missing_gff = 0

    for _, row in master.iterrows():
        sid = row["sample_id"]
        genus = row["genus_final"]
        genus_dir_name = genus if genus_counts.get(genus, 0) >= args.min_genus_n else "_minor_genera"
        fgroup = functional_group(genus)

        faa_src = row.get("merged_faa")
        gff_src = row.get("merged_gff")

        for root, group_name, manifest_rows in [
            (genus_root, genus_dir_name, genus_manifest_rows),
            (func_root, fgroup, func_manifest_rows),
        ]:
            faa_dst = os.path.join(root, group_name, "genomes_faa", f"{sid}.faa")
            gff_dst = os.path.join(root, group_name, "genomes_gff", f"{sid}.gff")
            ok_faa = place_file(faa_src, faa_dst, args.mode)
            ok_gff = place_file(gff_src, gff_dst, args.mode)
            manifest_rows.append({
                "sample_id": sid, "group": group_name, "genus_final": genus,
                "species_final": row.get("species_final"),
                "faa_ok": ok_faa, "gff_ok": ok_gff,
            })

        if not (faa_src and os.path.isfile(str(faa_src))):
            n_missing_faa += 1
        if not (gff_src and os.path.isfile(str(gff_src))):
            n_missing_gff += 1

    for root, rows in [(genus_root, genus_manifest_rows), (func_root, func_manifest_rows)]:
        df = pd.DataFrame(rows)
        for group_name, sub in df.groupby("group"):
            manifest_path = os.path.join(root, group_name, "group_manifest.tsv")
            sub.to_csv(manifest_path, sep="\t", index=False, encoding="utf-8-sig")

    print(f"완료. 결과 위치: {args.outdir}")
    if n_missing_faa or n_missing_gff:
        print(f"[경고] 원본 파일 누락 - faa {n_missing_faa}건, gff {n_missing_gff}건 "
              f"(genomes_faa/genomes_gff 경로 확인 필요)")

    print("\n[by_genus 그룹별 시료 수]")
    print(pd.DataFrame(genus_manifest_rows)["group"].value_counts().to_string())

    print("\n[by_functional_group 그룹별 시료 수]")
    print(pd.DataFrame(func_manifest_rows)["group"].value_counts().to_string())


if __name__ == "__main__":
    main()
