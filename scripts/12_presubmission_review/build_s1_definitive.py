#!/usr/bin/env python3
"""
build_s1_definitive.py

최종 정의판: 220개 전체를 GTDB-Tk 결과로 채우되, 본문/Table S4에 이미 구체적
ANI 수치와 함께 인용된 12건(Part A 7건 + Part B 5건)은 해당 문서와 정확히
일치하는 값으로 덮어써서(override), 본문과 Table S1 간 불일치를 방지한다.

사용법:
    python3 build_s1_definitive.py \
        --s1-tsv Supplementary_Table_S1_genome_QC_metrics.tsv \
        --gtdbtk-summary /mnt/f/WGS_Consolidated/gtdbtk_result_241/classify/gtdbtk.bac120.summary.tsv \
        --outdir .
"""
import argparse
import re
import pandas as pd

# 본문/Table S4에 구체적 ANI 수치와 함께 이미 인용된 12건 - 반드시 이 값으로 고정
FIXED_OVERRIDES = {
    "HN00167634_C4-2": ("Bacillus", "Bacillus paranthracis"),
    "HN00171167_F3062": ("Bacillus", "Bacillus aerophilus"),
    "HN00180670_CHKJ1127": ("Enterococcus", "Enterococcus lactis"),
    "HN00180670_CHKJ1223-2": ("Enterococcus", "Enterococcus lactis"),
    "HN00180670_MK1-12": ("Enterococcus", "Enterococcus lactis"),
    "HN00196390_CHKJ1122": ("Enterococcus", "Enterococcus lactis"),
    "HN00251139_BTH25001": ("Enterococcus", "Enterococcus lactis"),
    "HN00251139_BJS25005": ("Thalassorhabdomicrobium", "Thalassorhabdomicrobium marinisediminis"),
    "HN00251139_BJS25013": ("Marihabitans", "Marihabitans asiaticum"),
    "HN00251139_BJS25032": ("Marihabitans", "Marihabitans asiaticum"),
    "HN00251139_KJS25025": ("Nesterenkonia", "Nesterenkonia koreensis"),
    "HN00280011_C3-36": ("Kocuria", "Kocuria atrinae"),
}


def normalize_gtdb_token(token):
    """GTDB 다계통 분리 표기(_A, _B 등)를 NCBI/ICNP 스타일로 정규화 (genus, species 공통)"""
    return re.sub(r'_[A-Z]+$', '', token)


def parse_gtdb_classification(classification_str):
    parts = dict(p.split('__', 1) for p in classification_str.split(';') if '__' in p)
    genus_raw = parts.get('g', '')
    species_raw = parts.get('s', '')

    genus_norm = normalize_gtdb_token(genus_raw) if genus_raw else None

    species_norm = None
    if species_raw:
        sp_tokens = species_raw.split(' ', 1)
        if len(sp_tokens) == 2:
            epithet_norm = normalize_gtdb_token(sp_tokens[1])
            species_norm = f"{genus_norm} {epithet_norm}"

    return genus_norm, species_norm


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--s1-tsv", required=True)
    ap.add_argument("--gtdbtk-summary", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    s1 = pd.read_csv(args.s1_tsv, sep="\t")
    gtdb = pd.read_csv(args.gtdbtk_summary, sep="\t")

    gtdb_map = {}
    for _, row in gtdb.iterrows():
        sid = row['user_genome']
        genus, species = parse_gtdb_classification(row['classification'])
        gtdb_map[sid] = (genus, species)

    n_species_level = 0
    n_genus_only = 0
    n_override = 0

    new_genus_col = []
    new_species_col = []
    for _, row in s1.iterrows():
        sid = row['sample_id']

        if sid in FIXED_OVERRIDES:
            genus, species = FIXED_OVERRIDES[sid]
            new_genus_col.append(genus)
            new_species_col.append(species)
            n_species_level += 1
            n_override += 1
            continue

        genus, species = gtdb_map.get(sid, (None, None))
        has_real_species = bool(species) and len(species.split()) >= 2
        if has_real_species:
            new_genus_col.append(genus)
            new_species_col.append(species)
            n_species_level += 1
        elif genus:
            new_genus_col.append(genus)
            new_species_col.append(f"{genus} sp.")
            n_genus_only += 1
        else:
            new_genus_col.append('unresolved')
            new_species_col.append('')
            n_genus_only += 1

    s1['genus_final_validated'] = new_genus_col
    s1['species_final_validated'] = new_species_col

    print(f"Species-level 해결: {n_species_level}/{len(s1)} (이 중 본문 인용 12건 override 적용: {n_override})")
    print(f"Genus-level만 해결: {n_genus_only}/{len(s1)}")

    from collections import Counter
    genus_counts = Counter(new_genus_col)
    print(f"\nEnterococcus: {genus_counts.get('Enterococcus', 0)}")
    print(f"Lactococcus: {genus_counts.get('Lactococcus', 0)}")
    named = [g for g in genus_counts if g != 'unresolved']
    print(f"고유 속 개수(unresolved 제외): {len(named)}")

    # 12건 최종 검증 출력
    print("\n=== 12건 override 최종 검증 ===")
    for sid, (g, sp) in FIXED_OVERRIDES.items():
        actual_row = s1[s1['sample_id'] == sid]
        actual_sp = actual_row.iloc[0]['species_final_validated']
        status = "OK" if actual_sp == sp else "!!"
        print(f"{sid}: {actual_sp} [{status}]")

    out_path = f"{args.outdir}/Supplementary_Table_S1_genome_QC_metrics_DEFINITIVE.tsv"
    s1.to_csv(out_path, sep="\t", index=False)
    print(f"\n저장 완료: {out_path}")


if __name__ == "__main__":
    main()
