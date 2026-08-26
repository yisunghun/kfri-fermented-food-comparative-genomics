#!/usr/bin/env python3
"""
resolve_species_via_ani.py

BLAST 기반 species_normalized 중 종 수준까지 판정되지 않은 것들
(예: 'Lactobacillus sp.', 'Bacillus sp.', 'unknown(...)')을
ANI 매트릭스 상의 근연 시료(>= ani_threshold, 기본 95%)를 근거로
종명을 확정합니다.

판정 로직 (샘플 X가 애매한 경우):
  1) ANI 매트릭스에서 X와 ani_threshold 이상인 다른 시료들을 찾음
  2) 그 중 species_normalized가 '명확한'(=애매하지 않은) 시료들의 종명을 수집
  3) 종명이 1개로 수렴 -> species_final = 그 종명, resolution_method = ANI_neighbor
     종명이 여러 개로 갈림 -> species_final = 원래값 유지, resolution_method = ANI_conflict (수동 검토 필요)
     명확한 이웃이 아예 없음 -> species_final = 원래값 유지, resolution_method = unresolved
  4) 원래부터 명확했던 시료는 species_final = species_normalized, resolution_method = original_call

'애매함' 판정: species가 'Genus species' 두 단어 형태가 아니거나,
두번째 단어가 'sp'/'sp.'이거나, 'unknown'으로 시작하는 경우.

의존성: 같은 폴더의 analyze_ani.py (build_matrix, sample_id_from_fna_path 재사용)

사용법:
  python3 resolve_species_via_ani.py \
      --ani ani_result.tsv \
      --master master_table_normalized.tsv \
      --species-col species_normalized \
      --out master_table_final.tsv \
      --ani-threshold 95.0
"""
import argparse
import sys

import numpy as np
import pandas as pd

try:
    from analyze_ani import build_matrix, sample_id_from_fna_path  # noqa: F401
except ImportError:
    print("오류: analyze_ani.py를 같은 폴더에서 찾을 수 없습니다. "
          "resolve_species_via_ani.py와 같은 디렉토리에 두세요.", file=sys.stderr)
    raise


def is_ambiguous(species: str) -> bool:
    if not isinstance(species, str) or species.strip() == "":
        return True
    if species.startswith("unknown"):
        return True
    parts = species.split()
    if len(parts) < 2:
        return True
    epithet = parts[1].rstrip(".").lower()
    if epithet == "sp":
        return True
    return False


def resolve(master: pd.DataFrame, ani_matrix: pd.DataFrame, species_col: str, ani_threshold: float):
    species_map = master.set_index("sample_id")[species_col].to_dict()

    final_species, methods, neighbor_ids, neighbor_anis = [], [], [], []

    for sid in master["sample_id"]:
        sp = species_map.get(sid)
        if not is_ambiguous(sp):
            final_species.append(sp)
            methods.append("original_call")
            neighbor_ids.append("")
            neighbor_anis.append(np.nan)
            continue

        if sid not in ani_matrix.index:
            final_species.append(sp)
            methods.append("unresolved(no_ani_row)")
            neighbor_ids.append("")
            neighbor_anis.append(np.nan)
            continue

        row = ani_matrix.loc[sid].drop(labels=[sid], errors="ignore")
        close_neighbors = row[row >= ani_threshold].sort_values(ascending=False)

        candidates = {}  # resolved_species -> (best_ani, neighbor_sid)
        for neighbor_sid, ani_val in close_neighbors.items():
            neighbor_sp = species_map.get(neighbor_sid)
            if neighbor_sp and not is_ambiguous(neighbor_sp):
                if neighbor_sp not in candidates:
                    candidates[neighbor_sp] = (ani_val, neighbor_sid)

        if len(candidates) == 0:
            final_species.append(sp)
            methods.append("unresolved")
            neighbor_ids.append("")
            neighbor_anis.append(np.nan)
        elif len(candidates) == 1:
            resolved_sp, (best_ani, best_sid) = next(iter(candidates.items()))
            final_species.append(resolved_sp)
            methods.append("ANI_neighbor")
            neighbor_ids.append(best_sid)
            neighbor_anis.append(best_ani)
        else:
            # 후보가 여러 종으로 갈림 -> 원래값 유지하고 충돌로 표시, 최고 ANI 이웃 정보만 참고용으로 기록
            best_sp = max(candidates.items(), key=lambda kv: kv[1][0])
            final_species.append(sp)
            methods.append(f"ANI_conflict({','.join(candidates.keys())})")
            neighbor_ids.append(best_sp[1][1])
            neighbor_anis.append(best_sp[1][0])

    master = master.copy()
    master["species_final"] = final_species
    master["resolution_method"] = methods
    master["ani_neighbor_sample"] = neighbor_ids
    master["ani_neighbor_value"] = neighbor_anis
    master["genus_final"] = master["species_final"].apply(
        lambda s: s.split()[0] if isinstance(s, str) and s and not is_ambiguous(s) else "unresolved"
    )
    return master


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--ani", required=True)
    ap.add_argument("--master", required=True)
    ap.add_argument("--species-col", default="species_normalized")
    ap.add_argument("--out", required=True)
    ap.add_argument("--ani-threshold", type=float, default=95.0)
    args = ap.parse_args()

    master = pd.read_csv(args.master, sep="\t")
    sample_ids = master["sample_id"].tolist()
    ani_matrix = build_matrix(args.ani, sample_ids)

    n_ambiguous_before = master[args.species_col].apply(is_ambiguous).sum()
    print(f"보정 전 애매한 판정: {n_ambiguous_before}건 / 총 {len(master)}건")

    result = resolve(master, ani_matrix, args.species_col, args.ani_threshold)
    result.to_csv(args.out, sep="\t", index=False, encoding="utf-8-sig")

    print(f"\n[resolution_method 분포]")
    print(result["resolution_method"].apply(
        lambda m: m.split("(")[0]
    ).value_counts().to_string())

    n_ambiguous_after = result["species_final"].apply(is_ambiguous).sum()
    print(f"\n보정 후에도 애매한 판정(unresolved/conflict): {n_ambiguous_after}건")
    print(f"저장 완료: {args.out}")

    conflicts = result[result["resolution_method"].str.startswith("ANI_conflict", na=False)]
    if not conflicts.empty:
        print(f"\n[충돌(ANI_conflict) 상세 - 수동 검토 권장]")
        print(conflicts[["sample_id", args.species_col, "resolution_method"]].to_string(index=False))
