#!/usr/bin/env python3
"""
select_dereplicated_representatives.py

ani_matrix.csv에서 지정한 임계값(기본 99.9%) 이상으로 연결된 시료들을
"클론 그룹"으로 묶고, 그룹당 대표 시료 1개만 남긴 목록을 만든다.
(단일 시료로만 이루어진 "그룹"은 그대로 전부 유지됨)

사용법:
    python3 select_dereplicated_representatives.py \
        --ani-matrix ani_analysis/ani_matrix.csv \
        --threshold 99.9 \
        --outdir dereplication
"""
import argparse
import os
import numpy as np
import pandas as pd
import networkx as nx


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ani-matrix", required=True)
    ap.add_argument("--threshold", type=float, default=99.9)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    ani = pd.read_csv(args.ani_matrix, index_col=0)
    ani.index = ani.index.astype(str)
    ani.columns = ani.columns.astype(str)
    all_samples = list(ani.index)

    G = nx.Graph()
    G.add_nodes_from(all_samples)  # 모든 시료를 일단 노드로 (짝 없는 시료도 그룹 크기 1로 유지)

    n = len(ani)
    vals = ani.values
    for i in range(n):
        for j in range(i + 1, n):
            v = vals[i, j]
            if pd.notna(v) and v >= args.threshold:
                G.add_edge(all_samples[i], all_samples[j])

    components = list(nx.connected_components(G))
    print(f"임계값 ANI >= {args.threshold}% 기준: 전체 {len(all_samples)}개 시료 -> {len(components)}개 그룹")

    rows = []
    representatives = []
    for gid, comp in enumerate(sorted(components, key=lambda c: -len(c)), 1):
        comp_sorted = sorted(comp)
        rep = comp_sorted[0]  # 알파벳/문자열 순 첫 번째를 대표로 (재현 가능하도록 결정론적 선택)
        representatives.append(rep)
        for s in comp_sorted:
            rows.append({"sample_id": s, "clone_group_id": gid, "group_size": len(comp_sorted),
                         "is_representative": s == rep})

    groups_df = pd.DataFrame(rows)
    groups_df.to_csv(os.path.join(args.outdir, "clone_groups.tsv"), sep="\t", index=False)

    with open(os.path.join(args.outdir, "representative_samples.txt"), "w") as f:
        for s in representatives:
            f.write(s + "\n")

    n_multi = sum(1 for c in components if len(c) > 1)
    n_removed = len(all_samples) - len(representatives)
    print(f"다중 시료 그룹(클론 그룹): {n_multi}개")
    print(f"대표 시료로 축소: {len(all_samples)}개 -> {len(representatives)}개 (제거됨: {n_removed}개)")
    print(f"\n결과 파일:")
    print(f"  {os.path.join(args.outdir, 'clone_groups.tsv')}")
    print(f"  {os.path.join(args.outdir, 'representative_samples.txt')}")


if __name__ == "__main__":
    main()
