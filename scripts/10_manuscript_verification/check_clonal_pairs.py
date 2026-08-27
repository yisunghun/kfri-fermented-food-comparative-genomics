#!/usr/bin/env python3
"""
check_clonal_pairs.py

이미 계산되어 있는 220x220 ANI 행렬(ani_matrix.csv)을 이용해, ANI가 매우
높은(=사실상 동일 균주로 의심되는) 시료 쌍을 찾는다. 논문 통계(LAB vs
Bacillus-group 등)에서 이런 쌍이 존재하면 통계적 독립성 가정이 깨질 수 있음
(pseudoreplication) - 리뷰 지적사항 대응용.

사용법:
    python3 check_clonal_pairs.py \
        --ani-matrix ani_analysis/ani_matrix.csv \
        --master master_table_qc.tsv \
        --outdir clonal_check \
        --thresholds 99.9 99.5 99.0
"""
import argparse
import os
import pandas as pd
import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ani-matrix", required=True)
    ap.add_argument("--master", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--thresholds", type=float, nargs="+", default=[99.9, 99.5, 99.0])
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    ani = pd.read_csv(args.ani_matrix, index_col=0)
    ani.index = ani.index.astype(str)
    ani.columns = ani.columns.astype(str)

    master = pd.read_csv(args.master, sep="\t")
    master = master.set_index("sample_id")

    # 상삼각행렬만 사용 (자기자신 제외, 중복쌍 제외)
    n = len(ani)
    mask = np.triu(np.ones((n, n), dtype=bool), k=1)
    pairs = []
    ani_vals = ani.values
    idx = ani.index.tolist()
    for i in range(n):
        for j in range(n):
            if not mask[i, j]:
                continue
            v = ani_vals[i, j]
            if pd.isna(v):
                continue
            pairs.append((idx[i], idx[j], v))

    pairs_df = pd.DataFrame(pairs, columns=["sample_a", "sample_b", "ani"])
    pairs_df = pairs_df.sort_values("ani", ascending=False)

    def annotate(df):
        df = df.copy()
        for col, src in [("genus_a", "genus_final"), ("species_a", "species_final"),
                          ("functional_group_a", "functional_group")]:
            df[col] = df["sample_a"].map(master[src]) if src in master.columns else np.nan
        for col, src in [("genus_b", "genus_final"), ("species_b", "species_final"),
                          ("functional_group_b", "functional_group")]:
            df[col] = df["sample_b"].map(master[src]) if src in master.columns else np.nan
        return df

    pairs_df = annotate(pairs_df)
    pairs_df.to_csv(os.path.join(args.outdir, "all_pairs_ani_sorted.csv"), index=False)

    print("=== 임계값별 초근접(잠재적 클론) 쌍 개수 ===")
    for th in sorted(args.thresholds, reverse=True):
        n_pairs = (pairs_df["ani"] >= th).sum()
        involved = set(pairs_df.loc[pairs_df["ani"] >= th, "sample_a"]) | \
                   set(pairs_df.loc[pairs_df["ani"] >= th, "sample_b"])
        print(f"  ANI >= {th}%: {n_pairs}쌍 (관련 시료 {len(involved)}개)")

    top_th = min(args.thresholds)
    close_pairs = pairs_df[pairs_df["ani"] >= top_th].copy()
    out_path = os.path.join(args.outdir, f"close_pairs_ge{top_th}.csv")
    close_pairs.to_csv(out_path, index=False)

    print(f"\n=== ANI >= {top_th}% 상세 목록 (상위 30개) ===")
    cols = ["sample_a", "sample_b", "ani", "genus_a", "species_a", "genus_b", "species_b"]
    print(close_pairs[cols].head(30).to_string(index=False))

    print(f"\n전체 상세 목록: {out_path}")
    print(f"전체 쌍(정렬됨): {os.path.join(args.outdir, 'all_pairs_ani_sorted.csv')}")

    # 가장 엄격한 임계값(99.9%) 기준으로 커넥티드 컴포넌트(클론 그룹) 찾기
    strictest = max(args.thresholds)
    strict_pairs = pairs_df[pairs_df["ani"] >= strictest]
    if len(strict_pairs) > 0:
        import networkx as nx
        G = nx.Graph()
        G.add_edges_from(zip(strict_pairs["sample_a"], strict_pairs["sample_b"]))
        components = [c for c in nx.connected_components(G) if len(c) > 1]
        print(f"\n=== ANI >= {strictest}% 기준 클론 그룹(연결된 시료 묶음): {len(components)}개 ===")
        for i, comp in enumerate(components, 1):
            print(f"  그룹 {i}: {sorted(comp)}")
            for s in comp:
                if s in master.index:
                    print(f"    {s}: {master.loc[s].get('species_final', '?')}, "
                          f"functional_group={master.loc[s].get('functional_group', '?')}")


if __name__ == "__main__":
    main()
