#!/usr/bin/env python3
"""
analyze_ani.py

fastANI all-vs-all 결과(ani_result.tsv)와 master_table.tsv를 합쳐서:
  1) 대칭 ANI 매트릭스 생성 (CSV로 저장)
  2) 계층적 클러스터링 덴드로그램 (속/genus 색상 표시, PDF로 저장)
  3) BLAST 기반 species_guess와 ANI 기반(>=95% = 동종) 클러스터링 결과 불일치 지점 리포트
     - 같은 species_guess인데 ANI가 종 기준(95%)보다 낮은 경우
     - species_guess는 다른데 ANI가 종 기준 이상으로 매우 가까운 경우

사용법:
  python3 analyze_ani.py \
      --ani ani_result.tsv \
      --master master_table.tsv \
      --outdir ./ani_analysis

필요 패키지: pip install pandas scipy matplotlib numpy
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform

SPECIES_ANI_THRESHOLD = 95.0  # 통상적으로 사용되는 종(species) 구분 기준 ANI


def setup_korean_font():
    """설치된 한글 폰트를 찾아 matplotlib에 지정. 없으면 경고만 출력하고 계속 진행."""
    import matplotlib.font_manager as fm

    candidates = ["NanumGothic", "Nanum Gothic", "Malgun Gothic", "AppleGothic",
                  "Noto Sans CJK KR", "Noto Sans KR", "UnDotum", "Batang"]
    available = {f.name for f in fm.fontManager.ttflist}
    for c in candidates:
        if c in available:
            plt.rcParams["font.family"] = c
            plt.rcParams["axes.unicode_minus"] = False
            print(f"[한글 폰트] '{c}' 사용", file=sys.stderr)
            return c
    print("[경고] 한글 지원 폰트를 찾지 못했습니다. 그래프의 한글이 깨져 보일 수 있습니다.\n"
          "       설치: sudo apt-get install -y fonts-nanum && rm -rf ~/.cache/matplotlib",
          file=sys.stderr)
    return None


def sample_id_from_fna_path(path: str) -> str:
    return os.path.splitext(os.path.basename(path))[0]


def build_matrix(ani_tsv: str, sample_ids: list[str]) -> pd.DataFrame:
    cols = ["query", "reference", "ani", "frag_matched", "frag_total"]
    df = pd.read_csv(ani_tsv, sep="\t", header=None, names=cols)
    df["query"] = df["query"].apply(sample_id_from_fna_path)
    df["reference"] = df["reference"].apply(sample_id_from_fna_path)

    mat = pd.DataFrame(np.nan, index=sample_ids, columns=sample_ids, dtype=float)
    for _, row in df.iterrows():
        if row["query"] in mat.index and row["reference"] in mat.columns:
            mat.loc[row["query"], row["reference"]] = row["ani"]

    # fastANI는 비대칭 결과를 줄 수 있음(양방향 중 하나만 threshold 넘을 때)
    # -> query->ref, ref->query 두 값의 nanmean으로 대칭화 (둘 다 NaN이면 NaN 유지)
    stacked = np.stack([mat.to_numpy(), mat.to_numpy().T])
    with np.errstate(invalid="ignore"):
        combined = np.nanmean(stacked, axis=0)
    np.fill_diagonal(combined, 100.0)
    return pd.DataFrame(combined, index=sample_ids, columns=sample_ids)


def plot_dendrogram(ani_matrix: pd.DataFrame, labels_df: pd.DataFrame, out_pdf: str):
    # 거리 = 100 - ANI (fastANI 미검출 쌍은 80% 미만으로 간주해 낮은 값으로 채움 -> 먼 거리)
    dist_arr = (100.0 - ani_matrix.fillna(70.0)).to_numpy(copy=True)
    dist_arr = np.clip(dist_arr, 0.0, None)  # 부동소수점 오차로 인한 미세 음수 방지
    np.fill_diagonal(dist_arr, 0.0)
    condensed = squareform(dist_arr, checks=False)
    Z = linkage(condensed, method="average")

    genus_map = labels_df.set_index("sample_id")["genus"].to_dict()
    unique_genera = sorted(set(genus_map.values()))
    cmap = plt.colormaps["tab20"].resampled(max(len(unique_genera), 1))
    color_map = {g: cmap(i) for i, g in enumerate(unique_genera)}

    fig, ax = plt.subplots(figsize=(14, max(8, len(ani_matrix) * 0.12)))
    dendrogram(
        Z,
        labels=ani_matrix.index.tolist(),
        orientation="left",
        ax=ax,
        leaf_font_size=6,
    )
    # y축 라벨에 genus 색상 입히기
    for tick_label in ax.get_ymajorticklabels():
        sid = tick_label.get_text()
        genus = genus_map.get(sid, "unknown")
        tick_label.set_color(color_map.get(genus, "black"))

    ax.set_title("ANI-based Hierarchical Clustering (label color = species_guess genus)")
    ax.set_xlabel("Distance (100 - ANI)")
    fig.tight_layout()
    fig.savefig(out_pdf)
    print(f"덴드로그램 저장: {out_pdf}")


def cross_check(ani_matrix: pd.DataFrame, labels_df: pd.DataFrame, out_csv: str):
    labels = labels_df.set_index("sample_id")["species_guess"].to_dict()
    rows = []
    ids = ani_matrix.index.tolist()
    for i in range(len(ids)):
        for j in range(i + 1, len(ids)):
            a, b = ids[i], ids[j]
            ani = ani_matrix.loc[a, b]
            if pd.isna(ani):
                continue
            sp_a, sp_b = labels.get(a, "unknown"), labels.get(b, "unknown")
            same_species_label = (sp_a == sp_b) and sp_a != "unknown"
            same_species_ani = ani >= SPECIES_ANI_THRESHOLD
            if same_species_label != same_species_ani:
                rows.append({
                    "sample_a": a, "sample_b": b, "ANI": round(ani, 2),
                    "species_guess_a": sp_a, "species_guess_b": sp_b,
                    "issue": ("same species_guess but ANI below species threshold"
                              if same_species_label else
                              "different species_guess but ANI above species threshold (very close relatives)"),
                })
    out_df = pd.DataFrame(rows).sort_values("ANI", ascending=False) if rows else pd.DataFrame(
        columns=["sample_a", "sample_b", "ANI", "species_guess_a", "species_guess_b", "issue"])
    # .tsv를 엑셀에서 더블클릭으로 열면 탭 구분이 자동 인식 안 되는 경우가 많아 .csv(콤마구분)로 저장
    out_df.to_csv(out_csv, index=False, encoding="utf-8-sig")
    print(f"불일치 리포트({len(out_df)}건): {out_csv}")
    return out_df


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--ani", required=True, help="fastANI 결과 tsv (run_fastani.sh 출력)")
    ap.add_argument("--master", required=True, help="master_table.tsv 경로")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--species-col", default="species_guess",
                     help="genus/species 라벨로 사용할 컬럼명 (정규화 후엔 species_normalized 지정)")
    args = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # noqa: E402  (Agg 백엔드 설정 후 import)
    globals()["plt"] = plt
    setup_korean_font()

    os.makedirs(args.outdir, exist_ok=True)

    master = pd.read_csv(args.master, sep="\t")
    if args.species_col not in master.columns:
        raise SystemExit(f"오류: '{args.species_col}' 컬럼이 master 테이블에 없습니다. "
                          f"사용 가능한 컬럼: {list(master.columns)}")
    master["genus"] = master[args.species_col].apply(
        lambda s: str(s).split()[0] if isinstance(s, str) and s and not s.startswith("unknown") else "unknown"
    )
    master["species_guess"] = master[args.species_col]  # 아래 함수들과의 호환을 위해 통일

    sample_ids = master["sample_id"].tolist()
    print(f"master_table 기준 {len(sample_ids)}개 시료")

    ani_matrix = build_matrix(args.ani, sample_ids)
    matrix_csv = os.path.join(args.outdir, "ani_matrix.csv")
    ani_matrix.to_csv(matrix_csv, encoding="utf-8-sig")
    print(f"ANI 매트릭스 저장: {matrix_csv}")

    plot_dendrogram(ani_matrix, master, os.path.join(args.outdir, "ani_dendrogram.pdf"))
    cross_check(ani_matrix, master, os.path.join(args.outdir, "ani_species_mismatch.csv"))
