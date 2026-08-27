#!/usr/bin/env python3
"""
resolve_ani_mismatches.py

reference_ani_validation.tsv에서 '원래 배정된 종'과의 ANI가 낮게 나온(또는
계산조차 안 된) 시료들을, 이번엔 원래 배정된 종 하나만이 아니라 확보해둔
"72개 reference genome 전체"와 비교하여 실제로 가장 가까운 종이 무엇인지
찾는다. genus 단위 재배정까지 필요한지 판단하기 위한 스크립트.

사용법:
    python3 resolve_ani_mismatches.py \
        --validation-tsv reference_ani_validation.tsv \
        --genomes-dir genomes_fna \
        --reference-dir reference_genomes/fasta \
        --outdir mismatch_resolution \
        --ani-threshold 95 \
        --threads 24
"""
import argparse
import glob
import os
import subprocess
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--validation-tsv", required=True)
    ap.add_argument("--genomes-dir", required=True)
    ap.add_argument("--reference-dir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--ani-threshold", type=float, default=95.0)
    ap.add_argument("--threads", type=int, default=24)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    val = pd.read_csv(args.validation_tsv, sep="\t")

    # 재검토 대상 = (OK인데 ANI<threshold) 또는 (BELOW_ANI_DETECTION_THRESHOLD)
    # UNRESOLVED_SPECIES_SKIPPED(원래부터 genus만 판정)는 제외 - 이미 정상 처리됨
    problem = val[
        ((val["status"] == "OK") & (val["ani"] < args.ani_threshold))
        | (val["status"] == "BELOW_ANI_DETECTION_THRESHOLD")
    ].copy()

    print(f"재검토 대상: {len(problem)}개")

    # 쿼리 리스트 파일
    query_list_path = os.path.join(args.outdir, "query_list.txt")
    with open(query_list_path, "w") as f:
        n_found = 0
        for sid in problem["sample_id"]:
            fpath = os.path.join(args.genomes_dir, f"{sid}.fna")
            if os.path.exists(fpath):
                f.write(fpath + "\n")
                n_found += 1
    print(f"쿼리 게놈 확인됨: {n_found}/{len(problem)}")

    # 레퍼런스 리스트 파일 (72개 전체)
    ref_list_path = os.path.join(args.outdir, "ref_list.txt")
    ref_files = sorted(glob.glob(os.path.join(args.reference_dir, "*.fna")))
    with open(ref_list_path, "w") as f:
        for rf in ref_files:
            f.write(rf + "\n")
    print(f"비교 대상 reference: {len(ref_files)}개")

    # fastANI many-vs-many
    out_path = os.path.join(args.outdir, "all_vs_all.ani.tsv")
    cmd = [
        "fastANI", "--ql", query_list_path, "--rl", ref_list_path,
        "-o", out_path, "-t", str(args.threads),
    ]
    print("fastANI 실행 중 (쿼리 x 레퍼런스 전체 조합)...")
    subprocess.run(cmd, check=True)

    cols = ["query", "reference", "ani", "aligned_fragments", "total_fragments"]
    ani_df = pd.read_csv(out_path, sep="\t", header=None, names=cols)
    ani_df["sample_id"] = ani_df["query"].apply(lambda p: os.path.basename(p).replace(".fna", ""))
    ani_df["candidate_species"] = ani_df["reference"].apply(
        lambda p: os.path.basename(p).replace(".fna", "").replace("_", " ")
    )

    # 시료별 top3 후보 추출
    results = []
    for sid, grp in ani_df.groupby("sample_id"):
        grp_sorted = grp.sort_values("ani", ascending=False)
        top3 = grp_sorted.head(3)
        orig_row = problem[problem["sample_id"] == sid].iloc[0]
        for rank, (_, r) in enumerate(top3.iterrows(), 1):
            results.append({
                "sample_id": sid,
                "originally_assigned_species": orig_row["reference_species"],
                "originally_assigned_ani": orig_row.get("ani", None),
                "rank": rank,
                "candidate_species": r["candidate_species"],
                "candidate_ani": r["ani"],
            })

    out_df = pd.DataFrame(results)
    out_df.to_csv(os.path.join(args.outdir, "mismatch_resolution_top3.tsv"), sep="\t", index=False)

    # 요약: rank1 후보가 threshold를 넘는지
    rank1 = out_df[out_df["rank"] == 1].copy()
    rank1["resolved"] = rank1["candidate_ani"] >= args.ani_threshold
    rank1["genus_changed"] = rank1.apply(
        lambda r: r["candidate_species"].split()[0] != str(r["originally_assigned_species"]).split()[0],
        axis=1,
    )

    print("\n=== 요약 ===")
    print(f"72개 reference 중 최선의 매치로 ANI>={args.ani_threshold}% 재확인됨: {rank1['resolved'].sum()}개")
    print(f"72개 reference 어느 것과 비교해도 ANI<{args.ani_threshold}%(72개 목록 밖의 종일 가능성): {(~rank1['resolved']).sum()}개")
    print(f"원래 배정과 다른 속(genus)으로 최선 매치가 나온 경우: {rank1['genus_changed'].sum()}개")

    rank1.to_csv(os.path.join(args.outdir, "mismatch_resolution_summary.tsv"), sep="\t", index=False)
    print(f"\n상세: {os.path.join(args.outdir, 'mismatch_resolution_top3.tsv')}")
    print(f"요약: {os.path.join(args.outdir, 'mismatch_resolution_summary.tsv')}")

    print("\n[속(genus)이 바뀌는 것으로 보이는 시료 목록]")
    changed = rank1[rank1["genus_changed"]]
    if len(changed) > 0:
        print(changed[["sample_id", "originally_assigned_species", "candidate_species", "candidate_ani"]].to_string(index=False))
    else:
        print("없음")

    print("\n[72개 reference 어느 것과도 안 맞는(목록 밖 종으로 추정) 시료 목록]")
    unresolved = rank1[~rank1["resolved"]]
    if len(unresolved) > 0:
        print(unresolved[["sample_id", "originally_assigned_species", "candidate_species", "candidate_ani"]].to_string(index=False))
    else:
        print("없음")


if __name__ == "__main__":
    main()
