#!/usr/bin/env python3
"""
validate_species_via_reference_ani.py

기존 파이프라인의 species 검증은 "220개 study genome 서로간" ANI였기 때문에
독립적인 검증이 아니라는 리뷰 지적을 반영, 각 isolate를 자신이 최종 판정된
species의 실제 NCBI reference/representative genome과 1:1 fastANI 비교하여
"reference genome 대비 ANI >= 95%"라는 독립적 근거를 만든다.

요구사항:
    - fetch_reference_genomes.sh 를 먼저 실행해서 reference_genomes/fasta/ 에
      "Genus_species.fna" 형태의 참조 게놈들이 준비되어 있어야 함
    - conda env: compgenomics (fastANI 설치되어 있어야 함)

사용법:
    python3 validate_species_via_reference_ani.py \
        --master master_table_qc.tsv \
        --genomes-dir genomes_fna \
        --reference-dir reference_genomes/fasta \
        --outdir reference_ani_validation \
        --threads 8
"""
import argparse
import os
import subprocess
import sys
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed


def run_fastani_pair(sample_id, query_fasta, ref_fasta, ref_species, tmp_dir):
    """단일 isolate vs 단일 reference genome fastANI 실행, 결과 dict 반환"""
    out_path = os.path.join(tmp_dir, f"{sample_id}.ani.out")
    cmd = ["fastANI", "-q", query_fasta, "-r", ref_fasta, "-o", out_path]
    result = {
        "sample_id": sample_id,
        "reference_species": ref_species,
        "reference_fasta": os.path.basename(ref_fasta),
        "ani": None,
        "aligned_fragments": None,
        "total_fragments": None,
        "status": "NO_OUTPUT",
    }
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        if proc.returncode != 0:
            result["status"] = f"FASTANI_ERROR: {proc.stderr.strip()[:200]}"
            return result
        if os.path.exists(out_path) and os.path.getsize(out_path) > 0:
            with open(out_path) as f:
                line = f.readline().strip().split("\t")
            # fastANI 출력: query ref ANI aligned_frags total_frags
            result["ani"] = float(line[2])
            result["aligned_fragments"] = int(line[3])
            result["total_fragments"] = int(line[4])
            result["status"] = "OK"
        else:
            # ANI 계산 실패 (너무 먼 genome pair인 경우 fastANI가 결과를 안 냄)
            result["status"] = "BELOW_ANI_DETECTION_THRESHOLD"
    except subprocess.TimeoutExpired:
        result["status"] = "TIMEOUT"
    except Exception as e:
        result["status"] = f"EXCEPTION: {e}"
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--master", required=True, help="master_table_qc.tsv")
    ap.add_argument("--genomes-dir", required=True, help="isolate별 consolidated .fna 폴더 (파일명: <sample_id>.fna)")
    ap.add_argument("--reference-dir", required=True, help="fetch_reference_genomes.sh로 받은 reference fasta 폴더")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--ani-threshold", type=float, default=95.0)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    tmp_dir = os.path.join(args.outdir, "tmp_ani_out")
    os.makedirs(tmp_dir, exist_ok=True)

    master = pd.read_csv(args.master, sep="\t")
    if "sample_id" not in master.columns or "species_final" not in master.columns:
        sys.exit("master 파일에 sample_id, species_final 컬럼이 필요합니다.")

    jobs = []
    skipped = []
    for _, row in master.iterrows():
        sample_id = row["sample_id"]
        species = str(row["species_final"]).strip()

        if species.endswith("sp.") or species.lower().startswith("unknown") or species == "" or species == "nan":
            skipped.append((sample_id, species, "UNRESOLVED_SPECIES_SKIPPED"))
            continue

        query_fasta = os.path.join(args.genomes_dir, f"{sample_id}.fna")
        ref_fasta = os.path.join(args.reference_dir, species.replace(" ", "_") + ".fna")

        if not os.path.exists(query_fasta):
            skipped.append((sample_id, species, "QUERY_FASTA_MISSING"))
            continue
        if not os.path.exists(ref_fasta):
            skipped.append((sample_id, species, "REFERENCE_FASTA_MISSING"))
            continue

        jobs.append((sample_id, query_fasta, ref_fasta, species))

    print(f"총 {len(master)}개 시료 중 {len(jobs)}개 검증 실행, {len(skipped)}개 스킵")

    results = []
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        futures = {
            executor.submit(run_fastani_pair, sid, qf, rf, sp, tmp_dir): sid
            for sid, qf, rf, sp in jobs
        }
        for i, future in enumerate(as_completed(futures), 1):
            results.append(future.result())
            if i % 20 == 0:
                print(f"  진행: {i}/{len(jobs)}")

    df = pd.DataFrame(results)
    df["pass_95pct"] = df["ani"].apply(lambda x: (x is not None) and (x >= args.ani_threshold))

    skipped_df = pd.DataFrame(skipped, columns=["sample_id", "reference_species", "status"])
    skipped_df["ani"] = None
    skipped_df["pass_95pct"] = False

    full_df = pd.concat([df, skipped_df], ignore_index=True, sort=False)
    full_df.to_csv(os.path.join(args.outdir, "reference_ani_validation.tsv"), sep="\t", index=False)

    n_pass = full_df["pass_95pct"].sum()
    n_tested = df[df["status"] == "OK"].shape[0]
    n_fail_low_ani = df[(df["status"] == "OK") & (df["ani"] < args.ani_threshold)].shape[0]

    print("\n=== 요약 ===")
    print(f"검증 가능(reference 확보): {len(jobs)}개")
    print(f"  - fastANI 정상 계산: {n_tested}개")
    print(f"    - ANI >= {args.ani_threshold}% (검증 통과): {n_pass}개")
    print(f"    - ANI < {args.ani_threshold}% (재검토 필요): {n_fail_low_ani}개")
    print(f"  - fastANI 계산 실패(너무 먼 genome 등): {len(jobs) - n_tested}개")
    print(f"스킵됨(reference/query 없음 등): {len(skipped)}개")
    print(f"\n상세 결과: {os.path.join(args.outdir, 'reference_ani_validation.tsv')}")

    if n_fail_low_ani > 0:
        print(f"\n[재검토 필요 목록]")
        low = df[(df["status"] == "OK") & (df["ani"] < args.ani_threshold)]
        print(low[["sample_id", "reference_species", "ani"]].to_string(index=False))


if __name__ == "__main__":
    main()
