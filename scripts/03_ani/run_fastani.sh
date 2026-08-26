#!/usr/bin/env bash
###############################################################################
# run_fastani.sh
#
# genomes_fna/*.fna 전체에 대해 all-vs-all ANI(Average Nucleotide Identity)를
# fastANI로 계산합니다.
#
# 사용법:
#   chmod +x run_fastani.sh
#   ./run_fastani.sh "/mnt/f/WGS_Consolidated/genomes_fna" "/mnt/f/WGS_Consolidated/ani_out" [threads]
#
# 사전 조건: conda activate compgenomics (fastANI 설치된 환경)
###############################################################################

set -euo pipefail

FNA_DIR="${1:?사용법: $0 <genomes_fna 경로> <출력경로> [threads]}"
OUT_DIR="${2:?사용법: $0 <genomes_fna 경로> <출력경로> [threads]}"
THREADS="${3:-4}"

mkdir -p "$OUT_DIR"
GENOME_LIST="$OUT_DIR/genome_list.txt"
ANI_OUT="$OUT_DIR/ani_result.tsv"

find "$FNA_DIR" -maxdepth 1 -type f -iname "*.fna" | sort > "$GENOME_LIST"

N_GENOMES=$(wc -l < "$GENOME_LIST")
echo "총 ${N_GENOMES}개 genome 발견. all-vs-all ANI 계산 시작..."
echo "(참고: ${N_GENOMES}x${N_GENOMES} 조합이라 시료 수가 많으면 시간이 꽤 걸립니다)"

fastANI \
    --ql "$GENOME_LIST" \
    --rl "$GENOME_LIST" \
    -o "$ANI_OUT" \
    -t "$THREADS"

echo "완료: $ANI_OUT"
echo "  컬럼: query_genome, reference_genome, ANI(%), count_bidirectional_fragments, total_query_fragments"
echo "  (참고: ANI가 대략 80% 미만인 쌍은 fastANI가 아예 결과를 출력하지 않습니다 - 정상 동작입니다)"
