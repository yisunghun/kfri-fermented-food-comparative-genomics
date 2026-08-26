#!/usr/bin/env bash
###############################################################################
# merge_gbk_files.sh
#
# consolidate_wgs.sh와 같은 소스 트리를 사용하여, 각 시료의 contig별 .gbk
# 파일들을 시료 단위 하나의 multi-record GenBank 파일로 병합합니다.
# (antiSMASH 입력용 - GenBank는 여러 LOCUS 레코드를 한 파일에 이어붙이는 것을
#  표준으로 지원하므로 단순 cat으로 병합 가능)
#
# 사용법:
#   chmod +x merge_gbk_files.sh
#   ./merge_gbk_files.sh "/mnt/f/WGS_Results/#Whole_Genome_Sequencing_Macrogen" "/mnt/f/WGS_Consolidated"
###############################################################################

set -uo pipefail

SRC_ROOT="${1:?사용법: $0 <소스루트경로> <출력루트경로(WGS_Consolidated)>}"
OUT_ROOT="${2:?사용법: $0 <소스루트경로> <출력루트경로(WGS_Consolidated)>}"

mkdir -p "$OUT_ROOT/genomes_gbk"
LOG_FILE="$OUT_ROOT/merge_gbk.log"
> "$LOG_FILE"

N_OK=0
N_WARN=0

find "$SRC_ROOT" -type f -iname "consensus.fasta" | while read -r CONSENSUS; do
    ASSEMBLY_DIR=$(dirname "$CONSENSUS")
    SAMPLE_TOP_DIR=$(dirname "$ASSEMBLY_DIR")
    SAMPLE_ID=$(basename "$SAMPLE_TOP_DIR")

    mapfile -t CONTIG_DIRS < <(find "$ASSEMBLY_DIR" -mindepth 1 -maxdepth 2 -type d -iname "contig*" | sort)

    OUT_GBK="$OUT_ROOT/genomes_gbk/${SAMPLE_ID}.gbk"
    > "$OUT_GBK"

    if [ "${#CONTIG_DIRS[@]}" -eq 0 ]; then
        echo "[WARN] contig 폴더 없음: $SAMPLE_ID" >> "$LOG_FILE"
        rm -f "$OUT_GBK"
        continue
    fi

    for CDIR in "${CONTIG_DIRS[@]}"; do
        CNAME=$(basename "$CDIR")
        if [ -f "$CDIR/$CNAME.gbk" ]; then
            cat "$CDIR/$CNAME.gbk" >> "$OUT_GBK"
        else
            echo "[WARN] gbk 없음: $SAMPLE_ID/$CNAME" >> "$LOG_FILE"
        fi
    done

    if [ ! -s "$OUT_GBK" ]; then
        rm -f "$OUT_GBK"
        echo "[WARN] 최종 gbk 비어있음(전부 실패): $SAMPLE_ID" >> "$LOG_FILE"
    fi
done

N_DONE=$(find "$OUT_ROOT/genomes_gbk" -name "*.gbk" | wc -l)
echo "완료. genomes_gbk/ 에 ${N_DONE}개 시료 생성됨."
echo "경고 로그: $LOG_FILE"
