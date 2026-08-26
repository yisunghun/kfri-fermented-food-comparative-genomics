#!/usr/bin/env bash
###############################################################################
# run_panaroo_all_groups.sh
#
# organize_by_group.py가 만든 grouped/by_genus, grouped/by_functional_group
# 하위의 각 그룹 폴더(genomes_gff/*.gff)에 대해 Panaroo pangenome을 실행합니다.
#
# 기본적으로 서로 다른 속이 섞인 그룹(_minor_genera, unresolved, Unresolved,
# Other_Environmental)은 제외합니다 (pangenome이 생물학적으로 의미 없음).
#
# 사용법:
#   chmod +x run_panaroo_all_groups.sh
#   ./run_panaroo_all_groups.sh <grouped 루트> <출력루트> [threads] [min_genomes] [추가제외그룹,콤마구분]
#
# 예:
#   ./run_panaroo_all_groups.sh "/mnt/f/WGS_Consolidated/grouped" "/mnt/f/WGS_Consolidated/pangenome" 8 3
#
# 사전 조건: conda activate compgenomics (panaroo 설치된 환경)
###############################################################################

set -uo pipefail

GROUPED_ROOT="${1:?사용법: $0 <grouped루트> <출력루트> [threads] [min_genomes] [추가제외그룹]}"
OUT_ROOT="${2:?사용법: $0 <grouped루트> <출력루트> [threads] [min_genomes] [추가제외그룹]}"
THREADS="${3:-4}"
MIN_GENOMES="${4:-3}"
EXTRA_EXCLUDE="${5:-}"

DEFAULT_EXCLUDE="_minor_genera unresolved Unresolved Other_Environmental"
EXCLUDE_LIST="$DEFAULT_EXCLUDE $(echo "$EXTRA_EXCLUDE" | tr ',' ' ')"

mkdir -p "$OUT_ROOT"
LOG="$OUT_ROOT/panaroo_run.log"
> "$LOG"

is_excluded() {
    local name="$1"
    for ex in $EXCLUDE_LIST; do
        [ "$name" = "$ex" ] && return 0
    done
    return 1
}

echo "제외 그룹: $EXCLUDE_LIST" | tee -a "$LOG"

find "$GROUPED_ROOT" -type d -name "genomes_gff" | sort | while read -r GFF_DIR; do
    GROUP_DIR=$(dirname "$GFF_DIR")
    GROUP_NAME=$(basename "$GROUP_DIR")
    CATEGORY=$(basename "$(dirname "$GROUP_DIR")")   # by_genus 또는 by_functional_group

    if is_excluded "$GROUP_NAME"; then
        echo "[제외] $CATEGORY/$GROUP_NAME" | tee -a "$LOG"
        continue
    fi

    N_GFF=$(find "$GFF_DIR" -maxdepth 1 -iname "*.gff" | wc -l)
    if [ "$N_GFF" -lt "$MIN_GENOMES" ]; then
        echo "[SKIP] $CATEGORY/$GROUP_NAME (n=$N_GFF < $MIN_GENOMES)" | tee -a "$LOG"
        continue
    fi

    OUT_DIR="$OUT_ROOT/$CATEGORY/$GROUP_NAME"
    if [ -f "$OUT_DIR/gene_presence_absence.csv" ]; then
        echo "[이미완료] $CATEGORY/$GROUP_NAME -> 스킵 (재실행하려면 $OUT_DIR 삭제 후 재시도)" | tee -a "$LOG"
        continue
    fi
    mkdir -p "$OUT_DIR"

    echo "[실행] $CATEGORY/$GROUP_NAME (n=$N_GFF genomes) -> $OUT_DIR" | tee -a "$LOG"
    START_TS=$(date +%s)
    if panaroo -i "$GFF_DIR"/*.gff -o "$OUT_DIR" --clean-mode strict -a core -t "$THREADS" >> "$LOG" 2>&1; then
        END_TS=$(date +%s)
        echo "  완료 ($((END_TS - START_TS))초 소요)" | tee -a "$LOG"
    else
        echo "  [오류] 실패 - 로그($LOG) 확인 필요" | tee -a "$LOG"
    fi
done

echo "전체 완료. 상세 로그: $LOG"
echo "각 그룹 결과: $OUT_ROOT/<category>/<group>/ 안의 summary_statistics.txt, gene_presence_absence.csv 등"
