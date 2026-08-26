#!/usr/bin/env bash
###############################################################################
# run_abricate_batch.sh
#
# genomes_fna/*.fna 전체(또는 필터링된 시료)에 대해 abricate로
# CARD(항생제내성유전자), VFDB(병원성인자) 스크리닝을 실행하고,
# 각 DB별로 abricate --summary 매트릭스를 생성합니다.
#
# 사용법:
#   chmod +x run_abricate_batch.sh
#   ./run_abricate_batch.sh <genomes_fna 경로> <출력루트>
#
# 사전 조건: conda activate abricate_env (abricate --setupdb 완료된 상태)
###############################################################################

set -uo pipefail

FNA_DIR="${1:?사용법: $0 <genomes_fna경로> <출력루트>}"
OUT_ROOT="${2:?사용법: $0 <genomes_fna경로> <출력루트>}"
DBS="card vfdb"

for DB in $DBS; do
    mkdir -p "$OUT_ROOT/$DB"
done
LOG="$OUT_ROOT/abricate_run.log"
> "$LOG"

N_FNA=$(find "$FNA_DIR" -maxdepth 1 -iname "*.fna" | wc -l)
echo "총 ${N_FNA}개 genome 발견" | tee -a "$LOG"

for DB in $DBS; do
    echo "===== DB: $DB =====" | tee -a "$LOG"
    N=0
    find "$FNA_DIR" -maxdepth 1 -iname "*.fna" | sort | while read -r FNA; do
        SID=$(basename "$FNA" .fna)
        OUT_TAB="$OUT_ROOT/$DB/${SID}.tab"
        if [ -f "$OUT_TAB" ]; then
            continue
        fi
        abricate --db "$DB" --quiet "$FNA" > "$OUT_TAB" 2>> "$LOG"
        N=$((N+1))
        if [ $((N % 20)) -eq 0 ]; then
            echo "  [$DB] ${N}개 처리..." | tee -a "$LOG"
        fi
    done

    # 모든 개별 결과를 하나의 요약 매트릭스로 (시료 x 유전자, 존재시 %coverage/identity 등 표시)
    SUMMARY_OUT="$OUT_ROOT/${DB}_summary.tsv"
    abricate --summary "$OUT_ROOT/$DB"/*.tab > "$SUMMARY_OUT" 2>> "$LOG"
    echo "  요약 매트릭스 저장: $SUMMARY_OUT" | tee -a "$LOG"
done

echo "전체 완료." | tee -a "$LOG"
