#!/usr/bin/env bash
###############################################################################
# run_antismash_batch.sh
#
# genomes_gbk/*.gbk 각각에 대해 antiSMASH를 실행합니다.
# master_table_qc.tsv의 functional_group으로 필터링 가능 (기본: LAB,Bacillus_group).
#
# 사용법:
#   chmod +x run_antismash_batch.sh
#   ./run_antismash_batch.sh <gbk루트> <출력루트> <master_table_qc.tsv> [threads] [그룹,콤마구분]
#
# 예:
#   ./run_antismash_batch.sh "/mnt/f/WGS_Consolidated/genomes_gbk" \
#       "/mnt/f/WGS_Consolidated/antismash_out" \
#       "/mnt/f/WGS_Consolidated/master_table_qc.tsv" \
#       8 "LAB,Bacillus_group"
#
# 사전 조건: conda activate antismash_env
###############################################################################

set -uo pipefail

GBK_ROOT="${1:?사용법: $0 <gbk루트> <출력루트> <master_table_qc.tsv> [threads] [그룹리스트]}"
OUT_ROOT="${2:?사용법: $0 <gbk루트> <출력루트> <master_table_qc.tsv> [threads] [그룹리스트]}"
MASTER_TSV="${3:?사용법: $0 <gbk루트> <출력루트> <master_table_qc.tsv> [threads] [그룹리스트]}"
THREADS="${4:-4}"
GROUP_FILTER="${5:-LAB,Bacillus_group}"

mkdir -p "$OUT_ROOT"
LOG="$OUT_ROOT/antismash_run.log"
> "$LOG"

# functional_group을 genus_final로부터 즉석 계산하여, 지정된 그룹에 속하는 sample_id 목록 추출
# (pandas 등 외부 의존성 없이 awk만으로 처리 - antismash_env에 pandas가 없을 수 있음)
SAMPLE_LIST=$(awk -F'\t' -v groups="$GROUP_FILTER" '
BEGIN {
    n = split(groups, wanted_arr, ",")
    for (i = 1; i <= n; i++) wanted[wanted_arr[i]] = 1

    split("Lactiplantibacillus,Levilactobacillus,Latilactobacillus,Lactobacillus,Lactococcus,Leuconostoc,Weissella,Pediococcus,Enterococcus,Tetragenococcus,Lacticaseibacillus,Limosilactobacillus,Lentilactobacillus,Loigolactobacillus,Fructilactobacillus", lab_arr, ",")
    for (i in lab_arr) LAB[lab_arr[i]] = 1

    split("Bacillus,Paenibacillus,Oceanobacillus,Virgibacillus,Priestia,Rossellomorea,Shouchella,Halobacillus", bac_arr, ",")
    for (i in bac_arr) BAC[bac_arr[i]] = 1
}
NR==1 {
    # sample_id는 항상 1번째 컬럼(BOM 때문에 헤더 텍스트 매칭이 깨질 수 있어 위치로 고정)
    sid_col = 1
    for (i=1; i<=NF; i++) { if ($i=="genus_final") genus_col=i }
    if (genus_col == 0) {
        print "ERROR: genus_final column not found in header" > "/dev/stderr"
        exit 1
    }
    next
}
{
    g = $genus_col
    if (g == "unresolved") fg = "Unresolved"
    else if (g in LAB) fg = "LAB"
    else if (g in BAC) fg = "Bacillus_group"
    else fg = "Other_Environmental"

    if (fg in wanted) print $sid_col
}
' "$MASTER_TSV")

N_TOTAL=$(echo "$SAMPLE_LIST" | grep -c .)
echo "대상 그룹: $GROUP_FILTER -> $N_TOTAL개 시료" | tee -a "$LOG"

N_DONE=0
N_SKIP=0
N_RUN=0
printf '%s\n' "$SAMPLE_LIST" | while IFS= read -r SID; do
    [ -z "$SID" ] && continue
    GBK="$GBK_ROOT/${SID}.gbk"
    if [ ! -f "$GBK" ]; then
        echo "[없음] $SID (gbk 파일 없음)" | tee -a "$LOG"
        continue
    fi

    OUT_DIR="$OUT_ROOT/$SID"
    if [ -d "$OUT_DIR" ] && [ -f "$OUT_DIR/index.html" ]; then
        N_SKIP=$((N_SKIP+1))
        continue
    fi

    echo "[실행] $SID" | tee -a "$LOG"
    START=$(date +%s)
    if antismash --cpus "$THREADS" --output-dir "$OUT_DIR" --genefinding-tool none "$GBK" >> "$LOG" 2>&1; then
        END=$(date +%s)
        N_RUN=$((N_RUN+1))
        echo "  완료 ($((END-START))초)" | tee -a "$LOG"
    else
        echo "  [오류] 실패 - 로그 확인" | tee -a "$LOG"
    fi
done

echo "전체 완료. 새로 실행: $N_RUN, 기존완료 스킵: $N_SKIP / 총 $N_TOTAL" | tee -a "$LOG"
