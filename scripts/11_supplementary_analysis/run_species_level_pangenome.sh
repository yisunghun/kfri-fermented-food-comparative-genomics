#!/usr/bin/env bash
###############################################################################
# run_species_level_pangenome.sh
#
# 이미 속(genus) 단위로 정리되어 있는 grouped/by_genus/<Genus>/genomes_gff/
# 폴더에서, 지정한 species_final에 해당하는 isolate들만 골라 별도 폴더로
# 모은 뒤 Panaroo를 실행한다 (genus 대신 species 단위 pangenome).
#
# 사용법:
#   conda activate panaroo_env
#   ./run_species_level_pangenome.sh <master_table_qc.tsv> <grouped_by_genus_root> \
#       <species_pangenome_out_root> <threads> "<Species name 1>" "<Species name 2>" ...
#
#   예)
#   ./run_species_level_pangenome.sh \
#       /mnt/f/WGS_Consolidated/master_table_qc.tsv \
#       /mnt/f/WGS_Consolidated/grouped/by_genus \
#       /mnt/f/WGS_Consolidated/species_pangenome \
#       24 \
#       "Bacillus velezensis" "Lactiplantibacillus plantarum"
###############################################################################
set -euo pipefail

MASTER="${1:?사용법: $0 <master.tsv> <by_genus_root> <out_root> <threads> \"Species 1\" [\"Species 2\" ...]}"
GENUS_ROOT="${2:?}"
OUT_ROOT="${3:?}"
THREADS="${4:-24}"
shift 4
SPECIES_LIST=("$@")

if [ ${#SPECIES_LIST[@]} -eq 0 ]; then
    echo "[오류] 최소 1개 이상의 species 이름을 따옴표로 감싸서 지정하세요." >&2
    exit 1
fi

mkdir -p "$OUT_ROOT"

for SPECIES in "${SPECIES_LIST[@]}"; do
    SAFE_NAME=$(echo "$SPECIES" | tr ' ' '_')
    GENUS=$(echo "$SPECIES" | awk '{print $1}')
    SRC_GFF_DIR="$GENUS_ROOT/$GENUS/genomes_gff"

    if [ ! -d "$SRC_GFF_DIR" ]; then
        echo "[경고] $SRC_GFF_DIR 가 없습니다. '$SPECIES' 건너뜀." >&2
        continue
    fi

    DEST_DIR="$OUT_ROOT/${SAFE_NAME}/genomes_gff"
    mkdir -p "$DEST_DIR"

    echo "=== '$SPECIES' 준비 중 (속 폴더: $SRC_GFF_DIR) ==="

    # 사전 점검: sample_id / species_final 컬럼을 실제로 찾을 수 있는지 먼저 확인
    HEADER_CHECK=$(head -1 "$MASTER" | sed 's/\xef\xbb\xbf//' | tr '\t' '\n' | grep -c -E '^(sample_id|species_final)$')
    if [ "$HEADER_CHECK" -lt 2 ]; then
        echo "  [오류] $MASTER 헤더에서 sample_id 또는 species_final 컬럼을 못 찾았습니다. 건너뜀." >&2
        continue
    fi

    # master_table_qc.tsv에서 이 species_final에 해당하는 sample_id 목록 추출
    SAMPLE_IDS=$(awk -F'\t' -v sp="$SPECIES" '
        NR==1 {
            sub(/^\xef\xbb\xbf/, "", $0)   # UTF-8 BOM 제거 (master_table_qc.tsv 첫 컬럼에 붙어있는 경우 대응)
            for(i=1;i<=NF;i++){
                if($i=="sample_id") sid_col=i
                if($i=="species_final") sp_col=i
            }
            next
        }
        sid_col == 0 || sp_col == 0 { next }  # 컬럼을 못 찾았으면 안전하게 스킵
        $sp_col == sp { print $sid_col }
    ' "$MASTER")

    N_FOUND=0
    N_MISSING=0
    for SID in $SAMPLE_IDS; do
        SRC="$SRC_GFF_DIR/${SID}.gff"
        if [ -f "$SRC" ]; then
            cp "$SRC" "$DEST_DIR/"
            N_FOUND=$((N_FOUND + 1))
        else
            echo "  [경고] $SID 의 gff 파일을 $SRC_GFF_DIR 에서 못 찾음"
            N_MISSING=$((N_MISSING + 1))
        fi
    done
    echo "  -> ${N_FOUND}개 확보 (누락 ${N_MISSING}개)"

    if [ "$N_FOUND" -lt 3 ]; then
        echo "  [경고] 시료가 3개 미만이라 Panaroo 의미가 없습니다. 건너뜀."
        continue
    fi

    PAN_OUT="$OUT_ROOT/${SAFE_NAME}/panaroo_out"
    mkdir -p "$PAN_OUT"
    echo "  Panaroo 실행 중 (${N_FOUND}개 genome)..."
    panaroo -i "$DEST_DIR"/*.gff -o "$PAN_OUT" --clean-mode strict -a core -t "$THREADS"

    echo "  --- '$SPECIES' 결과 요약 ---"
    if [ -f "$PAN_OUT/summary_statistics.txt" ]; then
        cat "$PAN_OUT/summary_statistics.txt"
    fi
    echo ""
done

echo "완료. 각 species의 결과는 $OUT_ROOT/<Species_name>/panaroo_out/ 에 저장되었습니다."
