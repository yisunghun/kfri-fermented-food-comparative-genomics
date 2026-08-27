#!/usr/bin/env bash
###############################################################################
# rerun_lactococcus_enterococcus_pangenome.sh (v2 - 실제 폴더 구조 반영)
#
# 독립 reference ANI 검증(3.3절)으로 아래 2개 isolate가
#   Lactococcus lactis -> Enterococcus lactis
# 로 재동정되었음이 확인됨:
#   - HN00180670_CHKJ1127
#   - HN00196390_CHKJ1122
#
# 실제 폴더 구조: <grouped_root>/<Genus>/genomes_gff/<sample_id>.gff
#                  <grouped_root>/<Genus>/genomes_faa/<sample_id>.faa
#                  <grouped_root>/<Genus>/group_manifest.tsv
#
# 이 2개 isolate의 gff/faa 파일을 Lactococcus -> Enterococcus로 옮기고,
# group_manifest.tsv도 함께 갱신한 뒤, 두 속에 대해서만 Panaroo를 재실행하여
# Table 5(genus-level pangenome summary)를 보정된 구성(n=7/n=23) 기준으로
# 갱신한다.
#
# 사용법:
#   conda activate panaroo_env
#   ./rerun_lactococcus_enterococcus_pangenome.sh \
#       /mnt/f/WGS_Consolidated/grouped/by_genus \
#       /mnt/f/WGS_Consolidated/pangenome \
#       24
###############################################################################
set -euo pipefail

GROUPED_ROOT="${1:?사용법: $0 <grouped_root> <pangenome_out_root> <threads>}"
PAN_OUT_ROOT="${2:?사용법: $0 <grouped_root> <pangenome_out_root> <threads>}"
THREADS="${3:-24}"

REASSIGNED_SAMPLES=("HN00180670_CHKJ1127" "HN00196390_CHKJ1122")

LACTO_DIR="$GROUPED_ROOT/Lactococcus"
ENTERO_DIR="$GROUPED_ROOT/Enterococcus"

for d in "$LACTO_DIR/genomes_gff" "$ENTERO_DIR/genomes_gff"; do
    if [ ! -d "$d" ]; then
        echo "[오류] $d 폴더를 찾을 수 없습니다." >&2
        exit 1
    fi
done

echo "=== 1) 재동정된 isolate GFF/FAA 파일 이동 ==="
for SID in "${REASSIGNED_SAMPLES[@]}"; do
    for EXT in gff faa; do
        SRC="$LACTO_DIR/genomes_${EXT}/${SID}.${EXT}"
        DST_DIR="$ENTERO_DIR/genomes_${EXT}"
        DST="$DST_DIR/${SID}.${EXT}"
        if [ -f "$SRC" ]; then
            mkdir -p "$DST_DIR"
            mv "$SRC" "$DST"
            echo "이동: $SRC -> $DST"
        elif [ -f "$DST" ]; then
            echo "[스킵] $SID.$EXT 는 이미 Enterococcus 폴더에 있음"
        else
            echo "[경고] $SID.$EXT 를 못 찾음 (Lactococcus/Enterococcus 어느 쪽에도 없음)" >&2
        fi
    done
done

echo ""
echo "=== 2) group_manifest.tsv 갱신 ==="
for SID in "${REASSIGNED_SAMPLES[@]}"; do
    LACTO_MANIFEST="$LACTO_DIR/group_manifest.tsv"
    ENTERO_MANIFEST="$ENTERO_DIR/group_manifest.tsv"

    if [ -f "$LACTO_MANIFEST" ] && grep -q "^${SID}"$'\t' "$LACTO_MANIFEST"; then
        ROW=$(grep "^${SID}"$'\t' "$LACTO_MANIFEST")
        NEW_ROW=$(echo "$ROW" | awk -F'\t' -v OFS='\t' '{ $2="Enterococcus"; $3="Enterococcus"; $4="Enterococcus lactis"; print }')
        grep -v "^${SID}"$'\t' "$LACTO_MANIFEST" > "${LACTO_MANIFEST}.tmp" && mv "${LACTO_MANIFEST}.tmp" "$LACTO_MANIFEST"
        echo "$NEW_ROW" >> "$ENTERO_MANIFEST"
        echo "manifest 갱신: $SID -> Enterococcus 로 이동"
    else
        echo "[참고] $SID 는 이미 Lactococcus manifest에 없음 (이미 처리됐거나 이름 확인 필요)"
    fi
done

N_LACTO=$(find "$LACTO_DIR/genomes_gff" -name "*.gff" | wc -l)
N_ENTERO=$(find "$ENTERO_DIR/genomes_gff" -name "*.gff" | wc -l)
echo ""
echo "이동 후 Lactococcus: ${N_LACTO}개 (기대값 7), Enterococcus: ${N_ENTERO}개 (기대값 23)"

echo ""
echo "=== 3) Panaroo 재실행 ==="
for GENUS in Lactococcus Enterococcus; do
    IN_DIR="$GROUPED_ROOT/$GENUS/genomes_gff"
    OUT_DIR="$PAN_OUT_ROOT/${GENUS}_corrected"
    mkdir -p "$OUT_DIR"
    N_GFF=$(find "$IN_DIR" -name "*.gff" | wc -l)
    echo "--- $GENUS (${N_GFF}개 genome) ---"
    panaroo -i "$IN_DIR"/*.gff -o "$OUT_DIR" --clean-mode strict -a core -t "$THREADS"
done

echo ""
echo "=== 4) 결과 요약 ==="
for GENUS in Lactococcus Enterococcus; do
    OUT_DIR="$PAN_OUT_ROOT/${GENUS}_corrected"
    SUMMARY="$OUT_DIR/summary_statistics.txt"
    if [ -f "$SUMMARY" ]; then
        echo "--- $GENUS ---"
        cat "$SUMMARY"
        echo ""
    fi
done

echo "완료. 위 요약 내용을 그대로 복사해서 알려주시면 Table 5를 갱신하겠습니다."
