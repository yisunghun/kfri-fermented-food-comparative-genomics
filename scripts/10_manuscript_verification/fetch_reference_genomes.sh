#!/usr/bin/env bash
###############################################################################
# fetch_reference_genomes.sh
#
# master_table_qc.tsv에서 species_final(genus sp. 등 미해결 제외) 고유 목록을
# 추출하여, 각 종에 대해 NCBI의 대표/reference genome을 1개씩 다운로드한다.
# 이후 각 isolate를 "자신이 판정된 종의 실제 reference genome"과 ANI 비교하여
# 독립적인 taxonomy 검증 근거를 만든다 (기존처럼 220개 시료끼리 비교하는 것과
# 다르게, 외부 기준과 비교하는 것이 핵심).
#
# 사용법:
#   conda activate ncbi_datasets
#   ./fetch_reference_genomes.sh master_table_qc.tsv ./reference_genomes
###############################################################################
set -uo pipefail

MASTER_TSV="${1:?사용법: $0 <master_table_qc.tsv> <출력폴더>}"
OUT_DIR="${2:?사용법: $0 <master_table_qc.tsv> <출력폴더>}"

mkdir -p "$OUT_DIR/zips" "$OUT_DIR/fasta"
LOG="$OUT_DIR/fetch.log"
> "$LOG"

# species_final 컬럼에서 "Genus species" 형태(سp. 등 미해결 제외)만 고유 추출
SPECIES_LIST=$(awk -F'\t' '
NR==1 { for(i=1;i<=NF;i++) if($i=="species_final") col=i; next }
{
    sp = $col
    if (sp ~ /sp\.$/ || sp ~ /^unknown/ || sp == "") next
    print sp
}
' "$MASTER_TSV" | sort -u)

N_TOTAL=$(echo "$SPECIES_LIST" | grep -c .)
echo "고유 종 목록: ${N_TOTAL}개" | tee -a "$LOG"

N_OK=0
N_FAIL=0
: > "$OUT_DIR/species_accession_map.tsv"
echo -e "species\taccession\tstatus" >> "$OUT_DIR/species_accession_map.tsv"

echo "$SPECIES_LIST" | while IFS= read -r SPECIES; do
    [ -z "$SPECIES" ] && continue
    SAFE_NAME=$(echo "$SPECIES" | tr ' ' '_')
    OUT_FASTA="$OUT_DIR/fasta/${SAFE_NAME}.fna"
    if [ -f "$OUT_FASTA" ]; then
        echo "[스킵] $SPECIES (이미 있음)" | tee -a "$LOG"
        continue
    fi

    ZIP_PATH="$OUT_DIR/zips/${SAFE_NAME}.zip"
    echo "[다운로드 시도] $SPECIES" | tee -a "$LOG"
    if datasets download genome taxon "$SPECIES" --reference --assembly-source RefSeq \
        --include genome --filename "$ZIP_PATH" >> "$LOG" 2>&1; then
        mkdir -p "$OUT_DIR/tmp_${SAFE_NAME}"
        unzip -oq "$ZIP_PATH" -d "$OUT_DIR/tmp_${SAFE_NAME}"
        FASTA_SRC=$(find "$OUT_DIR/tmp_${SAFE_NAME}" -name "*.fna" | head -1)
        ACCESSION=$(find "$OUT_DIR/tmp_${SAFE_NAME}" -name "*.fna" | head -1 | xargs -I{} basename {} .fna)
        if [ -n "$FASTA_SRC" ]; then
            cp "$FASTA_SRC" "$OUT_FASTA"
            echo -e "${SPECIES}\t${ACCESSION}\tOK" >> "$OUT_DIR/species_accession_map.tsv"
            echo "  성공: $ACCESSION" | tee -a "$LOG"
        else
            echo -e "${SPECIES}\tNA\tFASTA_NOT_FOUND" >> "$OUT_DIR/species_accession_map.tsv"
            echo "  [경고] fasta 파일을 못 찾음" | tee -a "$LOG"
        fi
        rm -rf "$OUT_DIR/tmp_${SAFE_NAME}" "$ZIP_PATH"
    else
        echo -e "${SPECIES}\tNA\tDOWNLOAD_FAILED" >> "$OUT_DIR/species_accession_map.tsv"
        echo "  [실패] $SPECIES - RefSeq reference 없음, GenBank로 재시도" | tee -a "$LOG"
        # RefSeq에 지정 reference가 없는 종은 GenBank 아무 assembly나 1개 시도
        if datasets download genome taxon "$SPECIES" --assembly-source GenBank --assembly-level complete \
            --include genome --filename "$ZIP_PATH" >> "$LOG" 2>&1; then
            mkdir -p "$OUT_DIR/tmp_${SAFE_NAME}"
            unzip -oq "$ZIP_PATH" -d "$OUT_DIR/tmp_${SAFE_NAME}"
            FASTA_SRC=$(find "$OUT_DIR/tmp_${SAFE_NAME}" -name "*.fna" | head -1)
            if [ -n "$FASTA_SRC" ]; then
                cp "$FASTA_SRC" "$OUT_FASTA"
                echo "  GenBank로 성공" | tee -a "$LOG"
            fi
            rm -rf "$OUT_DIR/tmp_${SAFE_NAME}" "$ZIP_PATH"
        fi
    fi
done

N_FETCHED=$(find "$OUT_DIR/fasta" -name "*.fna" | wc -l)
echo "완료. ${N_FETCHED}/${N_TOTAL}개 종의 reference genome 확보." | tee -a "$LOG"
echo "매핑 표: $OUT_DIR/species_accession_map.tsv"
echo "실패한 종은 로그에서 확인 후, 필요시 NCBI 웹에서 수동으로 확인해주세요."
