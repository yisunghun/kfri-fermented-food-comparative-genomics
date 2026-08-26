#!/usr/bin/env bash
###############################################################################
# consolidate_wgs.sh (v2)
#
# Macrogen WGS 결과 폴더(연도별 > Analysis_Data_Done > 시료폴더 > contig단위 어노테이션)를
# 비교유전체 분석(GTDB-Tk, fastANI, Roary/Panaroo 등)에 바로 쓸 수 있는 형태로 병합합니다.
#
# v2 변경점: 종명 판정을 PDF 파일명 대신 *_BLAST.xlsx의 'Result' 시트 기반으로 변경.
#            (extract_species_from_blast.py 호출, 염색체 contig 기준 대표 종명 + 오염 의심 플래그)
#
# 가정하는 구조:
#   <SRC_ROOT>/<연도>/Analysis_Data_Done/<시료폴더>/<AssemblyDir>/consensus.fasta
#   <SRC_ROOT>/<연도>/Analysis_Data_Done/<시료폴더>/<AssemblyDir>/<AssemblyDir>_BLAST.xlsx
#   <SRC_ROOT>/<연도>/Analysis_Data_Done/<시료폴더>/<AssemblyDir>/<AssemblyDir>/contig*/contig*.{faa,ffn,gff}
#
# 필요 파일: extract_species_from_blast.py 를 이 스크립트와 같은 위치(혹은 PATH)에 두세요.
#           pip install pandas openpyxl 필요.
#
# 사용법:
#   chmod +x consolidate_wgs.sh
#   ./consolidate_wgs.sh "/mnt/f/WGS_Results/#Whole_Genome_Sequencing_Macrogen" "/mnt/f/WGS_Consolidated"
#
# 처음엔 연도 하나만 넣어서 테스트 권장:
#   ./consolidate_wgs.sh "/mnt/f/WGS_Results/#Whole_Genome_Sequencing_Macrogen/2021" "/mnt/f/WGS_Consolidated_test"
###############################################################################

set -uo pipefail

SRC_ROOT="${1:?사용법: $0 <소스루트경로> <출력루트경로> [오염판정최소길이bp]}"
OUT_ROOT="${2:?사용법: $0 <소스루트경로> <출력루트경로> [오염판정최소길이bp]}"
CONTAM_MIN_LEN="${3:-20000}"   # 이 길이(bp) 미만 contig는 오염 판정에서 제외 (기본 20000)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY_EXTRACTOR="$SCRIPT_DIR/extract_species_from_blast.py"

if [ ! -f "$PY_EXTRACTOR" ]; then
    echo "오류: extract_species_from_blast.py 를 $SCRIPT_DIR 에서 찾을 수 없습니다." >&2
    exit 1
fi

mkdir -p "$OUT_ROOT/genomes_fna" "$OUT_ROOT/genomes_faa" "$OUT_ROOT/genomes_ffn" "$OUT_ROOT/genomes_gff"
MASTER_TSV="$OUT_ROOT/master_table.tsv"
LOG_FILE="$OUT_ROOT/consolidate.log"
echo -e "sample_id\tshort_id\tyear\tspecies_guess\tn_unique_hits\tcontam_flag\tn_contigs\tconsensus_fasta_src\tblast_xlsx_src\tmerged_faa\tmerged_ffn\tmerged_gff" > "$MASTER_TSV"
> "$LOG_FILE"

# consensus.fasta를 기준으로 시료를 식별 (시료 1개당 정확히 1개씩 있다고 가정)
find "$SRC_ROOT" -type f -iname "consensus.fasta" | while read -r CONSENSUS; do
    ASSEMBLY_DIR=$(dirname "$CONSENSUS")
    SAMPLE_TOP_DIR=$(dirname "$ASSEMBLY_DIR")
    SAMPLE_ID=$(basename "$SAMPLE_TOP_DIR")            # 예: HN00157641_AMT60212 (고유성 보장용 ID)
    SHORT_ID=$(basename "$ASSEMBLY_DIR")                # 예: AMT60212
    YEAR=$(basename "$(dirname "$(dirname "$SAMPLE_TOP_DIR")")")

    # 종명 판정: <ASSEMBLY_DIR>/<SHORT_ID>_BLAST.xlsx 의 Result 시트 기반
    BLAST_XLSX="$ASSEMBLY_DIR/${SHORT_ID}_BLAST.xlsx"
    if [ -f "$BLAST_XLSX" ]; then
        RESULT_LINE=$(python3 "$PY_EXTRACTOR" "$BLAST_XLSX" --format tsv --contam-min-length "$CONTAM_MIN_LEN" 2>>"$LOG_FILE")
        SPECIES_GUESS=$(echo "$RESULT_LINE" | cut -f1)
        N_UNIQUE=$(echo "$RESULT_LINE" | cut -f2)
        CONTAM_FLAG=$(echo "$RESULT_LINE" | cut -f3)
        if [ "$CONTAM_FLAG" = "yes" ]; then
            echo "[WARN] Possible contamination/mixed strain (species disagreement among contigs): $SAMPLE_ID ($SPECIES_GUESS, ${N_UNIQUE} species detected)" >> "$LOG_FILE"
        fi
    else
        SPECIES_GUESS="unknown(no_blast_xlsx)"
        N_UNIQUE=0
        CONTAM_FLAG="unknown"
        echo "[WARN] BLAST xlsx not found: $BLAST_XLSX" >> "$LOG_FILE"
    fi

    # contig* 어노테이션 폴더 탐색 (ASSEMBLY_DIR 하위 최대 2단계까지)
    mapfile -t CONTIG_DIRS < <(find "$ASSEMBLY_DIR" -mindepth 1 -maxdepth 2 -type d -iname "contig*" | sort)

    OUT_FAA="$OUT_ROOT/genomes_faa/${SAMPLE_ID}.faa"
    OUT_FFN="$OUT_ROOT/genomes_ffn/${SAMPLE_ID}.ffn"
    OUT_GFF="$OUT_ROOT/genomes_gff/${SAMPLE_ID}.gff"
    OUT_FNA="$OUT_ROOT/genomes_fna/${SAMPLE_ID}.fna"
    > "$OUT_FAA"; > "$OUT_FFN"; > "$OUT_GFF"

    if [ "${#CONTIG_DIRS[@]}" -gt 0 ]; then
        for CDIR in "${CONTIG_DIRS[@]}"; do
            CNAME=$(basename "$CDIR")
            # 헤더에 시료ID_contig이름을 prefix로 붙여 서로 다른 균주간 유전자 ID 충돌 방지
            [ -f "$CDIR/$CNAME.faa" ] && sed "s/^>/>${SAMPLE_ID}_${CNAME}_/" "$CDIR/$CNAME.faa" >> "$OUT_FAA"
            [ -f "$CDIR/$CNAME.ffn" ] && sed "s/^>/>${SAMPLE_ID}_${CNAME}_/" "$CDIR/$CNAME.ffn" >> "$OUT_FFN"
            [ -f "$CDIR/$CNAME.gff" ] && grep -v "^##FASTA" "$CDIR/$CNAME.gff" >> "$OUT_GFF"
        done
    else
        echo "[WARN] No contig folders found (structure may differ, manual check needed): $ASSEMBLY_DIR" >> "$LOG_FILE"
    fi

    # Roary/Panaroo 호환을 위해 gff 끝에 시퀀스(##FASTA) 첨부.
    # 주의: 최상위 consensus.fasta의 헤더가 각 contig gff의 seqid와 어긋날 수 있어
    #       (서로 다른 소스), 반드시 같은 contig 폴더의 .fna를 그 gff와 짝지어 붙인다.
    #       (같은 어노테이션 실행 결과물이라 seqid 일치가 보장됨)
    if [ -s "$OUT_GFF" ] && [ "${#CONTIG_DIRS[@]}" -gt 0 ]; then
        echo "##FASTA" >> "$OUT_GFF"
        for CDIR in "${CONTIG_DIRS[@]}"; do
            CNAME=$(basename "$CDIR")
            [ -f "$CDIR/$CNAME.fna" ] && cat "$CDIR/$CNAME.fna" >> "$OUT_GFF"
        done
    fi

    cp "$CONSENSUS" "$OUT_FNA"

    N_CONTIGS="${#CONTIG_DIRS[@]}"
    echo -e "${SAMPLE_ID}\t${SHORT_ID}\t${YEAR}\t${SPECIES_GUESS}\t${N_UNIQUE}\t${CONTAM_FLAG}\t${N_CONTIGS}\t${CONSENSUS}\t${BLAST_XLSX}\t${OUT_FAA}\t${OUT_FFN}\t${OUT_GFF}" >> "$MASTER_TSV"
done

echo "완료."
echo "  - 마스터 테이블: $MASTER_TSV"
echo "  - 경고/이슈 로그: $LOG_FILE (오염 의심 시료, 구조 이상 시료 등)"
echo "  - 병합된 genome 파일: $OUT_ROOT/genomes_{fna,faa,ffn,gff}/"
