#!/usr/bin/env bash
###############################################################################
# characterize_kpneumoniae_outlier.sh
#
# 3.8절에서 발견한 CARD=26, VFDB=62의 극단적 이상치 시료가 실제로
# master_table_qc.tsv 상에서 species_final="Klebsiella pneumoniae"로 판정된
# 그 시료인지 자동으로 찾아, 아래 두 가지를 추가로 확인한다:
#   1) MLST (다좌위서열형별, torstenseemann/mlst 도구, PubMLST 스킴 이용)
#   2) PlasmidFinder DB로 plasmid replicon 스크리닝 (abricate)
#      -> 저항성/독성 유전자가 plasmid에 실려 있을 가능성(수평전달 위험) 확인
#
# 사전 설치 (최초 1회):
#   mamba create -n mlst_env -c bioconda -c conda-forge mlst -y
#
# 사용법:
#   conda activate abricate_env   # plasmidfinder DB 스크리닝용
#   ./characterize_kpneumoniae_outlier.sh \
#       master_table_qc.tsv genomes_fna ./kpneumoniae_result
###############################################################################
set -euo pipefail

MASTER="${1:?사용법: $0 <master_table_qc.tsv> <genomes_fna_dir> <outdir>}"
GENOMES_DIR="${2:?}"
OUTDIR="${3:?}"
mkdir -p "$OUTDIR"

# BOM 제거 후 species_final == "Klebsiella pneumoniae"인 sample_id 탐색
SAMPLE_ID=$(awk -F'\t' '
    NR==1 {
        sub(/^\xef\xbb\xbf/, "", $0)
        for(i=1;i<=NF;i++){
            if($i=="sample_id") sid_col=i
            if($i=="species_final") sp_col=i
        }
        next
    }
    $sp_col == "Klebsiella pneumoniae" { print $sid_col }
' "$MASTER")

N_FOUND=$(echo "$SAMPLE_ID" | grep -c . || true)
if [ "$N_FOUND" -ne 1 ]; then
    echo "[오류] Klebsiella pneumoniae로 판정된 시료가 ${N_FOUND}개 발견됨 (1개여야 정상)." >&2
    echo "발견된 시료: $SAMPLE_ID" >&2
    exit 1
fi

echo "=== 대상 시료: $SAMPLE_ID ==="
FASTA="$GENOMES_DIR/${SAMPLE_ID}.fna"
if [ ! -f "$FASTA" ]; then
    echo "[오류] $FASTA 를 찾을 수 없습니다." >&2
    exit 1
fi

echo ""
echo "=== 1) MLST 실행 ==="
if command -v mlst >/dev/null 2>&1; then
    mlst "$FASTA" | tee "$OUTDIR/${SAMPLE_ID}_mlst.tsv"
else
    echo "[경고] mlst 명령을 찾을 수 없습니다. 'conda activate mlst_env' 후 다시 실행하세요." >&2
fi

echo ""
echo "=== 2) PlasmidFinder 스크리닝 ==="
if command -v abricate >/dev/null 2>&1; then
    abricate --db plasmidfinder "$FASTA" | tee "$OUTDIR/${SAMPLE_ID}_plasmidfinder.tsv"
else
    echo "[경고] abricate 명령을 찾을 수 없습니다. 'conda activate abricate_env' 후 다시 실행하세요." >&2
fi

echo ""
echo "완료. 결과 저장 위치: $OUTDIR/"
