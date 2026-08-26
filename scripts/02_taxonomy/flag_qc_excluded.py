#!/usr/bin/env python3
"""
flag_qc_excluded.py

diagnose_pangenome_genomes.py에서 발견된, Panaroo 유전자 추출률이 비정상적으로
낮았던 시료들을 master_table_final.tsv에 'qc_status' 컬럼으로 기록합니다.
(원인: 조사 결과 원본 gff/fasta 자체는 정상 - seqid 일치, 서열 길이 정상,
 ID 중복 없음 - 이지만 Panaroo 처리 시 유전자 대부분이 유실됨. 정확한 원인은
 미상이나 전체 220개 중 5개(2.3%)로 비중이 작아 이번 분석에서는 제외하고 진행.)

사용법:
  python3 flag_qc_excluded.py \
      --master master_table_final.tsv \
      --out master_table_qc.tsv
"""
import argparse

import pandas as pd

QC_EXCLUDED_SAMPLES = {
    "HN00222446_Sea08-36": "Panaroo gene recovery rate 0.02% (1 of 5276 genes recognized; cause unknown)",
    "HN00251139_BGIL25041": "Panaroo gene recovery rate 0% (0 of 3952 genes recognized; cause unknown; prior contamination flag)",
    "HN00251139_BGIL25015": "Panaroo gene recovery rate 1% (44 of 4240 genes recognized; cause unknown)",
    "HN00251139_BJS25031": "Panaroo gene recovery rate 23% (1186 of 5076 genes recognized; cause unknown; prior contamination flag)",
    "HN00280011_KCUT25001": "Panaroo gene recovery rate 45% (1744 of 3874 genes recognized; cause unknown)",
}

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--master", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    df = pd.read_csv(args.master, sep="\t")
    df["qc_status"] = df["sample_id"].apply(
        lambda s: "excluded" if s in QC_EXCLUDED_SAMPLES else "ok"
    )
    df["qc_note"] = df["sample_id"].apply(lambda s: QC_EXCLUDED_SAMPLES.get(s, ""))

    df.to_csv(args.out, sep="\t", index=False, encoding="utf-8-sig")
    print(f"저장 완료: {args.out}")
    print(f"QC 제외 시료: {(df['qc_status'] == 'excluded').sum()}건 / 총 {len(df)}건")
