#!/usr/bin/env python3
"""
extract_species_from_blast.py

시료별 *_BLAST.xlsx 파일의 'Result' 시트에서 genome 대표 종명을 추출합니다.

Result 시트 구조 (확인된 실제 형태):
  - 실제 컬럼 헤더는 엑셀상 2번째 행에 있음 (병합 셀 때문에 pandas 기본 읽기론 밀림)
    → Name, Q_Length, Q_Start, Q_End, Q_Coverage, Description, Accession,
      S_Length, S_Start, S_End, S_Coverage, Bit, E-value,
      I_Match/Total, I_Pct.(%), G_Match/Total, G_Pct.(%)
  - 한 시료 안에서도 contig마다 서로 다른 BLAST hit이 나옴
    (예: contig1=염색체, contig2~4=plasmid)
  - 종명은 Description 컬럼 텍스트 안에 accession과 함께 섞여 있음
    예: "CP031702.1 Lactobacillus plantarum strain IDCC3501 chromosome, complete genome"

판정 로직:
  1) Description에 'plasmid'가 포함된 행은 제외 (염색체 판정용)
  2) 남은 행 중 Q_Length(contig 길이)가 가장 큰 행 = 대표 염색체 hit
  3) 그 Description에서 Genus + species 두 단어를 파싱해 종명으로 사용
  4) 부가로 모든 contig의 hit 종명이 서로 다르면 '오염/혼합 가능성' 경고

사용법:
  # 단일 파일 점검
  python3 extract_species_from_blast.py --inspect "AMT60212_BLAST.xlsx"

  # 단일 파일에서 종명만 추출 (stdout에 종명만 출력, 배치 스크립트에서 캡처하기 좋음)
  python3 extract_species_from_blast.py "AMT60212_BLAST.xlsx"

  # 상세 정보(모든 contig hit, 오염 의심 여부 등)까지 보고 싶을 때
  python3 extract_species_from_blast.py "AMT60212_BLAST.xlsx" --verbose

필요 패키지: pip install pandas openpyxl
"""
import argparse
import re
import sys

import pandas as pd

HEADER_ROW_INDEX = 1  # 엑셀상 실제 컬럼명이 있는 행 (0-indexed)


def load_result_sheet(path: str, sheet_name: str = "Result") -> pd.DataFrame:
    df = pd.read_excel(path, sheet_name=sheet_name, header=HEADER_ROW_INDEX)
    df = df.dropna(subset=["Description"]) if "Description" in df.columns else df
    return df


def parse_species(description: str) -> str:
    """Description 문자열에서 'Genus species' 두 단어를 추출."""
    if not isinstance(description, str):
        return "unknown"
    tokens = description.strip().split()
    # 첫 토큰이 accession 형태(예: CP031702.1)면 제거
    if tokens and re.match(r"^[A-Za-z0-9_]+\.\d+$", tokens[0]):
        tokens = tokens[1:]
    if len(tokens) >= 2 and tokens[0][:1].isupper():
        genus, species = tokens[0], tokens[1].rstrip(",")
        return f"{genus} {species}"
    return " ".join(tokens[:3]) if tokens else "unknown"


def pick_representative(df: pd.DataFrame):
    """plasmid 제외 후 가장 긴 contig(염색체로 추정)의 hit을 대표로 선택."""
    is_plasmid = df["Description"].str.contains("plasmid", case=False, na=False)
    chrom_df = df[~is_plasmid]
    target_df = chrom_df if not chrom_df.empty else df
    if "Q_Length" in target_df.columns:
        target_df = target_df.copy()
        target_df["Q_Length"] = pd.to_numeric(target_df["Q_Length"], errors="coerce")
        top_row = target_df.loc[target_df["Q_Length"].idxmax()]
    else:
        top_row = target_df.iloc[0]
    return top_row


def extract(path: str, sheet_name: str = "Result", verbose: bool = False,
            contam_min_length: int = 20000):
    """
    contam_min_length: 오염/혼합 판정 시 고려할 contig의 최소 길이(bp).
      이 값보다 작은 contig(대개 소형 plasmid)는 종 다양성이 있어도
      정상적인 수평전달 산물일 가능성이 높아 오염 판정에서 제외한다.
    """
    df = load_result_sheet(path, sheet_name)
    if df.empty:
        return "unknown(empty_result_sheet)", 0, []

    top_row = pick_representative(df)
    species = parse_species(top_row["Description"])

    df = df.copy()
    df["_species"] = df["Description"].apply(parse_species)
    df["_qlen"] = pd.to_numeric(df.get("Q_Length"), errors="coerce")

    # 오염 판정은 '염색체급' contig(길이 threshold 이상)만 대상으로 함
    major_df = df[df["_qlen"] >= contam_min_length]
    if major_df.empty:
        major_df = df  # 전부 소형이면 어쩔 수 없이 전체로 판정

    n_unique = major_df["_species"].nunique()
    per_contig = list(zip(df.get("Name", df.index), df["_species"], df["_qlen"]))

    if verbose:
        print(f"[대표 종명] {species}  (기준 contig: {top_row.get('Name', '?')}, "
              f"Q_Length={top_row.get('Q_Length', '?')})", file=sys.stderr)
        print(f"[염색체급 contig(>={contam_min_length}bp) 기준 종 다양성] 총 {n_unique}종 검출", file=sys.stderr)
        for name, sp, qlen in per_contig:
            tag = "" if (pd.notna(qlen) and qlen >= contam_min_length) else "  (소형/판정제외)"
            print(f"    - {name}: {sp}  [Q_Length={qlen}]{tag}", file=sys.stderr)
        if n_unique > 1:
            print("[경고] 염색체급 contig 간에도 종이 불일치 → 실제 오염/혼합 균주 의심, 수동 확인 권장",
                  file=sys.stderr)

    return species, n_unique, per_contig


def inspect(path: str, sheet_name: str = "Result"):
    xls = pd.ExcelFile(path)
    print(f"[시트 목록] {xls.sheet_names}")
    df_raw = pd.read_excel(path, sheet_name=sheet_name, header=HEADER_ROW_INDEX, nrows=5)
    print(f"\n[시트 '{sheet_name}' 컬럼 목록 (header row={HEADER_ROW_INDEX})]")
    for c in df_raw.columns:
        print(f"  - {c}")
    print("\n[상위 5행 미리보기]")
    print(df_raw.to_string())


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("xlsx", help="시료의 *_BLAST.xlsx 경로")
    ap.add_argument("--inspect", action="store_true", help="컬럼 구조만 확인하고 종료")
    ap.add_argument("--sheet-name", default="Result")
    ap.add_argument("--verbose", action="store_true", help="contig별 상세 정보 및 오염 의심 경고 출력(stderr)")
    ap.add_argument("--format", choices=["species", "tsv"], default="species",
                     help="species: 종명만 출력 / tsv: 종명\\tn_unique\\tcontam_flag 출력 (배치 스크립트용)")
    ap.add_argument("--contam-min-length", type=int, default=20000,
                     help="오염 판정 시 고려할 contig 최소 길이(bp), 기본 20000. "
                          "이보다 작은 contig(주로 소형 plasmid)는 종이 달라도 오염 판정에서 제외")
    args = ap.parse_args()

    if args.inspect:
        inspect(args.xlsx, args.sheet_name)
    else:
        try:
            species, n_unique, per_contig = extract(
                args.xlsx, args.sheet_name, args.verbose, args.contam_min_length
            )
        except Exception as e:
            species, n_unique = f"unknown(error:{type(e).__name__})", 0
        if args.format == "tsv":
            contam_flag = "yes" if n_unique > 1 else "no"
            print(f"{species}\t{n_unique}\t{contam_flag}")
        else:
            # 배치 스크립트에서 캡처하기 쉽도록 stdout에는 종명만 출력
            print(species)
