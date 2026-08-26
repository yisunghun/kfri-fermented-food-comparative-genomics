# 발효식품 유래 세균 220종 비교유전체 분석 파이프라인

한국 발효식품(김치, 장류 등)에서 유래한 약 220개 세균 전장유전체를 대상으로,
Macrogen 원본 산출물 → 종 판정 → ANI 클러스터링 → Pangenome 분석 → 기능(COG)
프로파일 비교 → 항생제내성/병원성 스크리닝 → 이차대사산물(BGC) 스크리닝 →
최종 통합 계통수+히트맵까지 이어지는 재현 가능한 프로토콜입니다.

이 문서는 각 단계의 **최종적으로 성공한 버전**만 기술합니다. 중간에 실패했던
시도가 있었던 부분은 "⚠ 겪었던 문제" 박스로 원인과 해결법을 짧게 남겨서,
동일한 프로젝트를 재현하는 사람이 같은 실수를 반복하지 않도록 했습니다.

---

## 0. 원본 데이터 구조

Macrogen이 시퀀싱/어노테이션 결과를 넘겨준 구조는 연도별 폴더로 되어 있습니다:

```
<SRC_ROOT>/<연도>/Analysis_Data_Done/<시료폴더>/<AssemblyDir>/
    consensus.fasta                     # 전체 genome 조립본(모든 contig)
    <short_id>_BLAST.xlsx               # contig별 BLAST hit ('Result' 시트)
    <AssemblyDir>/contig1/ contig2/ …   # contig 단위 개별 어노테이션(Prokka/tbl2asn류)
        contigN.gff  contigN.faa  contigN.ffn  contigN.fna  contigN.gbk

<시료폴더>/<short_id>_FunctionalAnnotation/FunctionalAnnotation/
    annotation_EggNOG.xlsx              # 'Eggnog_Count'(또는 'Eggnog Count') 시트
```

**파이프라인 전반에서 활용한 핵심 구조적 사실:** 각 contig의 `.gff`,
`.faa`, `.ffn`, `.fna`, `.gbk`는 **같은 어노테이션 실행에서 함께 생성**됐기
때문에, 특정 contig의 `.fna` 헤더는 그 contig 자신의 `.gff` 안 `seqid`와
**반드시 일치**합니다. 반면 최상위 `consensus.fasta`는 어셈블러가 별도로
만든 것이라 헤더가 일치한다는 보장이 없습니다. 이 차이가 프로젝트 중 가장
큰 버그의 원인이었고(§3.1 참고), 이번 프로젝트에서 얻은 가장 중요한 교훈입니다.

---

## 1. 환경설정

서로 다른 Python 버전/의존성 요구사항이 충돌해서 **conda/mamba 환경을 4개**
분리해서 사용했습니다.

```bash
# 1) 메인 분석 환경 (fastANI, pandas/numpy/scipy/matplotlib,
#    scikit-learn, FastTree, Biopython)
mamba create -n compgenomics --override-channels -c bioconda -c conda-forge \
    fastani pandas scipy matplotlib numpy openpyxl -y
conda activate compgenomics
pip install scikit-learn --break-system-packages
mamba install -c bioconda -c conda-forge fasttree biopython -y

# 2) Panaroo(pangenome) - 더 낮은 Python 버전 필요; compgenomics가
#    자동으로 잡는 Python 3.14와 충돌(python-edlib이 아직 3.14 미지원)
mamba create -n panaroo_env --override-channels -c bioconda -c conda-forge \
    python=3.10 panaroo -y

# 3) antiSMASH (이차대사산물 BGC 탐지)
mamba create -n antismash_env --override-channels -c bioconda -c conda-forge \
    antismash -y
conda activate antismash_env
download-antismash-databases
antismash --check-prereqs

# 4) abricate (CARD / VFDB 스크리닝)
mamba create -n abricate_env --override-channels -c bioconda -c conda-forge \
    abricate -y
conda activate abricate_env
abricate --setupdb
```

> ⚠ **겪었던 문제 — `mamba create`가 느리거나 불안정함**
> `--override-channels` 없이 실행하면 Anaconda의 기본 채널(`pkgs/main`,
> `pkgs/r`, 상용 채널)까지 같이 인덱싱해서 수 분이 걸리고, 가끔
> "environment specs not solvable" 오류로 아예 실패합니다. 항상
> `--override-channels`를 `-c bioconda -c conda-forge`와 함께 명시하세요.

matplotlib에서 한글 폰트가 필요한 경우(이 파이프라인의 최종 산출물은 영문이지만,
중간 점검용으로 한글 그래프를 그릴 때 필요):

```bash
sudo apt-get install -y fonts-nanum
rm -rf ~/.cache/matplotlib   # matplotlib이 새 폰트를 다시 스캔하도록 강제
```

---

## 2. 1단계 — contig 단위 어노테이션을 시료 단위로 병합

**스크립트:** `consolidate_wgs.sh`

시료별로 흩어진 `contigN.faa/.ffn/.gff`를 하나로 합치고, 조립본을
`genomes_fna/`로 모으고, BLAST 기반 종 판정과 오염 의심 플래그를 담은
`master_table.tsv`를 만듭니다.

```bash
./consolidate_wgs.sh \
    "/mnt/f/WGS_Results/#Whole_Genome_Sequencing_Macrogen" \
    "/mnt/f/WGS_Consolidated"
```

최종 버전의 핵심 설계 포인트:

- **종 판정은 PDF 보고서 파일명이 아니라 `*_BLAST.xlsx`에서 가져옵니다.**
  PDF는 편집된 요약본이라 실제 최상위 BLAST hit과 다를 수 있습니다.
- `Result` 시트의 실제 헤더 행은 **인덱스 1**(0-indexed)입니다. 엑셀의
  병합 셀 때문에 실제 컬럼명(`Name, Q_Length, …, Description, …`)이
  한 행 아래로 밀려있습니다.
- 종명은 **plasmid가 아닌 가장 긴 contig(=염색체)**의 `Description`
  필드에서 가져옵니다 — 단순히 첫 BLAST hit을 쓰지 않는 이유는, 작은
  plasmid contig는 수평전달이 흔해서 오염이 아닌데도 다른 종/속을
  가리킬 수 있기 때문입니다.
- **20kb 이상**(`--contam-min-length`, 조정 가능)인 contig들끼리만 종이
  갈릴 때만 `contam_flag=yes`로 표시합니다 — 작은 plasmid는 오염 판정에서
  제외합니다.
- 각 시료 `.gff`에 붙이는 `##FASTA` 블록은 **각 contig 자신의 `.fna`**를
  이어붙여 만듭니다. 최상위 `consensus.fasta`는 절대 쓰지 않습니다(아래
  문제 박스 참고).

```bash
# 종 추출은 다음 스크립트에 위임:
python3 extract_species_from_blast.py <시료>_BLAST.xlsx --verbose
```

> ⚠ **겪었던 문제 — GFF/FASTA 헤더 불일치로 Panaroo(및 antiSMASH도)가
> 깨질 뻔함**
> 원래 코드는 `##FASTA` 뒤에 최상위 `consensus.fasta`를 통째로 붙였습니다.
> 그 헤더(예: `>contig1`)가 각 contig 자신의 `.gff` 안 `seqid`(예:
> `gnl|MG|AMT60212_1_1`)와 **일치하지 않았습니다** — 어셈블러와
> 어노테이터가 서로 다른 명명 체계를 쓴 겁니다. 그 결과 Panaroo가 모든
> 유전자에 대해 `Invalid gene sequence!` 오류를 냈습니다. **해결:** 각
> contig 자신의 `.fna`를 그 contig의 `.gff`와 1:1로 짝지어 붙임으로써,
> 명명 체계가 어떻든 `seqid` 일치를 보장했습니다.

> ⚠ **겪었던 문제 — 오염 플래그가 너무 민감함**
> 처음엔 작은 plasmid까지 포함해서 모든 contig를 비교하다 보니 시료의
> 약 75%가 "오염"으로 잡혔습니다. 20kb 이상 contig로만 비교 범위를
> 제한하니 현실적인 수준으로 줄었습니다.

---

## 3. 2단계 — 분류학 정규화

**스크립트:** `normalize_taxonomy.py`

2020년 이전 `Lactobacillus` 속명을 재분류(Zheng et al. 2020) 이후의 신규
속명으로 매핑합니다 (예: `Lactobacillus plantarum` →
`Lactiplantibacillus plantarum`). BLAST 참조 DB가 신구 명명법을 혼용하고
있어서, 유전체적으로 완전히 동일한 hit인데도 이후 단계에서 가짜
ANI/분류 불일치로 보일 수 있기 때문에 필요합니다.

```bash
python3 normalize_taxonomy.py \
    --master master_table.tsv --out master_table_normalized.tsv
```

---

## 4. 3단계 — 평균염기서열동일성(ANI)

**스크립트:** `run_fastani.sh`, `analyze_ani.py`

```bash
./run_fastani.sh genomes_fna/ ani_out/ 8

python3 analyze_ani.py \
    --ani ani_out/ani_result.tsv \
    --master master_table_normalized.tsv \
    --outdir ani_analysis/ \
    --species-col species_normalized
```

산출물: 대칭 ANI 매트릭스(CSV), 속(genus) 색상으로 구분된 덴드로그램(PDF),
그리고 BLAST 기반 종 판정과 95%-ANI 종 기준 사이의 불일치 표(유용한
QC/발견 소재).

> ⚠ **겪었던 문제 — 음수 거리값이 `scipy.linkage`를 크래시시킴**
> fastANI의 (때때로 한 방향만 나오는) 결과를 대칭화하는 초기 코드가 직접
> 만든 수식을 썼는데, 이게 간혹 100%를 살짝 넘는 값을 만들어서 "거리"
> (`100 - ANI`)가 음수가 됐습니다. **해결:**
> `np.nanmean(np.stack([mat, mat.T]))`로 대칭화 — 더 간단하면서 정확합니다.

> ⚠ **겪었던 문제 — read-only 배열에 `np.fill_diagonal`**
> 최신 numpy/pandas에서는 `pandas.DataFrame.fillna(...).values`가
> read-only view를 반환할 수 있습니다. **해결:** in-place numpy 연산 전에
> `.to_numpy(copy=True)`를 호출.

> ⚠ **겪었던 문제 — 엑셀에서 열면 한글이 깨지거나 열 구분이 안 됨**
> 두 가지 별개 문제였습니다: (1) BOM 없는 UTF-8은 한국어 로케일 엑셀에서
> 잘못 해석됨 → `encoding="utf-8-sig"`로 저장. (2) `.tsv`를 엑셀에서
> 더블클릭하면 탭 구분자를 자동 인식 못하는 경우가 많음 → 협업자가 직접
> 엑셀로 열 파일은 콤마 구분 `.csv`를 우선 사용.

---

## 5. 4단계 — ANI 이웃 기반 애매한 종 판정 보정

**스크립트:** `resolve_species_via_ani.py`

BLAST 판정이 속 수준(`"... sp."`)에 그친 시료에 대해, 게놈 전체 ANI가
95% 이상이면서 종명이 명확한 이웃 시료를 찾아 그 종명을 채택합니다.
서로 다른 이웃이 충돌하면(같은 임계값 이상에서 종이 여러 개로 갈리면)
자동판정하지 않고 수동 검토용으로 남겨둡니다 — 이렇게 해서 자동화 단계를
보수적으로 유지했습니다.

```bash
python3 resolve_species_via_ani.py \
    --ani ani_out/ani_result.tsv \
    --master master_table_normalized.tsv \
    --species-col species_normalized \
    --out master_table_final.tsv \
    --ani-threshold 95.0
```

`species_final`, `genus_final`, `resolution_method`
(`original_call` / `ANI_neighbor` / `unresolved` / `ANI_conflict(...)`)
컬럼이 추가됩니다.

---

## 6. 5단계 — QC 제외 플래그

**스크립트:** `flag_qc_excluded.py`

220개 중 5개(2.3%)는 Panaroo에서 심각한 유전자 인식 실패(어노테이션된
CDS 중 약 0~45%만 실제로 인식됨)를 보였는데, 생각할 수 있는 모든
구조적 점검(seqid 일치, 서열 길이 vs CDS 좌표 정합성, 유전자 ID 중복
여부)을 통과했음에도 원인을 못 찾았습니다. 근본 원인 규명 없이 이
5개는 pangenome/계통수 분석에서 제외했고, 그 사유를 조용히 빼는 대신
명시적으로 기록해뒀습니다.

```bash
python3 flag_qc_excluded.py \
    --master master_table_final.tsv --out master_table_qc.tsv
```

> ⚠ **이 프로젝트에서 유일하게 미해결로 남은 이슈입니다.** 재현하실
> 때 이 5개 시료를 Panaroo/Prokka로 다시 돌려서 문제가 재현되는지,
> 아니면 일시적인 문제였는지 재점검해볼 가치가 있습니다.

---

## 7. 6단계 — 속(genus)/기능군별 재편

**스크립트:** `organize_by_group.py`

각 시료의 `.faa`/`.gff`를 그룹 단위 분석을 위해 두 계층 구조로 복사합니다:

- `by_genus/<Genus>/` — 시료 3개 이상인 속만 자체 폴더를 가짐; 그 미만은
  `_minor_genera/`로 몰아넣음.
- `by_functional_group/{LAB, Bacillus_group, Other_Environmental,
  Unresolved}/` — 대략적인 생태/기능 그룹핑 (LAB = 유산균류 속;
  Bacillus_group = 장류 발효에 흔한 포자형성 Bacillaceae/
  Paenibacillaceae 계열 속).

```bash
python3 organize_by_group.py \
    --master master_table_qc.tsv --outdir grouped/ \
    --min-genus-n 3 --mode copy
```

---

## 8. 7단계 — Pangenome 분석 (Panaroo)

**스크립트:** `run_panaroo_all_groups.sh` (환경: `panaroo_env`)

```bash
conda activate panaroo_env
./run_panaroo_all_groups.sh grouped/ pangenome/ 24 3
```

`by_genus/*/genomes_gff`와 `by_functional_group/*/genomes_gff` 폴더를
순회하되, 서로 다른 속이 섞인 그룹(`_minor_genera`, `unresolved`,
`Other_Environmental`)과 이미 완료된 그룹은 건너뜁니다.

```
panaroo -i <그룹>/genomes_gff/*.gff -o <출력> --clean-mode strict -a core -t 24
```

> ⚠ **겪었던 문제 — `core_gene_alignment.aln`이 아예 안 만들어짐**
> Panaroo의 정렬(alignment) 단계는 선택사항입니다. `-a core`(또는
> `-a pan`)를 안 주면 `gene_presence_absence.csv`는 멀쩡해 보여도
> alignment 파일 자체가 생성되지 않습니다. 계통수를 만들 계획이면
> 반드시 `-a core`를 넣으세요.

> ⚠ **겪었던 문제 — 넓은 분류군 그룹에서 `core_genes = 0`**
> 이건 버그가 **아닙니다**. Panaroo/Roary는 "core"를 고정된 서열
> 유사도 클러스터링 기준으로 전체 시료의 99~100%에 존재하는 유전자로
> 정의합니다. 실제로는 여러 종을 포함하는 속 전체(예: `Bacillus`)나,
> 여러 속을 아우르는 기능군(`LAB`, `Bacillus_group`)에 이 기준을
> 적용하면, 진짜 오솔로그라도 서열이 너무 갈라져서 같은 클러스터로 안
> 묶이는 경우가 흔하고, 심지어 genome 단 1개만 불완전/발산해도 "core"
> 카테고리 전체가 0으로 무너질 수 있습니다. **진단 방법:** 내장된
> 카테고리 하나만 믿지 말고, `gene_presence_absence.Rtab`에서 전체
> 존재-비율 분포를 직접 계산해서(§9) 100/99/95/90/80/50% 등 여러
> 기준으로 확인한 뒤에 "공통 유전체가 없다"고 결론 내려야 합니다.
> **후속 영향:** core gene이 0개인 속은 `core_gene_alignment.aln`도
> 없어서 FastTree 계통수를 만들 수 없습니다 — 이 경우 ANI 거리 기반
> 대체 트리를 씁니다(§13).

---

## 9. 7-1단계 — QC 진단 및 유연한 core/accessory 분석

**스크립트:** `diagnose_pangenome_genomes.py`, `analyze_gene_presence.py`,
`summarize_pangenomes.py`

```bash
python3 summarize_pangenomes.py \
    --pangenome-root pangenome/ --out pangenome/pangenome_summary.csv

python3 diagnose_pangenome_genomes.py \
    --grouped-root grouped/ --pangenome-root pangenome/ \
    --out pangenome/genome_qc_report.csv --ratio-threshold 0.5

python3 analyze_gene_presence.py \
    --pangenome-root pangenome/ --outdir pangenome/presence_analysis/
```

`diagnose_pangenome_genomes.py`가 §6의 문제 시료 5개를 찾아낸
스크립트입니다: 시료별로 입력 `.gff`에 실제 어노테이션된 CDS 개수와
Panaroo가 그 genome에 대해 인식한 유전자 패밀리 수를 비교합니다. 격차가
크면(기본 ratio < 0.5) QC 문제로 플래그합니다.

---

## 10. 8단계 — genome별 COG 기능 프로파일

**스크립트:** `extract_eggnog_summary.py`

```bash
python3 extract_eggnog_summary.py \
    --master master_table_qc.tsv \
    --out-count eggnog_cog_count_wide.tsv \
    --out-ratio eggnog_cog_ratio_wide.tsv
```

각 시료의 `annotation_EggNOG.xlsx` → `Eggnog_Count` 시트(Macrogen이 이미
집계해둔 `Eggnog, Description, Count, Ratio (%)`)를 읽어서 시료 x
COG카테고리 wide 매트릭스를 만듭니다. **시료 간 비교에는 Ratio(%) 표를
쓰세요** (genome 크기가 서로 다르므로).

> ⚠ **겪었던 문제 — 시트명이 연도마다 다름**
> 2021년 데이터는 시트명이 `Eggnog_Count`였는데, 2023년 데이터는
> `Eggnog Count`(언더스코어 대신 공백)였습니다. **해결:** 시트명을
> 정확히 일치시키는 대신, 소문자화+공백/언더스코어 제거 후 비교하도록
> 정규화.

> ⚠ **겪었던 문제 — 한 시료 시트 안 카테고리 중복 행이 배치 병합을
> 크래시시킴**
> 일부 `Eggnog_Count` 시트에서 같은 카테고리 코드가 두 행 이상으로
> 나왔습니다. `set_index("Eggnog")`로 단순하게 `pandas.Series`를 만들면
> 인덱스가 고유하지 않게 되고, 이런 Series를 여러 개 합칠 때 크래시가
> 났습니다. **해결:** 합치기 전에 시료별로 `groupby("Eggnog").sum()`.

> ⚠ **겪었던 문제 — 가짜 "카테고리" 컬럼**
> 원본 엑셀의 구분선(`-----------------------------------------------
> ----`)이 진짜 COG 코드인 것처럼 집계됐습니다. 이후 모든 비교
> 스크립트에서 작은 `NON_CATEGORY_COLS` 집합(`{"Total", "-",
> "----...----"}`)으로 명시적으로 제외합니다.

---

## 11. 9단계 — 기능 프로파일 통계 비교

**스크립트:** `compare_functional_groups.py` (LAB vs Bacillus_group,
2그룹), `genus_cog_heatmap.py` (전체 속, Kruskal-Wallis)

```bash
python3 compare_functional_groups.py \
    --ratio-tsv eggnog_cog_ratio_wide.tsv \
    --outdir functional_comparison/ \
    --group-a LAB --group-b Bacillus_group

python3 genus_cog_heatmap.py \
    --ratio-tsv eggnog_cog_ratio_wide.tsv \
    --outdir genus_cog_analysis/ --min-n 3
```

카테고리별 Mann-Whitney U 검정(2그룹) / Kruskal-Wallis(3그룹 이상)에
Benjamini-Hochberg FDR 보정을 적용하고, 그룹별 막대그래프, 220개 전체
COG 프로파일 PCA, 속 x 카테고리 히트맵을 만듭니다.

> ⚠ **겪었던 문제 — `NaN` p-value 하나가 모든 카테고리의 FDR 보정을
> 조용히 망가뜨림**
> 두 그룹 모두 분산이 0인 COG 카테고리(예: 세균에 아예 없는 카테고리
> `Y`)에서 `scipy.stats.mannwhitneyu`가 `p = NaN`을 반환했습니다. 직접
> 구현한 Benjamini-Hochberg 코드가 원본 p-value 배열 전체에
> `np.argsort`/`np.minimum.accumulate`를 적용하는데, `NaN` 하나가
> 배열 전체의 정렬 순서를 오염시켜서, 원본 p-value가 `1e-30`처럼
> 극히 작아도 보정된 p-value가 전부 "유의하지 않음"으로 나왔습니다.
> **해결:** BH 보정 전에 `NaN`을 `1.0`(비유의로 간주)으로 치환하고,
> 검정 호출 자체도 `try/except ValueError`로 감싸서 완전히 퇴화된
> 입력도 잡아냄.

---

## 12. 10단계 — 항생제내성/병원성 유전자 스크리닝

**스크립트:** `run_abricate_batch.sh` (환경: `abricate_env`),
`compare_resistance_virulence.py`

```bash
conda activate abricate_env
./run_abricate_batch.sh genomes_fna/ abricate_out/
# abricate_out/card_summary.tsv, vfdb_summary.tsv 생성됨

conda activate compgenomics
python3 compare_resistance_virulence.py \
    --summary-tsv abricate_out/card_summary.tsv \
    --master master_table_qc.tsv --db-label CARD \
    --outdir resistance_comparison/ --group-a LAB --group-b Bacillus_group
# vfdb_summary.tsv / --db-label VFDB로 동일하게 반복
```

유전자별 Fisher's exact test(존재/부재) + BH-FDR, 그리고 그룹간 총
"부담"(genome당 검출 유전자 수) Mann-Whitney U 비교.

> ⚠ **겪었던 문제 — pandas/matplotlib API 변경으로 실행 도중 크래시**
> 최신 pandas에서는 `DataFrame.applymap()`이 제거됐습니다(대신
> `.map()` 사용). `Axes.boxplot(..., labels=...)`도 최근 matplotlib에서
> `tick_labels=`로 이름이 바뀌었습니다. 둘 다 한 줄짜리 수정이지만,
> 이 코드를 그대로 가져다 쓰신다면 설치된 버전을 확인해보세요.

---

## 13. 11단계 — antiSMASH를 이용한 이차대사산물(BGC) 스크리닝

**스크립트:** `merge_gbk_files.sh`, `run_antismash_batch.sh` (환경:
`antismash_env`), `analyze_antismash_bgcs.py`

```bash
# 1) contig별 GenBank 파일을 시료 단위 multi-record 파일로 병합
#    (consolidate_wgs.sh와 같은 소스 트리 탐색 방식 재사용)
./merge_gbk_files.sh \
    "/mnt/f/WGS_Results/#Whole_Genome_Sequencing_Macrogen" \
    "/mnt/f/WGS_Consolidated"

# 2) genome별 antiSMASH 실행, 기존 어노테이션 재사용
conda activate antismash_env
./run_antismash_batch.sh genomes_gbk/ antismash_out/ master_table_qc.tsv \
    24 "LAB,Bacillus_group"
./run_antismash_batch.sh genomes_gbk/ antismash_out/ master_table_qc.tsv \
    24 "Other_Environmental,Unresolved"   # 나머지 시료

# 3) JSON 결과를 시료 x BGC종류 매트릭스 + 통계로 집계
conda activate compgenomics
python3 analyze_antismash_bgcs.py \
    --antismash-root antismash_out/ --master master_table_qc.tsv \
    --outdir antismash_analysis/ --group-a LAB --group-b Bacillus_group
```

`--genefinding-tool none`은 유전자를 새로 예측하는 대신 기존
Prokka/tbl2asn 유전자 콜을 재사용합니다 — 더 빠르고, 파이프라인 나머지
부분과 유전자 번호 체계가 일관되게 유지됩니다. genome당 실행시간은
약 10~70초 (220개 전체 기준, 24스레드로 약 2~3시간).

실제로 사용된 antiSMASH 7 JSON 스키마:
`records[].areas[].products` (BGC 영역별 문자열 리스트) — 이건 문서나
기억에 의존해 스키마를 가정하는 대신, 배치 파서를 작성하기 전에 실제
시료 1개의 출력을 직접 열어서 확인한 것입니다.

> ⚠ **겪었던 문제 — 배치 필터링 스크립트가 `antismash_env`에 없는
> pandas에 의존함**
> `run_antismash_batch.sh`의 그룹 필터링 로직이 원래 pandas를 import하는
> 파이썬 한 줄짜리를 호출했습니다. `antismash_env`엔 기본적으로 pandas가
> 없습니다. **해결:** 필터링 로직을 순수 `awk`로 재작성(외부 의존성 없음,
> 단순 TSV 필터링이라 더 빠르기도 함).

> ⚠ **겪었던 문제 — UTF-8 BOM이 awk 컬럼 탐색을 깨뜨려 조용히 전체
> 줄을 출력함**
> `master_table_qc.tsv`가 `encoding="utf-8-sig"`(BOM)로 저장돼 있었는데,
> 이게 헤더 첫 토큰(`sample_id`)에 눈에 안 보이는 접두 문자를
> 붙입니다. awk의 `$i == "sample_id"` 비교가 매칭이 안 돼서 컬럼
> 인덱스가 설정되지 않은 채(awk에서는 `0`) 남았고, `print $0`
> (전체 줄)이 ID 대신 출력됐습니다 — 이 때문에 종명, 파일 경로 등이
> 쉘의 단어분리(word-splitting)를 거치면서 각각 별개의 "시료"처럼
> 처리돼 뒤죽박죽인 출력이 나왔습니다. **해결:** `sample_id`를
> (이 파이프라인의 모든 표에서 구조상 항상 그런 대로) 1번째 컬럼으로
> 위치 고정하고, 그 컬럼에 한해서는 텍스트 매칭에 의존하지 않음.

> ⚠ **겪었던 문제(사소함, 미수정) — 모든 작업이 성공했는데도 최종
> 요약 카운터가 0/0으로 나옴**
> `printf '%s\n' "$LIST" | while read ...; do COUNTER=$((COUNTER+1));
> done` 형태는 파이프(`|`) 때문에 루프 본문이 서브쉘에서 실행되어,
> `COUNTER` 증가분이 루프 종료 후 부모 쉘에 반영되지 않습니다. 실제
> antiSMASH 실행은 전부 정상 완료됐습니다(시료별 로그로 확인). 마지막에
> 출력되는 "N개 신규/M개 스킵" 요약 줄만 틀렸습니다. 사소한 문제로
> 남겨뒀습니다 — 요약 줄을 믿는 대신
> `find <출력폴더> -name index.html | wc -l`로 실제 완료 개수를
> 확인하는 걸 권장합니다.

---

## 14. 12단계 — 통합 그림

### 14a. 속별 계통수 + COG 히트맵

**스크립트:** `genus_tree_heatmap.py`

```bash
for GENUS in Bacillus Enterococcus Lactiplantibacillus Weissella \
             Latilactobacillus Levilactobacillus Lactococcus Pediococcus \
             Leuconostoc Paenibacillus Staphylococcus Oceanobacillus; do
    python3 genus_tree_heatmap.py \
        --genus "$GENUS" \
        --pangenome-root pangenome/ \
        --ratio-tsv eggnog_cog_ratio_wide.tsv \
        --master master_table_qc.tsv \
        --ani-matrix ani_analysis/ani_matrix.csv \
        --outdir tree_heatmap/
done
```

`core_gene_alignment.aln`으로 FastTree 계통수(`-nt` 모드)를 만들고, 그
옆에 시료별 COG 히트맵을 나란히 그립니다. `--ani-matrix`는 core gene이
0개인 속(여기서는 `Bacillus`, `Lactococcus`)의 fallback용으로 필요하며,
이 경우 ANI 거리 기반 계층적 클러스터링 트리가 대신 그려지고 그림
제목에 그 사실이 명확히 표시됩니다.

### 14b. 전체 데이터셋 통합 그림

**스크립트:** `global_tree_heatmap.py`

```bash
python3 global_tree_heatmap.py \
    --ani-matrix ani_analysis/ani_matrix.csv \
    --master master_table_qc.tsv \
    --antismash-summary antismash_analysis/antismash_bgc_summary.tsv \
    --card-summary abricate_out/card_summary.tsv \
    --vfdb-summary abricate_out/vfdb_summary.tsv \
    --outdir global_summary/
```

220개 전체 genome을 ANI 기반 덴드로그램(기능군별 색상 tip 라벨)으로
그리고, 그 옆에 정규화된 히트맵으로 antiSMASH 총 BGC 개수, CARD 유전자
개수, VFDB 유전자 개수를 붙입니다. 이 그림 하나가 어느 개별 분석
결과만으로는 드러나지 않았을 아웃라이어 시료(CARD=26, VFDB=62로 다른
모든 시료를 압도하는 *Klebsiella pneumoniae* 오염 시료)를 찾아냈습니다.

---

## 15. 핵심 생물학적 발견 요약

1. **LAB와 Bacillus_group은 COG 프로파일, CARD, antiSMASH 결과 전반에
   걸쳐 일관되게 다층적인 기능적 분화를 보입니다**: LAB genome은
   탄수화물대사/번역/복제-복구가 강화되어 있고(빠르고 단순한 발효
   전략), Bacillus_group genome은 운동성, 신호전달, 훨씬 다양한
   이차대사산물 화학무기고(NRPS/PKS/철분획득계; genome당 BGC 중앙값
   13 vs 4, p≈3×10⁻²⁸)와 더 높은 내재적 항생제내성 유전자 부담을
   갖고 있습니다.
2. Bacillus_group 일부 genome은 *B. cereus* 계열 비용혈성 장독소
   유전자(`nheA/B/C`)를 보유하고 있습니다 — 식품안전 관점에서 수동
   후속 확인이 필요한 발견입니다.
3. 한 시료(`HN00179262_F4055`)가 BLAST 결과 *Klebsiella pneumoniae*로
   판정됐고, 통합 CARD/VFDB 히트맵에서 극단적인 아웃라이어로 두드러졌습니다
   — 원료/공정 어느 단계에서 유입됐는지 역추적할 가치가 있는 오염
   사례로 보입니다.

---

## 16. 스크립트 목록

| 단계 | 스크립트 | 환경 |
|---|---|---|
| 통합 | `consolidate_wgs.sh`, `extract_species_from_blast.py` | compgenomics |
| 분류학 | `normalize_taxonomy.py`, `resolve_species_via_ani.py`, `flag_qc_excluded.py` | compgenomics |
| ANI | `run_fastani.sh`, `analyze_ani.py` | compgenomics |
| 그룹화 | `organize_by_group.py` | compgenomics |
| Pangenome | `run_panaroo_all_groups.sh` | panaroo_env |
| Pangenome QC | `summarize_pangenomes.py`, `diagnose_pangenome_genomes.py`, `analyze_gene_presence.py` | compgenomics |
| COG 기능 | `extract_eggnog_summary.py`, `compare_functional_groups.py`, `genus_cog_heatmap.py` | compgenomics |
| 내성/병원성 | `run_abricate_batch.sh`, `compare_resistance_virulence.py` | abricate_env / compgenomics |
| 이차대사산물 | `merge_gbk_files.sh`, `run_antismash_batch.sh`, `analyze_antismash_bgcs.py` | antismash_env / compgenomics |
| 통합 그림 | `genus_tree_heatmap.py`, `global_tree_heatmap.py` | compgenomics |

---

## 17. 이 파이프라인 재현을 위한 일반 교훈

- **엑셀 헤더 행/시트명을 절대 가정하지 마세요** — Macrogen 자체 산출물
  포맷이 2021년과 2023년 배치 사이에 바뀌었습니다(BLAST 헤더 행 오프셋,
  eggNOG 시트명 공백 여부). 배치 파서를 짜기 전에 항상 실제 파일 1개를
  먼저 확인하세요.
- **어노테이션(GFF/GBK)과 서열(FASTA) 파일 간 ID 일치 여부가 downstream
  툴(Panaroo, antiSMASH 등)의 가장 흔한 조용한 실패 지점**입니다.
  multi-contig 어노테이션을 병합할 때는 항상 각 contig의 어노테이션을
  별도로 만들어진 전체 genome FASTA가 아니라 **그 contig 자신의 서열
  파일**과 짝지으세요.
- **서로 충돌하는 툴 의존성은 별도 conda 환경으로 분리**하세요 —
  한 환경에 억지로 다 넣으려 하지 마세요 (Panaroo는 Python ≤3.11이
  필요했는데, 나머지 파이프라인은 자동으로 Python 3.14를 잡았습니다).
- **분류학적으로 넓은 그룹을 비교할 때는 pangenome 툴의 내장 "core"
  정의를 액면 그대로 믿지 마세요** — "공통 유전자가 없다"고 결론 내리기
  전에 항상 전체 유전자 존재-빈도 분포를 다시 계산해서 여러 임계값으로
  확인하세요.
- **다중검정보정 코드는 반드시 NaN p-value를 명시적으로 처리**해야
  합니다(분산이 0인 그룹에서 나옴) — 안 그러면 퇴화된 검정 하나가
  배치의 나머지 모든 결과를 조용히 무효화시킬 수 있습니다.
- **엑셀에서 열릴 산출물은 `encoding="utf-8-sig"`로 저장**하고, 협업자가
  더블클릭으로 열 파일은 탭 `.tsv`보다 콤마 `.csv`를 우선하세요. 그리고
  UTF-8 BOM이 헤더 행에서 순진한 `awk`/쉘의 컬럼-위치 로직을 깨뜨릴 수
  있다는 걸 유의하세요 — ID 컬럼은 텍스트 매칭 대신 위치를 고정하는 게
  안전합니다.
