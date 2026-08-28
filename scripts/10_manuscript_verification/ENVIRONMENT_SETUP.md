# 환경 설정 안내 (투고 전 자체 검증 과정에서 추가된 것)

기존 파이프라인 환경(`compgenomics`, `panaroo_env`, `antismash_env`, `abricate_env`)에
더해, 원고를 투고 전에 자체적으로 재검증하는 과정에서 아래 2개 환경이 새로
필요해졌습니다.

## 1) ncbi_datasets (신규)

용도: `fetch_reference_genomes.sh` — 각 판정 종의 NCBI reference/representative
genome을 다운로드하여 독립적인 taxonomy 검증(원고 2.4, 3.3절)에 사용.

```bash
mamba create -n ncbi_datasets -c conda-forge ncbi-datasets-cli -y
```

## 2) mlst_env (신규)

용도: `characterize_kpneumoniae_outlier.sh`의 MLST 부분 — 이상치로 확인된
Klebsiella pneumoniae 시료의 sequence type 판정(원고 2.7, 3.8절).

```bash
mamba create -n mlst_env -c bioconda -c conda-forge mlst -y
```

같은 스크립트의 PlasmidFinder 부분은 기존 `abricate_env`에서 실행합니다
(추가 설치 불필요, `abricate --db plasmidfinder`로 DB만 자동 다운로드됨).

## 3) 기존 compgenomics 환경에 추가로 필요한 패키지

`check_clonal_pairs.py`, `select_dereplicated_representatives.py`가 클론 그룹을
그래프의 connected component로 찾는 데 `networkx`를 사용합니다. 기존
`compgenomics` 환경에 없다면 추가 설치해주세요:

```bash
conda activate compgenomics
pip install networkx --break-system-packages
```

(scipy, scikit-learn, pandas, numpy, matplotlib은 기존 파이프라인 스크립트들이
이미 사용하고 있어 `compgenomics` 환경에 원래부터 있어야 정상입니다.)

## 환경별 스크립트 매핑 요약

| 환경 | 스크립트 |
|---|---|
| `compgenomics` (+ networkx) | 대부분의 `.py` 스크립트 (아래 표 외 전부) |
| `ncbi_datasets` (신규) | `fetch_reference_genomes.sh` |
| `panaroo_env` | `rerun_lactococcus_enterococcus_pangenome.sh`, `run_species_level_pangenome.sh` |
| `abricate_env` | `characterize_kpneumoniae_outlier.sh` (PlasmidFinder 부분) |
| `mlst_env` (신규) | `characterize_kpneumoniae_outlier.sh` (MLST 부분) |
