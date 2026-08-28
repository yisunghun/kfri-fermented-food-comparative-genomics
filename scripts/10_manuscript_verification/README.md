# 10_manuscript_verification

투고 전, 저자 자체적으로 원고를 재검토하며 taxonomy 검증 방법과 통계적
견고성을 보완하기 위해 만든 스크립트들입니다. (외부 심사자의 리뷰 코멘트가
아니라, 투고 전 내부 품질 점검 과정에서 자체적으로 확인한 사항들입니다.)

## Taxonomy 독립 검증 (원고 2.4, 3.3절)

| 스크립트 | 목적 |
|---|---|
| `fetch_reference_genomes.sh` | 판정된 각 종의 NCBI reference/representative genome 다운로드 |
| `validate_species_via_reference_ani.py` | 각 isolate를 자기 종의 reference genome과 1:1 ANI 비교 (독립 검증) |
| `resolve_ani_mismatches.py` | 검증 실패한 isolate를 72개 reference 전체와 비교해 실제 근연종 탐색 |

기존 파이프라인의 species 검증은 220개 study genome끼리 서로 비교하는
방식이라 순환논리였음을 자체적으로 확인, 외부 NCBI reference genome과의
독립 비교로 보완.

## 클론/비독립 isolate 및 민감도 분석 (3.9, 3.10절)

| 스크립트 | 목적 |
|---|---|
| `check_clonal_pairs.py` | ANI 행렬에서 근접(잠재적 클론) 시료쌍 탐지 |
| `select_dereplicated_representatives.py` | 클론 그룹당 대표 시료 1개 선정 |
| `filter_to_dereplicated.py` | 기존 데이터 파일들을 대표 시료만 남기도록 필터링 |
| `genus_stratified_sensitivity.py` | isolate 대신 genus를 분석단위로 한 민감도 분석 (특정 속의 과대표성 문제 점검) |

## 통계 방법론 보완 (2.9, 3.5–3.7절)

| 스크립트 | 목적 |
|---|---|
| `normalize_by_genome_size.py` | CARD/VFDB/antiSMASH burden을 genome 크기(Mb)로 정규화 재검정 |
| `add_effect_sizes.py` | 모든 Mann-Whitney U 비교에 rank-biserial correlation(효과크기) 추가 |
| `cog_pca_clr.py` | COG 비율 데이터에 CLR 변환 후 PCA 재검증 (compositional data) |
| `check_nheABC_denominator.py` | nheABC 양성 시료와 B. cereus 종복합체의 교차 확인 (분모 재검토) |

## VFDB 방향성 재검증

| 스크립트 | 목적 |
|---|---|
| `check_vfdb_prevalence.py` | VFDB 양성 시료 비율 및 분포 확인 |
| `verify_vfdb_direction.py` | dereplicated/genus-stratified 데이터에서 VFDB 방향(LAB vs Bacillus-group)을 단측검정으로 명확히 재확인 |

초기에 "Bacillus-group이 VFDB burden이 더 높다"고 서술했던 것이 평균값만
비교한 데서 온 오류였음을 자체적으로 발견 (Mann-Whitney는 순위 기반 검정이라
평균 비교로 방향을 추론하면 안 됨). 위 스크립트들로 재검증하여 실제로는
LAB(주로 Enterococcus)이 순위 기반으로 유의하게 더 높다는 것으로 정정.

## 환경

`ENVIRONMENT_SETUP.md` 참고 (신규 conda 환경 `ncbi_datasets`, 기존
`compgenomics`에 `networkx` 패키지 추가 필요).
