# 11_supplementary_analysis

투고 전 자체 검증 과정에서 추가로 수행한 보충 분석 스크립트들입니다.
(외부 심사자의 리뷰 코멘트가 아니라, 투고 전 내부 품질 점검 과정에서
저자가 자체적으로 확인하고 보완한 사항들입니다.)

| 스크립트 | 목적 | 관련 절 |
|---|---|---|
| `assess_genome_qc.py` | 220개 genome의 N50/contig수/GC% 등 assembly QC 지표 산출 | 2.1, 3.1, Supp. Table S1 |
| `rerun_lactococcus_enterococcus_pangenome.sh` | Taxonomy 재동정(10_manuscript_verification) 반영 후 Lactococcus/Enterococcus pangenome 재계산 | 3.3, 3.4 |
| `identify_species_pangenome_candidates.py` | 종 수준 pangenome 분석에 적합한 후보 종(시료수 충분한 종) 탐색 | 3.4 |
| `run_species_level_pangenome.sh` | 지정 종(B. velezensis, L. plantarum)에 대해 species-level pangenome 실행 | 3.4, Table 6 |
| `cog_boxplot_figure.py` | COG 상위 10개 카테고리를 box plot으로 재작업 (Mann-Whitney 비모수 검정과 시각적으로 일치시킴) | 3.5, Fig. 1b |
| `characterize_kpneumoniae_outlier.sh` | K. pneumoniae 이상치 시료의 MLST(ST 판정) 및 plasmid replicon 확인 | 2.7, 3.8 |

## 09_integrative_figures/global_tree_heatmap.py (수정본)

같은 커밋에서 기존 파이프라인 스크립트도 아래 이유로 함께 수정:
- gridspec 서브플롯 경계에 균주명 라벨이 잘리는 문제 수정 (`clip_on=False`, 왼쪽 여백 확보)
- 라벨을 `short_id` 기준으로 축약 (`set_yticklabels()`로 안전하게 교체 — `.set_text()`만 쓰면
  matplotlib이 draw 시점에 라벨을 재생성하며 원래 값으로 되돌아가는 문제가 있어 수정)
- `constrained_layout=True`로 전면 재구성 (라벨-히트맵 겹침 문제 해결)
- 히트맵 간격을 라벨에 더 가깝게 조정 (`fig.get_layout_engine().set(...)`)

## 환경

`10_manuscript_verification/ENVIRONMENT_SETUP.md` 참고 (신규 conda 환경
`mlst_env` — `characterize_kpneumoniae_outlier.sh`의 MLST 부분에 필요).
PlasmidFinder 부분은 기존 `abricate_env`, pangenome 스크립트들은 기존
`panaroo_env`를 그대로 사용합니다.
