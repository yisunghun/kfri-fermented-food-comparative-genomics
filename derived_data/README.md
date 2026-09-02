# derived_data

Consolidated, analysis-ready summary tables underlying the figures and
statistics reported in the manuscript. These are outputs of the
pipeline scripts in `scripts/`, provided here for reproducibility and
reanalysis without needing to re-run the full pipeline on raw genome
data (which is not publicly deposited; see manuscript Data
availability statement).

| File | Description | Produced by |
|---|---|---|
| `Supplementary_Table_S1_genome_QC_metrics.tsv` | Per-genome assembly QC metrics and final validated taxonomy (220 isolates) | `scripts/12_presubmission_review/build_s1_definitive.py` |
| `ani_matrix.csv` | Pairwise whole-genome ANI matrix (all 220 genomes) | `scripts/03_ani/` |
| `eggnog_cog_ratio_wide_UNIFORM.tsv` | Per-genome COG category relative abundance (%), uniform eggNOG-mapper v2.1.15 re-annotation | `scripts/06_functional_cog/` |
| `eggnog_cog_count_wide_UNIFORM.tsv` | Per-genome COG category gene counts, uniform re-annotation | `scripts/06_functional_cog/` |
| `card_summary.tsv` | Per-genome CARD antibiotic-resistance gene screening results | `scripts/07_resistance_virulence/` |
| `vfdb_summary.tsv` | Per-genome VFDB virulence-gene screening results | `scripts/07_resistance_virulence/` |
| `antismash_bgc_summary.tsv` | Per-genome antiSMASH biosynthetic gene cluster counts | `scripts/08_secondary_metabolites/` |
| `mismatch_resolution_summary.tsv` | 72-genome reference-panel ANI comparison for isolates below the single-reference 95% ANI threshold (Section 3.3, Supplementary Table S4 Part A) | `scripts/10_manuscript_verification/resolve_ani_mismatches.py` |
| `gtdbtk.bac120.summary.tsv` | Full GTDB-Tk v2.4.1 classification output for all 220 genomes (GTDB release R226) | `scripts/02_taxonomy/` (GTDB-Tk `classify_wf`) |
