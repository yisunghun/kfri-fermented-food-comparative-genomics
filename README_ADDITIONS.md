## Review-response additions (new since initial submission)

Two new script folders were added to address issues identified during manuscript verification, 
alongside one bug-fixed file in the original `09_integrative_figures/` folder.

```
scripts/
  09_integrative_figures/
    global_tree_heatmap.py          UPDATED: fixed tip-label rendering (matplotlib
                                     redraw bug), tip-label/heatmap overlap, and
                                     switched to constrained_layout for reliable
                                     multi-panel spacing at n=220.

  10_manuscript_verification/    	Independent validation & sensitivity analyses
    fetch_reference_genomes.sh          Download one NCBI reference genome per
                                         assigned species (for independent ANI check)
    validate_species_via_reference_ani.py  1:1 fastANI of each isolate vs. its own
                                         assigned species' reference genome
    resolve_ani_mismatches.py           For isolates failing the above, compare
                                         against the full reference panel to find
                                         the best-supported alternative species/genus
    check_clonal_pairs.py               Identify near-identical (ANI >= threshold)
                                         isolate pairs/groups from the ANI matrix
    select_dereplicated_representatives.py  Pick one representative isolate per
                                         near-identical group
    filter_to_dereplicated.py           Filter existing data tables down to the
                                         dereplicated representative set
    genus_stratified_sensitivity.py     Re-test LAB vs. Bacillus-group comparisons
                                         with genus (not isolate) as the unit of
                                         replication
    normalize_by_genome_size.py         Recompute CARD/VFDB/antiSMASH burden
                                         normalized per Mb of genome size
    check_nheABC_denominator.py         Cross-tabulate nheABC-positive isolates
                                         against the B. cereus sensu lato species
                                         complex membership
    add_effect_sizes.py                 Rank-biserial correlation (effect size) for
                                         every Mann-Whitney comparison
    cog_pca_clr.py                      CLR-transformed PCA of COG profiles
                                         (compositional-data-appropriate re-analysis)
    check_vfdb_prevalence.py            Diagnostic: VFDB gene-count distribution and
                                         positivity rate per functional group
    verify_vfdb_direction.py            Diagnostic: one-sided Mann-Whitney direction
                                         check for VFDB burden (dereplicated and
                                         genus-stratified subsets)

  11_supplementary_analysis/  			Supplementary tables/figures for review response
    assess_genome_qc.py                 Per-genome assembly QC metrics (N50, L50,
                                         contig count, GC%) -> Supplementary Table S1
    identify_species_pangenome_candidates.py  Rank species by isolate count to pick
                                         species-level pangenome candidates
    run_species_level_pangenome.sh      Build a species-level (not genus-level)
                                         Panaroo pangenome for a chosen species
    rerun_lactococcus_enterococcus_pangenome.sh  Re-run Panaroo for Lactococcus/
                                         Enterococcus after the 2-isolate genus
                                         correction (Section 3.3)
    cog_boxplot_figure.py               Box-plot (median/IQR) version of the
                                         top-effect-size COG categories -> Fig. 1b
    characterize_kpneumoniae_outlier.sh  MLST + PlasmidFinder screening for the
                                         K. pneumoniae outlier isolate
```

See `docs/PROTOCOL.md` for the original 22-script pipeline; the scripts above were
added specifically during manuscript verification and are documented inline (each
script's docstring/header explains its purpose and exact usage).
