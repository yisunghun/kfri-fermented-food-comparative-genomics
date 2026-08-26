# Comparative Genomics Pipeline for Fermented-Food Bacterial Isolates (n=220)

A reproducible protocol for comparing ~220 whole-genome-sequenced bacterial
isolates from Korean fermented foods (kimchi, jang, etc.), from raw Macrogen
deliverables through species assignment, ANI clustering, pangenome analysis,
functional (COG) profiling, antibiotic-resistance/virulence screening,
secondary-metabolite (BGC) screening, and final integrative tree+heatmap
figures.

This document only describes the **final, working** version of each step.
Where earlier attempts failed, a short "⚠ Pitfall encountered" box explains
what went wrong and how it was fixed, so the same mistakes can be avoided.

---

## 0. Input Data Layout

Sequencing/annotation was delivered by Macrogen with the following structure
(one folder per sequencing year):

```
<SRC_ROOT>/<year>/Analysis_Data_Done/<sample_folder>/<assembly_dir>/
    consensus.fasta                     # whole-genome assembly (all contigs)
    <short_id>_BLAST.xlsx               # BLAST hit per contig (sheet 'Result')
    <assembly_dir>/contig1/ contig2/ …  # PER-CONTIG annotation (Prokka/tbl2asn-style)
        contigN.gff  contigN.faa  contigN.ffn  contigN.fna  contigN.gbk

<sample_folder>/<short_id>_FunctionalAnnotation/FunctionalAnnotation/
    annotation_EggNOG.xlsx              # sheet 'Eggnog_Count' (or 'Eggnog Count')
```

**Important structural fact used throughout the pipeline:** each contig's
`.gff`, `.faa`, `.ffn`, `.fna`, `.gbk` were produced *together* by the same
annotation run, so a contig's `.fna` header is **guaranteed** to match the
`seqid` used in that contig's own `.gff`. The top-level `consensus.fasta`,
however, was produced independently (by the assembler) and its headers are
**not** guaranteed to match. This distinction caused a major bug (see §3.1)
and is the single most important lesson from this project.

---

## 1. Environment Setup

Four separate conda/mamba environments were required because of conflicting
Python-version and dependency requirements between tools.

```bash
# 1) Main analysis environment (fastANI, pandas/numpy/scipy/matplotlib,
#    scikit-learn, FastTree, Biopython)
mamba create -n compgenomics --override-channels -c bioconda -c conda-forge \
    fastani pandas scipy matplotlib numpy openpyxl -y
conda activate compgenomics
pip install scikit-learn --break-system-packages
mamba install -c bioconda -c conda-forge fasttree biopython -y

# 2) Panaroo (pangenome) - needs an OLDER Python; conflicts with compgenomics'
#    auto-resolved Python 3.14 (python-edlib doesn't support 3.14 yet)
mamba create -n panaroo_env --override-channels -c bioconda -c conda-forge \
    python=3.10 panaroo -y

# 3) antiSMASH (secondary metabolite BGC detection)
mamba create -n antismash_env --override-channels -c bioconda -c conda-forge \
    antismash -y
conda activate antismash_env
download-antismash-databases
antismash --check-prereqs

# 4) abricate (CARD / VFDB screening)
mamba create -n abricate_env --override-channels -c bioconda -c conda-forge \
    abricate -y
conda activate abricate_env
abricate --setupdb
```

> ⚠ **Pitfall — slow/unstable `mamba create`**
> Without `--override-channels`, mamba also indexes Anaconda's default
> `pkgs/main` / `pkgs/r` channels (commercial, unnecessary here), which can
> take several minutes and occasionally fails to resolve
> ("environment specs not solvable"). Always pass `--override-channels`
> with `-c bioconda -c conda-forge` explicitly.

Korean-font rendering for matplotlib (only needed if you keep Korean text in
plots — this pipeline's final outputs are English):

```bash
sudo apt-get install -y fonts-nanum
rm -rf ~/.cache/matplotlib   # force matplotlib to re-scan fonts
```

---

## 2. Step 1 — Consolidate per-contig annotation into per-genome files

**Script:** `consolidate_wgs.sh`

For each sample, merges all `contigN.faa/.ffn/.gff` files into one
per-genome file, copies the assembly to a flat `genomes_fna/` directory, and
builds `master_table.tsv` (one row per sample) with a BLAST-derived species
call and a contamination flag.

```bash
./consolidate_wgs.sh \
    "/mnt/f/WGS_Results/#Whole_Genome_Sequencing_Macrogen" \
    "/mnt/f/WGS_Consolidated"
```

Key design points in the final version:

- **Species call comes from `*_BLAST.xlsx`, not from the sample's PDF
  report filename.** The PDF is a curated summary and can disagree with the
  raw top BLAST hit.
- The `Result` sheet's real header row is **row index 1** (0-indexed),
  not row 0 — Excel's merged header cells push the actual column names
  (`Name, Q_Length, …, Description, …`) down one row.
- Species is read from the `Description` field of the **longest
  non-plasmid contig** (chromosome), not just the first BLAST hit — small
  plasmid contigs are frequently horizontally-transferred and can point to a
  different species/genus without indicating contamination.
- A genome is flagged `contam_flag=yes` only if contigs **≥20 kb**
  (`--contam-min-length`, tunable) disagree on species — this excludes
  small plasmids from the contamination call.
- The merged `##FASTA` block appended to each sample's `.gff` is built by
  concatenating **each contig's own `.fna`**, never the top-level
  `consensus.fasta` (see the pitfall box below).

```bash
# species extraction is delegated to:
python3 extract_species_from_blast.py <sample>_BLAST.xlsx --verbose
```

> ⚠ **Pitfall — GFF/FASTA header mismatch broke Panaroo (and would have
> broken antiSMASH too)**
> Original code appended the top-level `consensus.fasta` after `##FASTA`.
> Its headers (e.g. `>contig1`) did **not** match the `seqid` used inside
> each contig's own `.gff` (e.g. `gnl|MG|AMT60212_1_1`) — different naming
> systems from the assembler vs. the annotator. Panaroo then failed with
> `Invalid gene sequence!` for every gene. **Fix:** append each contig's own
> `.fna` (paired 1:1 with its `.gff`), guaranteeing `seqid` consistency
> regardless of naming convention.

> ⚠ **Pitfall — over-sensitive contamination flag**
> Initially every contig (including tiny plasmids) was compared, flagging
> ~75% of genomes as "contaminated". Restricting the comparison to contigs
> ≥20 kb dropped this to a realistic level.

---

## 3. Step 2 — Taxonomy normalization

**Script:** `normalize_taxonomy.py`

Maps pre-2020 `Lactobacillus` species names to their post-reclassification
genera (Zheng et al. 2020), e.g. `Lactobacillus plantarum` →
`Lactiplantibacillus plantarum`. Necessary because the BLAST reference
database mixes old and new nomenclature for genomically-identical hits,
which otherwise looks like a false ANI/taxonomy mismatch downstream.

```bash
python3 normalize_taxonomy.py \
    --master master_table.tsv --out master_table_normalized.tsv
```

---

## 4. Step 3 — Average Nucleotide Identity (ANI)

**Scripts:** `run_fastani.sh`, `analyze_ani.py`

```bash
./run_fastani.sh genomes_fna/ ani_out/ 8

python3 analyze_ani.py \
    --ani ani_out/ani_result.tsv \
    --master master_table_normalized.tsv \
    --outdir ani_analysis/ \
    --species-col species_normalized
```

Produces: symmetric ANI matrix (CSV), a genus-colored dendrogram (PDF), and
a table of disagreements between the BLAST-based species call and the
95%-ANI species threshold (useful QC/discovery artifact).

> ⚠ **Pitfall — negative distances crashed `scipy.linkage`**
> A first attempt to symmetrize the (sometimes one-directional) fastANI
> output used a hand-rolled formula that produced values slightly above
> 100%, which then produced negative "distances" (`100 - ANI`). **Fix:**
> symmetrize with `np.nanmean(np.stack([mat, mat.T]))`, which is both
> simpler and correct.

> ⚠ **Pitfall — `np.fill_diagonal` on a read-only array**
> `pandas.DataFrame.fillna(...).values` can return a read-only view in
> newer numpy/pandas. **Fix:** `.to_numpy(copy=True)` before any in-place
> numpy operation.

> ⚠ **Pitfall — garbled Korean text / broken columns when opened in Excel**
> Two separate issues: (1) UTF-8 without a BOM is misread by Excel under a
> Korean locale → save with `encoding="utf-8-sig"`. (2) Double-clicking a
> `.tsv` in Excel does not reliably auto-detect the tab delimiter →
> prefer comma-separated `.csv` for anything a collaborator will open
> directly in Excel.

---

## 5. Step 4 — Resolve ambiguous species calls via ANI neighbors

**Script:** `resolve_species_via_ani.py`

For samples whose BLAST call was only genus-level (`"... sp."`), looks for
a genome-wide ANI ≥95% neighbor with an unambiguous species call and adopts
it. Conflicting neighbors (multiple different resolved species above
threshold) are left unresolved and flagged for manual review — this keeps
the automated step conservative.

```bash
python3 resolve_species_via_ani.py \
    --ani ani_out/ani_result.tsv \
    --master master_table_normalized.tsv \
    --species-col species_normalized \
    --out master_table_final.tsv \
    --ani-threshold 95.0
```

Adds `species_final`, `genus_final`, `resolution_method`
(`original_call` / `ANI_neighbor` / `unresolved` / `ANI_conflict(...)`).

---

## 6. Step 5 — QC exclusion flag

**Script:** `flag_qc_excluded.py`

Five genomes (2.3% of 220) had a severe gene-recovery failure in Panaroo
(≈0–45% of their annotated CDS were actually usable) despite passing every
structural check we could think of (matching `seqid`, correct sequence
length vs. CDS coordinates, no duplicate gene IDs). The root cause was not
identified; these five were excluded from pangenome/tree analyses and the
reason is recorded for transparency rather than silently dropped.

```bash
python3 flag_qc_excluded.py \
    --master master_table_final.tsv --out master_table_qc.tsv
```

> ⚠ **This is the one open issue in the whole pipeline.** If reproducing
> this project, it is worth revisiting these 5 samples with a fresh
> Panaroo/Prokka run to see if the problem is reproducible or was transient.

---

## 7. Step 6 — Reorganize into genus / functional groups

**Script:** `organize_by_group.py`

Copies each sample's `.faa`/`.gff` into two parallel hierarchies for
group-level analyses:

- `by_genus/<Genus>/` — only genera with ≥3 samples get their own folder;
  smaller ones are pooled into `_minor_genera/`.
- `by_functional_group/{LAB, Bacillus_group, Other_Environmental,
  Unresolved}/` — a coarse ecological/functional grouping (LAB = lactic
  acid bacteria genera; Bacillus_group = spore-forming Bacillaceae/
  Paenibacillaceae genera typical of *jang*-type fermentation).

```bash
python3 organize_by_group.py \
    --master master_table_qc.tsv --outdir grouped/ \
    --min-genus-n 3 --mode copy
```

---

## 8. Step 7 — Pangenome analysis (Panaroo)

**Script:** `run_panaroo_all_groups.sh` (env: `panaroo_env`)

```bash
conda activate panaroo_env
./run_panaroo_all_groups.sh grouped/ pangenome/ 24 3
```

Loops over every `by_genus/*/genomes_gff` and `by_functional_group/*/
genomes_gff` folder, skipping mixed-genus bins (`_minor_genera`,
`unresolved`, `Other_Environmental`) and any group already completed.

```
panaroo -i <group>/genomes_gff/*.gff -o <out> --clean-mode strict -a core -t 24
```

> ⚠ **Pitfall — `core_gene_alignment.aln` was never produced**
> Panaroo's alignment step is opt-in: without `-a core` (or `-a pan`), no
> alignment file is written at all, even though `gene_presence_absence.csv`
> looks complete. Always pass `-a core` if you plan to build a tree.

> ⚠ **Pitfall — `core_genes = 0` for taxonomically broad groups**
> This is *not a bug*. Panaroo/Roary define "core" as present in
> 99–100% of genomes using a fixed sequence-identity clustering threshold.
> Applied across an entire genus that actually contains several distinct
> species (e.g., `Bacillus`), or across a functional group spanning many
> genera (`LAB`, `Bacillus_group`), true orthologs are often too divergent
> in sequence to cluster together, and even one incomplete/divergent genome
> can zero out the whole "core" category. **Diagnosis approach:** don't
> trust the built-in category alone — recompute the full presence-frequency
> distribution directly from `gene_presence_absence.Rtab` (see §9) with
> flexible thresholds (100/99/95/90/80/50%) before concluding there is no
> shared genome.
> **Downstream consequence:** genera with 0 core genes also have no
> `core_gene_alignment.aln`, so no FastTree phylogeny can be built for them
> — an ANI-distance-based fallback tree is used instead (§13).

---

## 9. Step 7b — QC diagnostics and flexible core/accessory analysis

**Scripts:** `diagnose_pangenome_genomes.py`, `analyze_gene_presence.py`,
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

`diagnose_pangenome_genomes.py` is what identified the 5 problem genomes in
§6: it compares, per sample, the number of CDS actually annotated in the
input `.gff` against the number of gene families Panaroo recognized for
that genome. A large gap (ratio < 0.5 by default) flags a likely QC issue.

---

## 10. Step 8 — COG functional profile per genome

**Scripts:** `extract_eggnog_summary.py`

```bash
python3 extract_eggnog_summary.py \
    --master master_table_qc.tsv \
    --out-count eggnog_cog_count_wide.tsv \
    --out-ratio eggnog_cog_ratio_wide.tsv
```

Reads each sample's `annotation_EggNOG.xlsx` → `Eggnog_Count` sheet
(already aggregated by Macrogen: `Eggnog, Description, Count, Ratio (%)`)
and builds a sample × COG-category wide matrix. **Use the Ratio(%) table**
for cross-sample comparison (genome sizes differ).

> ⚠ **Pitfall — sheet name changed across sequencing years**
> 2021 data used the sheet name `Eggnog_Count`; 2023 data used
> `Eggnog Count` (space instead of underscore). **Fix:** normalize
> (lower-case, strip spaces/underscores) before matching sheet names,
> rather than requiring an exact string match.

> ⚠ **Pitfall — duplicate category rows within one sample's sheet crashed
> the batch merge**
> Some `Eggnog_Count` sheets had the same category code appear on more than
> one row. Naively building a `pandas.Series` from `set_index("Eggnog")`
> produced a non-unique index, which crashed when many such Series were
> combined into one DataFrame. **Fix:** `groupby("Eggnog").sum()` per
> sample before combining.

> ⚠ **Pitfall — a spurious "category" column**
> A separator line (`---------------------------------------------------
> ----`) in the source Excel was picked up as if it were a real COG code.
> All downstream comparison scripts explicitly exclude a small
> `NON_CATEGORY_COLS` set (`{"Total", "-", "----...----"}`).

---

## 11. Step 9 — Statistical comparison of functional profiles

**Scripts:** `compare_functional_groups.py` (LAB vs. Bacillus_group, 2-way),
`genus_cog_heatmap.py` (all genera, Kruskal-Wallis)

```bash
python3 compare_functional_groups.py \
    --ratio-tsv eggnog_cog_ratio_wide.tsv \
    --outdir functional_comparison/ \
    --group-a LAB --group-b Bacillus_group

python3 genus_cog_heatmap.py \
    --ratio-tsv eggnog_cog_ratio_wide.tsv \
    --outdir genus_cog_analysis/ --min-n 3
```

Per-category Mann-Whitney U test (2 groups) / Kruskal-Wallis (≥3 groups),
with Benjamini-Hochberg FDR correction, plus a grouped bar chart, a PCA of
all 220 genomes' COG profiles, and a genus × category heatmap.

> ⚠ **Pitfall — a single `NaN` p-value silently broke FDR correction for
> every category**
> A COG category with zero variance in both groups (e.g., category `Y`,
> which never occurs in bacteria) made `scipy.stats.mannwhitneyu` return
> `p = NaN`. The custom Benjamini-Hochberg implementation used
> `np.argsort`/`np.minimum.accumulate` over the raw p-value array; a single
> `NaN` corrupted the sort order for the *entire* array, so every corrected
> p-value came out non-significant even when raw p-values were as small as
> `1e-30`. **Fix:** replace `NaN` with `1.0` (treat as non-significant)
> before running the BH correction, and wrap the test call in
> `try/except ValueError` to catch fully-degenerate inputs.

---

## 12. Step 10 — Antibiotic resistance & virulence screening

**Scripts:** `run_abricate_batch.sh` (env: `abricate_env`),
`compare_resistance_virulence.py`

```bash
conda activate abricate_env
./run_abricate_batch.sh genomes_fna/ abricate_out/
# produces abricate_out/card_summary.tsv and vfdb_summary.tsv

conda activate compgenomics
python3 compare_resistance_virulence.py \
    --summary-tsv abricate_out/card_summary.tsv \
    --master master_table_qc.tsv --db-label CARD \
    --outdir resistance_comparison/ --group-a LAB --group-b Bacillus_group
# repeat with vfdb_summary.tsv / --db-label VFDB
```

Per-gene Fisher's exact test (presence/absence) + BH-FDR, plus a
Mann-Whitney U comparison of total gene "burden" (genes detected per
genome) between groups.

> ⚠ **Pitfall — pandas/matplotlib API changes broke the script mid-run**
> `DataFrame.applymap()` is removed in newer pandas (`Series`/`DataFrame`
> `.map()` should be used instead). `Axes.boxplot(..., labels=...)` was
> renamed to `tick_labels=` in a recent matplotlib. Both are one-line fixes,
> but worth checking your installed versions if copying this code as-is.

---

## 13. Step 11 — Secondary metabolite (BGC) screening with antiSMASH

**Scripts:** `merge_gbk_files.sh`, `run_antismash_batch.sh` (env:
`antismash_env`), `analyze_antismash_bgcs.py`

```bash
# 1) Merge per-contig GenBank files into one multi-record file per sample
#    (reuses the same source tree walk as consolidate_wgs.sh)
./merge_gbk_files.sh \
    "/mnt/f/WGS_Results/#Whole_Genome_Sequencing_Macrogen" \
    "/mnt/f/WGS_Consolidated"

# 2) Run antiSMASH per genome, reusing existing annotation
conda activate antismash_env
./run_antismash_batch.sh genomes_gbk/ antismash_out/ master_table_qc.tsv \
    24 "LAB,Bacillus_group"
./run_antismash_batch.sh genomes_gbk/ antismash_out/ master_table_qc.tsv \
    24 "Other_Environmental,Unresolved"   # remaining samples

# 3) Aggregate JSON output into a sample x BGC-type matrix + statistics
conda activate compgenomics
python3 analyze_antismash_bgcs.py \
    --antismash-root antismash_out/ --master master_table_qc.tsv \
    --outdir antismash_analysis/ --group-a LAB --group-b Bacillus_group
```

`--genefinding-tool none` reuses the existing Prokka/tbl2asn gene calls
instead of re-predicting genes, which is both faster and keeps gene
numbering consistent with the rest of the pipeline. Runtime was ~10–70 s
per genome (≈2–3 h total for 220 genomes on 24 threads).

The antiSMASH 7 JSON schema actually used:
`records[].areas[].products` (a list of strings per BGC region) —
confirmed by inspecting one sample's output before writing the batch
parser, rather than assuming a schema from memory or documentation.

> ⚠ **Pitfall — batch filtering script depended on pandas, which wasn't
> installed in `antismash_env`**
> The group-filter logic in `run_antismash_batch.sh` originally called out
> to a Python one-liner that imported pandas. `antismash_env` doesn't ship
> pandas by default. **Fix:** rewrote the filter purely in `awk` (no
> external dependency), which also runs faster for a simple TSV filter.

> ⚠ **Pitfall — UTF-8 BOM broke the awk column lookup, silently printing
> whole rows**
> `master_table_qc.tsv` was saved with `encoding="utf-8-sig"` (BOM), which
> attaches an invisible prefix to the very first header token
> (`sample_id`). awk's `$i == "sample_id"` comparison then never matched,
> leaving the column index unset (`0` in awk), and `print $0` (the *entire*
> line) was emitted instead of just the ID — this produced garbage-looking
† output (species names, file paths, etc., each treated as a separate
"sample" by the shell's word-splitting). **Fix:** hardcode `sample_id` as
column 1 by position (it always is, by construction of every table in this
pipeline) rather than relying on a text match for that one column.

> ⚠ **Pitfall (cosmetic, unfixed) — final summary counters showed 0/0
> despite all jobs succeeding**
> `printf '%s\n' "$LIST" | while read ...; do COUNTER=$((COUNTER+1)); done`
> runs the loop body in a subshell (because of the pipe), so increments to
> `COUNTER` are lost once the loop exits. All actual antiSMASH runs
> completed correctly (verified via per-sample logs); only the trailing
> "N new / M skipped" summary line was wrong. Left as a known cosmetic
> issue — use `find <out_dir> -name index.html | wc -l` to get an accurate
> completed-sample count instead of trusting the printed summary.

---

## 14. Step 12 — Integrative figures

### 14a. Per-genus tree + COG heatmap

**Script:** `genus_tree_heatmap.py`

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

Builds a FastTree phylogeny from `core_gene_alignment.aln` (`-nt` mode) and
draws it beside a per-sample COG heatmap. `--ani-matrix` is required as a
fallback for genera with zero core genes (here: `Bacillus`, `Lactococcus`),
in which case an ANI-distance hierarchical-clustering tree is drawn instead
and clearly labeled as such in the figure title.

### 14b. Whole-dataset integrative figure

**Script:** `global_tree_heatmap.py`

```bash
python3 global_tree_heatmap.py \
    --ani-matrix ani_analysis/ani_matrix.csv \
    --master master_table_qc.tsv \
    --antismash-summary antismash_analysis/antismash_bgc_summary.tsv \
    --card-summary abricate_out/card_summary.tsv \
    --vfdb-summary abricate_out/vfdb_summary.tsv \
    --outdir global_summary/
```

All 220 genomes, ANI-based dendrogram (tip labels colored by functional
group), with an adjacent normalized heatmap of: total antiSMASH BGC count,
CARD gene count, VFDB gene count. This single figure is what surfaced an
outlier sample (a *Klebsiella pneumoniae* contaminant with CARD=26,
VFDB=62, both far above any other sample) that would not have been obvious
from any single-analysis output alone.

---

## 15. Summary of Key Biological Findings

1. **LAB vs. Bacillus_group show a consistent, multi-layered functional
   divergence** across COG profile, CARD, and antiSMASH results: LAB
   genomes are enriched for carbohydrate metabolism / translation /
   replication-repair (fast, simple fermentation strategy); Bacillus_group
   genomes are enriched for motility, signal transduction, and a much
   larger secondary-metabolite arsenal (NRPS/PKS/siderophores; median 13
   vs. 4 BGCs per genome, p≈3×10⁻²⁸) plus a higher intrinsic
   antibiotic-resistance-gene burden.
2. A subset of Bacillus_group genomes carry the *B. cereus*-type
   non-hemolytic enterotoxin genes (`nheA/B/C`) — a food-safety-relevant
   finding worth manual follow-up.
3. One sample (`HN00179262_F4055`) was identified as *Klebsiella
   pneumoniae* via BLAST and stood out as an extreme outlier in the
   integrative CARD/VFDB heatmap — most plausibly a contamination event
   worth tracing back to its source material/process.

---

## 16. Script Inventory

| Stage | Script | Environment |
|---|---|---|
| Consolidation | `consolidate_wgs.sh`, `extract_species_from_blast.py` | compgenomics |
| Taxonomy | `normalize_taxonomy.py`, `resolve_species_via_ani.py`, `flag_qc_excluded.py` | compgenomics |
| ANI | `run_fastani.sh`, `analyze_ani.py` | compgenomics |
| Grouping | `organize_by_group.py` | compgenomics |
| Pangenome | `run_panaroo_all_groups.sh` | panaroo_env |
| Pangenome QC | `summarize_pangenomes.py`, `diagnose_pangenome_genomes.py`, `analyze_gene_presence.py` | compgenomics |
| COG function | `extract_eggnog_summary.py`, `compare_functional_groups.py`, `genus_cog_heatmap.py` | compgenomics |
| Resistance/virulence | `run_abricate_batch.sh`, `compare_resistance_virulence.py` | abricate_env / compgenomics |
| Secondary metabolites | `merge_gbk_files.sh`, `run_antismash_batch.sh`, `analyze_antismash_bgcs.py` | antismash_env / compgenomics |
| Integrative figures | `genus_tree_heatmap.py`, `global_tree_heatmap.py` | compgenomics |

---

## 17. General Lessons for Reproducing This Pipeline

- **Never assume a spreadsheet's header row or sheet name** — Macrogen's
  own output format changed between 2021 and 2023 batches (BLAST header
  row offset; eggNOG sheet-name spacing). Always inspect one real file
  before writing a batch parser.
- **Header/ID consistency between annotation (GFF/GBK) and sequence
  (FASTA) files is the single most common silent-failure point** for any
  downstream tool (Panaroo, antiSMASH). When merging multi-contig
  annotations, always pair each contig's annotation with *its own*
  sequence file rather than a separately-produced whole-genome FASTA.
- **Isolate conflicting tool dependencies into separate conda
  environments** rather than trying to force one shared environment
  (Panaroo needed Python ≤3.11; the rest of the pipeline auto-resolved to
  Python 3.14).
- **Don't trust a pangenome tool's built-in "core" definition at face
  value** when comparing across a taxonomically broad group — always
  recompute the full gene-presence-frequency distribution and check it
  with several thresholds before concluding "no shared genes."
- **Multiple-testing correction code must explicitly handle NaN p-values**
  (from zero-variance groups); otherwise a single degenerate test can
  silently invalidate every other result in the batch.
- **Save Excel-bound outputs with `encoding="utf-8-sig"`** and prefer
  comma-`.csv` over tab-`.tsv` for anything a collaborator will open by
  double-clicking in Excel — and be aware that a UTF-8 BOM can break naive
  `awk`/shell column-position logic on the header line unless the ID
  column's position is hardcoded rather than text-matched.
