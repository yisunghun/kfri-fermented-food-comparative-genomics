# Comparative Genomics of Bacterial Isolates from Korean Fermented Foods

Reproducible pipeline for comparative whole-genome analysis of ~220 bacterial
isolates from Korean fermented foods (kimchi, jang-type products), covering
species assignment (BLAST + ANI cross-validation), genus-level pangenome
analysis, COG functional profiling, antibiotic-resistance and virulence-gene
screening, secondary-metabolite biosynthetic gene cluster prediction, and
integrative tree/heatmap visualization.

This repository accompanies:

> Yi, S. Comparative Genomics of 220 Bacterial Isolates from Korean Fermented
> Foods Reveals Divergent Functional Strategies Between Lactic Acid Bacteria
> and Bacillus-Group Genera. *International Journal of Food Microbiology*
> (submitted).

For the full step-by-step protocol, environment setup, and a catalogue of
issues encountered during development (with fixes), see
[`docs/PROTOCOL.md`](docs/PROTOCOL.md) (English) and
[`docs/PROTOCOL_KR.md`](docs/PROTOCOL_KR.md) (Korean).

## Repository structure

```
scripts/
  01_consolidation/          Merge per-contig annotations into per-genome files; BLAST-based species calls
  02_taxonomy/                Lactobacillus nomenclature normalization; ANI-based species resolution; QC exclusion
  03_ani/                      fastANI all-vs-all + dendrogram/species cross-validation
  04_grouping/                 Reorganize genomes into genus / functional-group hierarchies
  05_pangenome/                Panaroo batch runner + pangenome QC and flexible presence-frequency analysis
  06_functional_cog/          eggNOG/COG extraction and group-level statistical comparison
  07_resistance_virulence/    CARD/VFDB screening (abricate) and group comparison
  08_secondary_metabolites/   GenBank merging, antiSMASH batch runner, BGC summary/comparison
  09_integrative_figures/     Per-genus and whole-dataset tree + heatmap figures

environments/                 Conda environment definitions (see below)
docs/                         Full written protocol (English + Korean)
```

## Environment setup

Four separate conda environments are used due to conflicting Python-version
requirements between tools:

```bash
mamba env create -f environments/compgenomics.yml     # fastANI, pandas, scipy, FastTree, Biopython, etc.
mamba env create -f environments/panaroo_env.yml       # Panaroo (needs Python <=3.11)
mamba env create -f environments/antismash_env.yml     # antiSMASH
mamba env create -f environments/abricate_env.yml      # abricate (CARD/VFDB screening)
```

`antiSMASH` additionally requires downloading its reference databases once:

```bash
conda activate antismash_env
download-antismash-databases
```

`abricate` requires a one-time database setup:

```bash
conda activate abricate_env
abricate --setupdb
```

## Usage

Each script is runnable standalone and documents its own arguments via
`--help` (Python scripts) or a usage comment block at the top (shell
scripts). Scripts are numbered by the order in which they were used in the
pipeline; see `docs/PROTOCOL.md` for the exact command-line invocations,
expected inputs/outputs, and known pitfalls for each step.

## Data availability

Raw sequence data and full genome assemblies are available from the
corresponding author upon reasonable request. Complete, closed genome
assemblies for a subset of nine isolates are publicly available in NCBI
GenBank (accessions listed in the manuscript, Table 4).

## Citation

If you use this pipeline, please cite the associated manuscript (see above)
and the underlying tools listed in the manuscript's Methods section and
References (fastANI, Panaroo, FastTree, eggNOG-mapper, CARD, VFDB, antiSMASH,
ABRicate, Prokka, BLAST+, pandas, SciPy, scikit-learn, Matplotlib).

## License

Code in this repository is released under the MIT License (see `LICENSE`).

## Contact

Sunghun Yi, Korea Food Research Institute \u2014 sunghunyi@kfri.re.kr
