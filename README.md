# Small-RNA annotation pipeline

An R/bash pipeline for annotating small-RNA sequencing reads against a reference
genome. Each read is assigned to a single feature class — mature miRNA (miRBase),
snoRNA, piRNA, tRNA, repetitive element, or genic context (exon/intron/UTR) — in
either sense or antisense orientation, and the results are summarised as count
matrices, per-read tables, publication-style figures and an end-of-pipeline
report with built-in sanity checks.

The pipeline is genome-agnostic: the reference assembly and the annotation
strategy are defined in one configuration file, and any assembly can be added by
building a feature database from six raw feature files.

## Features

- Mature-miRNA annotation against miRBase loci (default feature set).
- Parameterised sense-annotation rule (see [Strategies](#strategies)) with a
  built-in comparison mode that runs all three rules and reports how they differ.
- Automatic feature-DB detection: the runner rebuilds the database from raw
  files if it is missing or incomplete.
- Automatic sample detection from BAM file names; per-sample and combined
  outputs.
- End-of-pipeline report (`pipeline_summary.md`) with automated sanity checks on
  every table (verifies internal consistency against the raw read counts).

## Repository layout

```
.
├── config/genome.R          # single source of truth: genome, strategy, samples
├── scripts/
│   ├── run_smallRNA_annotation.sh   # pipeline runner (the entry point)
│   ├── R/
│   │   ├── 00_build_DB/     # step 00: build the feature DB from raw files
│   │   ├── 01_preprocess/   # step 01: format BAMs -> reads GRanges + counts
│   │   ├── 02_annotation/   # step 02: annotate reads (main analysis)
│   │   ├── 03_figures/      # steps 03-06, s01: figures + comparison analysis
│   │   └── 04_summary/      # step 10: pipeline summary + sanity checks
├── bams/                    # input: aligned BAM files (gitignored)
├── DB/                      # feature DB (built, gitignored)
└── output/                  # results (gitignored)
```

## Requirements

- R (>= 4.0) with the following packages:

  ```r
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install(c("GenomicRanges", "GenomicFeatures", "rtracklayer",
                         "Rsamtools", "GenomicAlignments"))
  install.packages(c("data.table", "ggplot2", "scales", "patchwork"))
  ```

- A Unix-like shell with `bash`.

## Inputs

1. **Aligned BAM files** in `bams/` (e.g. `bams/Cumulus-cells.bowtie2.bam`).
   Sample labels are derived from the file names.
2. **Raw feature files** (only needed the first time a genome is used) in
   `DB/original_data_<genome>/`:

   | File          | Source        | Description                          |
   |---------------|---------------|--------------------------------------|
   | `refGene.gtf` | UCSC          | protein-coding (NM_) + noncoding (NR_) transcripts |
   | `miRNA.gff3`  | miRBase v22   | mature miRNA loci                    |
   | `tRNAs.bed`   | tRNAscan-SE   | tRNA genes                           |
   | `RM.bed`      | RepeatMasker  | repeat elements                      |
   | `piRBase.bed` | piRBase       | piRNA loci                           |
   | `piRNAdb.gtf` | piRNAdb       | piRNA clusters                       |

   The runner detects whether `DB/rdata_<genome>/` already holds the 24 feature
   objects and builds it via step 00 automatically if not.

## Usage

```bash
scripts/run_smallRNA_annotation.sh [GENOME] [STRATEGY]
scripts/run_smallRNA_annotation.sh -h | --help   # print usage and exit
```

| Argument | Description |
|---|---|
| `GENOME`   | Reference assembly, default `mm39` (also `SMALLRNA_GENOME`). `hg38` is supported; other assemblies work once a feature DB is built. |
| `STRATEGY` | Sense-annotation rule, default `fully-contained` (also `SMALLRNA_STRATEGY`). |

Environment variables:

| Variable | Description |
|---|---|
| `SMALLRNA_GENOME` | Same as the `GENOME` argument. |
| `SMALLRNA_STRATEGY` | Same as the `STRATEGY` argument. |
| `SMALLRNA_FORCE_REBUILD_DB=1` | Force a feature-DB rebuild even if the existing DB is complete. |
| `SMALLRNA_SUBSTRATEGY` | Set by the runner inside comparison mode — do not set manually. |

Examples:

```bash
scripts/run_smallRNA_annotation.sh                 # mm39, fully-contained
scripts/run_smallRNA_annotation.sh hg38 union      # human, union rule
scripts/run_smallRNA_annotation.sh hg38 comparison # all three rules + comparison analysis
```

## Strategies

The sense-annotation rule defines which read/feature overlaps count as a sense
hit. Antisense annotation is always "any overlap" (`ignore.strand = TRUE`).

| Strategy | Rule | Output dir |
|---|---|---|
| `fully-contained` (default) | read is fully contained in the feature | `output/` |
| `union` | read within feature **or** feature within read | `output/` |
| `any` | any overlap (≥ 1 bp) | `output/` |
| `comparison` | runs all three rules, then the overlap-rule comparison analysis | `output/comparison/` |

For a single-strategy run the outputs are written to `output/`; a comparison run
writes shared step-1 results to `output/comparison/` and per-strategy subfolders
(`fully_contained/`, `union/`, `any/`) beneath it.

## Pipeline steps

| Step | Script | What it does |
|---|---|---|
| 00 | `00_build_DB/00_build_annotation_DB.R` | Builds the 24 feature GRanges objects from the six raw files into `DB/rdata_<genome>/`. |
| 01 | `01_preprocess/01_format_bam.R` | Formats BAMs into per-read GRanges, computes unique reads and counts (`Table_01a–d`, `Figure_01`). |
| 02 | `02_annotation/02_annotation_smallRNA.R` | Annotates each read to a single feature class, sense/antisense; count matrices + per-read CSV (`Table_02a/b`, `Table_02m`, `Table_02_<sample>_unique_reads_annotation.csv`). |
| 03 | `03_figures/03_abundance_by_category.R` | Abundance/per-locus analyses for mature miRNAs (`Table_03a/b`, `Figure_03a–c`). |
| 04 | `03_figures/04_figure_annotation.R` | Annotation-composition barplots (`Figure_04a–e`). |
| 05 | `03_figures/05_figure_size.R` | Read-size distributions (`Figure_05a/b`). |
| 06 | `03_figures/06_tRNA_snoRNA_position.windows.R` | Positional profiles over tRNA/snoRNA genes (`Figure_06`). |
| s01 | `03_figures/s01_compare_overlap_rules.R` | (Comparison mode only) quantifies how the three rules differ (`Table_s01a–e`, `Figure_s01a/b`). |
| 10 | `04_summary/10_pipeline_summary.R` | Writes `pipeline_summary.md` and runs automated sanity checks (`Table_10`). |

## Outputs

Everything is written below `output/` (or `output/comparison/`), organised into
`tables/`, `figures/` and `rdata/`:

- `tables/` — cross-sample count matrices (`Table_01a–d`, `Table_02a/b`,
  `Table_02m`, `Table_03a/b`, `Table_s01a–e`, `Table_10`), per-sample
  annotation tables (`Table_02_<sample>_unique_reads_annotation.csv`).
- `figures/` — PDF/PNG figures (`Figure_01`, `Figure_03a–c`, `Figure_04a–e`,
  `Figure_05a/b`, `Figure_06`, `Figure_s01a/b`).
- `rdata/` — intermediate R objects (reads GRanges, annotation results,
  position-window counts) for downstream analysis.
- `pipeline_summary.md` — end-of-pipeline report: run parameters, inventory of
  every output file and the results of the sanity checks.

## Configuration

All run-level settings live in `config/genome.R`:

- `genome` — assembly ID (`mm39` default; overridable via `SMALLRNA_GENOME`).
- `strategy` — default annotation rule (overridable via `SMALLRNA_STRATEGY`).
- `samples` — sample labels in display order for the figure scripts.
- `chr.style` — chromosome naming convention (`"chr"`), enforced on the DB build.

## License

This project is released under the MIT License. See [LICENSE](LICENSE).
