# 04_annotation — Pipeline Validation Report

Date: 2026-07-31
Scripts: `scripts/run_smallRNA_annotation.sh` (steps 1–5), `scripts/R/**`
Reference genome features: mm39 (`DB/`)

## 1. Datasets and library statistics (step 1, `Table1a/b/c/d`)

| Sample | Total reads | Unique reads | % singletons | Median count | Max count | Dominant size (unique) | Dominant size (reads) |
|---|---|---|---|---|---|---|---|
| Cumulus-cells | 67,905,011 | 2,180,678 | 79.4 | 1 | 7,304,281 | 17 nt | 20 nt |
| Granulosa-cells | 72,511,711 | 2,060,922 | 79.7 | 1 | 7,457,431 | 17 nt | 20 nt |

Read-size bins (all reads): 16–19 nt: 11.9%/15.1%, **20–24 nt: 87.2%/84.2%**, 25–33 nt: 0.9%/0.7%, 34–46 nt: <0.1% (Cumulus/Granulosa).

Both libraries are strongly miRNA-dominated by read count (20–22 nt peak), with a long tail of low-abundance 16–46 nt reads (piRNA/repeat fragments).

## 2. Annotation composition (step 2, `Table2a/b`, `read.annotation`)

All-reads % (unique-reads %):

| Category | Cumulus-cells | Granulosa-cells |
|---|---|---|
| primiRNA (miRNA loci) | **85.58%** (0.51%) | **82.33%** (0.44%) |
| piRNA | 7.90% (21.06%) | 11.90% (26.74%) |
| RM (repeat masked) | 1.97% (24.74%) | 2.21% (22.86%) |
| refGene.NM.intron | 0.83% (10.62%) | 0.70% (11.54%) |
| AS.piRNA | 1.29% (14.06%) | 0.88% (11.56%) |
| AS.RM | 0.76% (12.47%) | 0.45% (9.02%) |
| AS.refGene.NM.intron | 0.32% (4.10%) | 0.29% (4.46%) |
| other | 1.20% (10.91%) | 1.07% (11.07%) |
| tRNA | 0.01% (0.07%) | 0.01% (0.10%) |
| snoRNA | 0.01% (0.04%) | 0.03% (0.06%) |

Interpretation: a small set of ultra-abundant 20–22 nt miRNA sequences (hosted in pri-miRNA transcripts) accounts for ~83–86% of all reads in both samples, while the bulk of sequence *diversity* resides in piRNA loci, repeat-masked regions, and introns (sense ≫ antisense). Sense/antisense ratio (unique reads): 1.9–2.5.

Per-class sizes (all reads):
- **primiRNA (mature miRNA)**: 20 nt 66.6%, 21 nt 20.4%, 22 nt 3.8% → canonical miRNA profile.
- **piRNA-annotated**: 16–21 nt (20 nt 23.4%) — *not* canonical 26–32 nt (see caveat 1).
- **tRNA**: 24 nt 24.0%, 25 nt 22.8%, 26 nt 17.1% → tRNA-derived fragments.
- **snoRNA**: 18–21 nt → snoRNA-derived small RNAs (~20 nt).

## 3. Position profiles (step 5, `Fig2c`)

- **tRNA**: sense ≫ antisense (242k vs 1.5k reads, Cumulus; 235k vs 0.7k, Granulosa); reads cluster near the 5' end (peak window position 16–20), consistent with 5'-tRNA-derived fragments.
- **snoRNA**: sense ≫ antisense (234k vs 0.1k, Cumulus; 507k vs 0.1k, Granulosa), reads spread across the gene body.

## 4. Critical bug found and fixed

**Symptom:** only 24% of reads (16.3M of 67.9M) were annotated; the top-expressed miRNA reads were absent from all annotation tables.

**Root cause:** step 1 collapses reads by genomic position (`chr/start/end/strand` + sequence), so multi-mapping reads (97% of reads by count; up to 10 loci each) share the same read name. Step 2 and `mygeneFeature()` removed already-annotated reads **by name**, silently discarding all non-hitting positions of multi-mappers — including the most abundant reads (e.g., the 7.3M-copy dominant miRNA).

**Fix:** replace name-based removal with position-key removal (`chr:start:end:strand`) in:
- `scripts/R/02_annotation/2_annotation_smallRNA_RMlast.R`
- `scripts/R/02_annotation/mygeneFeature.R`

**After fix:** 100% of reads are accounted for in `read.annotation` and `read.size`; the top miRNA (7,304,281 copies) is annotated as `primiRNA`. Total accounted reads now equal the library totals exactly (67,905,011 / 72,511,711).

## 5. Caveats

1. **piRNA-annotated reads are 16–21 nt**, not canonical 26–32 nt. The `within` overlap with piRNA loci captures mostly other small-RNA species overlapping piRNA loci; treat the piRNA category as "reads inside piRNA loci" rather than verified piRNAs.
2. **Multi-mapper semantics:** reads are counted once per mapping position (position-based), so a sequence mapping to several loci is counted at each locus with its position-specific frequency.
3. **Minus-strand window indexing:** in `5_tRNA_snoRNA_position.widnows.R` positions are indexed from the 5' end on the plus strand but from the 3' end on the minus strand (inconsistent convention).
4. Size range is 16–46 nt by library design; very short (<16 nt) or long (>46 nt) small RNAs are not captured.

## 6. Output inventory

- `output/tables/` (26 files): Table1a–d (library statistics), Table2a/b (unique/all-reads counts), per-sample per-class `*.count.txt` and `*.size.count.txt`
- `output/figures/` (38 files): Figure1a/b (size vs count), Figure3a–e (annotation count/percentage, 2-D facet `sample ~ class`), Fig2a/b (per-class size barplots + percentages), Fig2c (tRNA/snoRNA position profiles + antisense)
- `output/rdata/` (10 files): `*.bam.gr.RData`, `*.bam.unique.gr.RData`, `*.bam.annotated.gr.RData` (carries `type`/`region`/`count`), `mm39.{tRNA,snoRNA}.20bp.gr.RData`, `{tRNA,snoRNA}.20bp.dis.all.df.RData`
- Logs: `scripts/logs/step1–5.log`, master run log (per pipeline launch)

## 7. Mature-miRNA strategy comparison (step 2b, `output_matmiRNA/`)

**Design.** A second annotation strategy (`2b_annotation_smallRNA_matmiRNA_RMlast.R`) assigns reads to *mature* miRNA loci (2,110 features, miRBase v22 `miRNA` entries, 16–27 nt) instead of *primary* miRNA transcripts (1,227 features, 46–147 nt). Same logic/priority as strategy A; only the feature set differs.

**Read totals.**

| Sample | primiRNA cat. (A) | matmiRNA cat. (B) | B/A | non-mature reads in pri-miRNA |
|---|---|---|---|---|
| Cumulus | 58.12 M (85.6%) | 55.69 M (82.0%) | 95.8% | 2.43 M (4.2%) |
| Granulosa | 59.70 M (82.3%) | 57.49 M (79.3%) | 96.3% | 2.20 M (3.7%) |

**Key findings.**

1. **matmiRNA captures ~96% of the primiRNA signal.** Only 3.7–4.2% of pri-miRNA reads lie within pri-miRNA transcripts but *outside* mature loci. These residual reads are enriched for 23–24 nt reads (primiRNA strategy: 23 nt = 3.0%/2.6%; matmiRNA: 0.2%/0.1%) — consistent with pri-miRNA hairpin/loop fragments and 5'-arm offset reads rather than processed miRNAs.
2. **Canonical mature sizes.** matmiRNA read sizes are 20 nt = 69%, 21 nt = 21%, 22 nt = 3.8% (both samples) — a textbook processed-miRNA size profile; max read 26 nt.
3. **Exact miRNA identity.** The matmiRNA strategy resolves the dominant signal to canonical, well-characterized mammalian miRNAs: let-7 family dominates (~50 M reads; let-7a-5p 20.3 M, let-7c-5p 12.3 M, let-7f-5p 11.4 M, let-7b-5p 6.1 M, let-7e-5p 2.3 M in Cumulus), followed by miR-103-3p, miR-140-3p, miR-322-3p, miR-503-5p, miR-21a-5p, miR-107-3p, miR-98-5p, miR-320-3p. let-7 dominance in somatic ovarian cells is well established — a strong external biological validation.
4. **Composition shift is small and specific.** Re-classified non-mature pri-miRNA reads (~3–4%) fall mostly into `piRNA` (+3.2 pp Cumulus, +3.0 pp Granulosa); all other categories are essentially unchanged. The extra reads therefore overlap piRNA-annotated positions rather than representing real miRNA signal.
5. **Tighter unique composition.** matmiRNA unique-read share (0.28%/0.24%) is lower than primiRNA (0.51%/0.44%): the mature-only set is confined to a smaller, more specific position set.

**Conclusion.** The matmiRNA strategy is the more biologically meaningful annotation for this small-RNA dataset: it (i) resolves the dominant signal to exact mature-miRNA identities with canonical 20–22 nt sizes, (ii) excludes pri-miRNA hairpin/loop fragments that inflate the primiRNA category, and (iii) retains ~96% of the primiRNA read total, so no real signal is lost. Figures: `output_matmiRNA/figures/Comparison1_read_annotation_primiRNA_vs_matmiRNA.pdf`, `Comparison2_miRNA_size_distribution.pdf`; tables: `Table2m_mature_miRNA_expression.csv`, `Table2n_strategy_comparison_read_annotation.csv`.
