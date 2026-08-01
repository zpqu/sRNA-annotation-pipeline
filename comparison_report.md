# 04_annotation — Comparison of Three Sense-Annotation Strategies

Date: 2026-08-01
Branch: `feature/union-containment-annotation`
Reference genome: mm39; samples: Cumulus-cells, Granulosa-cells

## 1. Strategies compared

| ID | Feature set | Sense overlap rule | Output dir |
|---|---|---|---|
| A | primary miRNA transcripts (miRBase pri-miRNA, 1227 feat.) | read fully contained in feature | `output/` |
| B | mature miRNA loci (miRBase v22, 2110 feat., 16–27 nt) | read fully contained in feature | `output_matmiRNA/` |
| C | mature miRNA loci (same as B) | **union containment**: read fully contained in feature OR feature fully contained in read | `output_union/` |

Antisense annotation is unchanged in all three (type="any", ignore.strand).

## 2. Composition and read totals (all reads, % of library)

| Category | A (primiRNA) | B (matmiRNA) | C (matmiRNA union) |
|---|---|---|---|
| miRNA-class | 85.6% / 82.3% | 82.0% / 79.3% | 85.3% / 82.1% |
| piRNA | 7.9% / 11.9% | 11.5% / 14.9% | 9.1% / 13.1% |
| RM | 2.0% / 2.2% | 2.0% / 2.2% | 1.3% / 1.5% |

(Cumulus / Granulosa.)

Union containment (C) returns the miRNA class to essentially the same
fraction as the pri-miRNA strategy (A): 85.3% vs 85.6% (Cumulus). The added
reads (2.26 M / 2.06 M) are drawn back from `piRNA` and `RM`.

## 3. miRNA-class read sizes (Cumulus; % of the miRNA category)

| Size | A (primiRNA) | B (matmiRNA) | C (matmiRNA union) |
|---|---|---|---|
| 20 nt | 66.6 | 69.4 | 66.7 |
| 21 nt | 20.4 | 21.1 | 20.3 |
| 22 nt | 3.8 | 3.8 | 3.7 |
| 23 nt | 3.0 | 0.2 | 3.0 |
| 24 nt | 0.8 | ~0 | 0.8 |
| **20–22 nt total** | **90.7** | **94.4** | **90.8** |

The union rule re-introduces 23–24 nt reads into the miRNA class, making C's
size profile almost identical to A's and diluting the canonical 20–22 nt
fraction from 94.4% (B) to 90.8% (C).

## 4. What the union rule adds (C vs B)

- Added reads: 2,263,403 (Cumulus, +4.1%) / 2,061,044 (Granulosa, +3.6%).
- Added-read sizes: 23 nt = 1.60 M, 24 nt = 0.48 M, 25 nt = 0.12 M, 26 nt
  = 23 k, plus a trace of 22/21 nt.
- **94.9% of the added reads were classified `primiRNA` in strategy A**
  (2.4% were `piRNA`/`RM`): they are reads that lie inside a pri-miRNA
  transcript and fully contain a mature locus but extend beyond its
  boundaries (a mature locus of 22 nt inside a 23–27 nt read). These are
  pri-miRNA processing intermediates / reads straddling the mature locus,
  not the final mature miRNA product (their sequence is not the mature
  miRNA, since it carries extra flanking nucleotides).
- Consequently C attributes ~93% of strategy A's "non-mature pri-miRNA"
  reads (4.2% of the pri-miRNA category) back to specific mature loci,
  which is why C ≈ A in total reads.

## 5. Sanity checks (all strategies)

- **Top identity**: all three agree the library is let-7-dominated
  (let-7a-5p, let-7c-5p, let-7f-5p, let-7b-5p, let-7e-5p); B and C rank
  mature miRNAs nearly identically (Spearman rho = 0.979 on the top 200).
- **5p/3p arm bias**: B 97.6:2.1, C 97.5:2.2 — unchanged, strongly
  5'-biased.
- **Strand specificity (C)**: matmiRNA sense/antisense = 7,362:1
  (Cumulus), 9,459:1 (Granulosa); all categories sense > antisense.
- **Top single read**: identical under all three strategies — chr13
  48,691,717–48,691,736, 7.30 M copies = mmu-let-7a-5p locus.

## 6. Verdict: which strategy is most biologically meaningful?

**B (matmiRNA-within) is the most biologically meaningful** for mature-miRNA
annotation, and the union-containment variant (C) is not an improvement:

1. **Product purity.** B confines every miRNA-class read to the interior of
   a mature miRNA locus, so 94.4% of its reads are the canonical 20–22 nt
   mature products (exact or truncated). C (like A) falls to ~90.8% by
   admitting 23–24 nt reads.
2. **Attribution precision.** The reads that C additionally counts merely
   *contain* a mature locus (94.9% were pri-miRNA fragments in A). Their
   sequence extends beyond the mature locus, so they are processing
   intermediates or boundary-straddling reads — attributing them to a
   specific mature miRNA is biologically imprecise and inflates counts of
   abundant loci (e.g. let-7a-5p +2.8%).
3. **What the union strategy is good for.** It is diagnostic, not
   quantitative: it shows that most of strategy A's "non-mature pri-miRNA"
   reads (4.2%) are reads straddling mature-locus boundaries (consistent
   with Drosha/Dicer processing intermediates), and that excluding them (as
   B does) is correct.

**Recommendation**: keep B (`output_matmiRNA/`, matmiRNA-within) as the
primary mature-miRNA annotation. C's added reads should not be counted as
mature miRNA expression.

## 7. Outputs

- `output_union/tables/`: Table3a (three-way read-annotation composition),
  Table3b (union-added matmiRNA reads), Table3c (mature-miRNA expression B vs C)
- `output_union/figures/`: Comparison1 (composition), Comparison2 (miRNA size)
- Scripts: `2c_annotation_smallRNA_unionContainment_RMlast.R` (step 2c),
  `7_compare_three_strategies.R` (step 7); `mygeneFeature.R` gained a
  `type = "within" | "union"` option (default unchanged).
