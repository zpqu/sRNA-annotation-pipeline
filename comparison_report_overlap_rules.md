# 04_annotation — Comparison of Three Sense-Overlap Rules on Mature miRNA Loci

Date: 2026-08-01
Branch: `feature/overlap-rule-comparison`
Reference genome: mm39; samples: Cumulus-cells, Granulosa-cells

## 1. Strategies compared

All three strategies use the **same feature set** — mature miRNA loci
(miRBase v22, 2110 features, 16–27 nt) — and differ **only** in the sense
overlap rule:

| ID | Feature set | Sense overlap rule | Output dir |
|---|---|---|---|
| A | mature miRNA loci | **any overlap** (≥1 bp: read overlaps locus, ignore.strand for antisense) | `output_any/` |
| B | mature miRNA loci | **read fully contained** in feature (within) | `output_matmiRNA/` |
| C | mature miRNA loci | **union containment** (read in feature OR feature in read) | `output_union/` |

Antisense annotation is identical in all three (type="any", ignore.strand).
Priority order identical in all three: matmiRNA > snoRNA > piRNA > tRNA >
gene > RM > other.

## 2. Composition — all reads, % of library (Cumulus / Granulosa)

| Category | A (any) | B (within) | C (union) |
|---|---|---|---|
| matmiRNA | 85.6 / 82.3 | 82.0 / 79.3 | 85.3 / 82.1 |
| piRNA | 10.1 / 13.9 | 11.5 / 14.9 | 9.1 / 13.1 |
| RM | 1.1 / 1.2 | 2.0 / 2.2 | 1.3 / 1.5 |
| AS.piRNA | 0.36 / 0.24 | 1.29 / 0.89 | 1.05 / 0.71 |
| snoRNA | 0.027 / 0.040 | 0.015 / 0.029 | 0.015 / 0.029 |
| tRNA | 0.001 / 0.002 | 0.009 / 0.011 | 0.006 / 0.008 |

Both relaxed rules (A any, C union) inflate the miRNA class back to ~85%,
i.e. close to the pri-miRNA-within strategy (85.6%), by absorbing reads from
`piRNA`, `RM`, `tRNA` and the antisense categories.

## 3. Category totals — Cumulus (dAB = A − B, dCB = C − B, percentage points)

| Category | A (any) | B (within) | C (union) | dAB | dCB |
|---|---|---|---|---|---|
| matmiRNA | 85.59 | 82.00 | 85.34 | **+3.59** | **+3.33** |
| piRNA | 10.09 | 11.47 | 9.06 | −1.38 | −2.41 |
| RM | 1.12 | 1.97 | 1.34 | −0.85 | −0.63 |
| AS.piRNA | 0.36 | 1.29 | 1.05 | −0.93 | −0.24 |

Read-count totals: matmiRNA A = 58.12 M, B = 55.69 M, C = 57.95 M.

## 4. Read movement (cross-tabulation, B → A and B → C)

**B → A (within → any), top transitions (both samples, reads):**
- piRNA → matmiRNA: 4,625,446  ← reads that only partially overlap a mature locus
- RM → piRNA: 1,506,670        ← cascading: piRNA-any captures reads that were RM in B
- snoRNA → snoRNA: 31,295 (stable core; A adds +18,569→ reads in Cumulus, incl. 27 nt)
- tRNA → piRNA: 11,799 of 12,604 (A's piRNA-any rule absorbs 94% of B's tRNA reads)

**B → C (within → union):** piRNA → matmiRNA 4,317,012; RM → piRNA 978,431;
tRNA signal also reduced (6,330 → 4,147 in Cumulus, −34%).

## 5. miRNA-class read sizes (Cumulus, % of matmiRNA category)

| Size | A (any) | B (within) | C (union) |
|---|---|---|---|
| 20 nt | 66.6 | 69.4 | 66.7 |
| 21 nt | 20.4 | 21.1 | 20.3 |
| 22 nt | 3.8 | 3.8 | 3.7 |
| 23 nt | 3.0 | 0.2 | 3.0 |
| 24 nt | 0.8 | ~0 | 0.8 |
| **20–22 nt total** | **90.7** | **94.4** | **90.8** |
| reads | 58.12 M | 55.69 M | 57.95 M |

The +3.6 pp of matmiRNA in A/C consist of 23–24 nt reads at mature-locus
boundaries (pri-miRNA processing intermediates / boundary-straddling reads).

## 6. Effect on other genome-feature categories (Cumulus)

| Category | B (within) | A (any) | change | C (union) | change |
|---|---|---|---|---|---|
| snoRNA reads | 10,011 | 18,569 | **+85%** | 9,991 | ~0% |
| tRNA reads | 6,330 | 642 | **−90%** | 4,147 | **−34%** |
| piRNA reads | 7.79 M | 6.85 M | −12% | 6.15 M | −21% |

The "any" rule cascades through the whole annotation: the piRNA-any rule
absorbs tRNA-derived fragments before they reach the tRNA category (tRNA
signal nearly destroyed), while snoRNA-any adds boundary/27-nt reads that
double the snoRNA class. Union (C) similarly erodes tRNA (−34%) and piRNA.

## 7. Mature miRNA quantification across the three rules

Top identity (let-7a-5p, Cumulus):

| Strategy | reads | unique genomic positions |
|---|---|---|
| A (any) | 20,901,780 | 199 |
| B (within) | 20,323,453 | **56** |
| C (union) | 20,888,922 | 89 |

- Spearman rank correlation vs B (top-200 expressed loci): A 0.969, C 0.979
  — ranking is robust to the overlap rule.
- A assigns reads to let-7a-5p from 199 distinct positions (3.6× more than B):
  reads that touch a locus edge by ≥1 bp are counted as that miRNA.
- 5p/3p arm ratios are nearly identical: A 97.4:2.3, B 97.6:2.1, C 97.5:2.2
  (Cumulus); 97.9:1.8 / 98.1:1.7 / 98.0:1.7 (Granulosa).
- Strand specificity of A is healthy: matmiRNA sense:antisense 51,754:1
  (Cumulus) / 59,638:1 (Granulosa).

## 8. Which strategy is most biologically meaningful?

**B (read-fully-within mature loci) is the most biologically meaningful.**

1. **Canonical size purity.** B's matmiRNA category is 94.4% canonical
   20–22 nt reads vs 90.7% (A) and 90.8% (C). The extra reads in A/C are
   23–24 nt processing intermediates that straddle locus boundaries — not
   mature miRNA products.
2. **Locus-level specificity.** A assigns reads from 3.6× more genomic
   positions to each miRNA (199 vs 56 for let-7a-5p): a read touching a
   locus by a single base is counted as that miRNA under "any".
3. **Integrity of other features.** "Any" is applied to every category and
   cascades: it destroys tRNA annotation (−90%), doubles snoRNA with
   boundary reads (+85%), shrinks piRNA (−12%) and halves RM. Union also
   erodes tRNA (−34%) and piRNA (−21%).
4. **Ranking is unaffected** (rho ≥ 0.969), so B's stricter rule loses no
   biological signal for differential expression while keeping absolute
   counts and per-category annotation clean.

Both A (any) and C (union) mainly re-add reads that the pri-miRNA-within
strategy would also call primiRNA — i.e. the relaxed rules blur the
boundary between mature miRNA products and their precursors. The strict
containment rule (B) is the cleanest definition of a mature miRNA.
