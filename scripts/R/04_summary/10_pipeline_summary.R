##Date: 2026
##Author: Zhipeng
## Step 10: end-of-pipeline summary. After all annotation/analysis steps have
## finished, this script
##   1) lists every output file (with a description and size), grouped by
##      output directory / producing step,
##   2) summarises the main results (sample stats, read sizes, annotation
##      composition, top mature miRNAs, per-category abundance),
##   3) runs a set of sanity checks (totals consistent, tables consistent,
##      logs error-free) and reports PASS / FAIL / WARN per check.
## It auto-detects the optional analyses: if output_union/, output_any/ or
## step-9 abundance tables (Table5a) are present (e.g. from a manual run of
## steps 2c/2d/7/8/9), they are described and checked too.
##
## Run from scripts/R/04_summary/ (same convention as the other pipeline
## scripts) or via run_smallRNA_annotation.sh.
##
## Outputs (step 10):
##   output_matmiRNA/pipeline_summary.md        this report
##   output_matmiRNA/tables/Table6_sanity_checks.csv

suppressMessages({library(data.table); library(GenomicRanges)})

source("../../../config/genome.R")

## -----------------------------------------------------------------------------
## small helpers
## -----------------------------------------------------------------------------
out = function(...) cat(sprintf(...), "\n")
rep.lines = character()
rep = function(...) rep.lines <<- c(rep.lines, sprintf(...))

read.tab = function(f){
  if (file.exists(f)) fread(f) else NULL
}
fmt = function(x) format(x, big.mark = ",", scientific = FALSE)
pct2 = function(x) round(100 * x, 2)

rep("## Small-RNA annotation pipeline - end-of-pipeline summary\n")
rep("- Reference genome : `%s`", genome)
rep("- Feature DB       : `%s` (built by scripts/R/00_build_DB/)", db.dir)
rep("- Samples          : %s", paste(samples, collapse = ", "))
rep("- Generated        : %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
rep("- Strategy dirs found: %s",
    paste(sprintf("`%s`", c("output", "output_matmiRNA",
                            if (dir.exists(file.path(project.root, "output_union"))) "output_union",
                            if (dir.exists(file.path(project.root, "output_any"))) "output_any")), collapse = ", "))
rep("- Optional step 9 (abundance) : %s", ifelse(file.exists(file.path(out.dir, "tables/Table5a_category_abundance_summary.csv")), "run", "not run"))
rep("")

## -----------------------------------------------------------------------------
## 1. Output-file inventory
## -----------------------------------------------------------------------------
rep("## 1. Output files\n")
rep("Each output file produced by the pipeline is listed below with a short description and its size on disk.\n")

SAM <- paste(sprintf("(%s)", samples), collapse = "|")
desc.rules <- list(
  ## step 1 tables
  list(re = "^Table1a_sample_summary\\.csv$", d = "Per-sample summary: total reads, unique reads, singletons, dominant read sizes"),
  list(re = "^Table1b_read_size_distribution\\.csv$", d = "Read size distribution (n_unique and n_reads per read size)"),
  list(re = "^Table1c_read_count_distribution\\.csv$", d = "Read count (expression) distribution, log2 bins"),
  list(re = "^Table1d_read_size_vs_count\\.csv$", d = "2D read-size x read-count grid"),
  ## step 2 / 2b consolidated count tables
  list(re = "^Table2a_annotation_count_unique_reads\\.csv$", d = "Annotation composition, unique reads (each unique read counted once)"),
  list(re = "^Table2b_annotation_count_all_reads\\.csv$", d = "Annotation composition, all reads (weighted by read frequency)"),
  ## step 2b unique-read annotation tables
  list(re = "^Table_unique_reads_annotation\\.csv$", d = "Per-unique-read annotation + abundance, all samples combined"),
  list(re = sprintf("^(%s)\\.unique_reads_annotation\\.csv$", SAM), d = "Per-unique-read annotation + abundance for this sample (sorted by count descending)"),
  ## step 2 / 2b per-sample count files
  list(re = sprintf("^(%s)\\.read\\.annotation\\.count\\.txt$", SAM), d = "Read category x all-read counts"),
  list(re = sprintf("^(%s)\\.read\\.size\\.count\\.txt$", SAM), d = "Read size x all-read counts"),
  list(re = sprintf("^(%s)\\.(primiRNA|matmiRNA|snoRNA|piRNA|tRNA)\\.annotation\\.count\\.txt$", SAM), d = "Category reads x genomic-context counts (CDS/5UTR/3UTR/intron/up1k/down1k/RM/intergenic)"),
  list(re = sprintf("^(%s)\\.(primiRNA|matmiRNA|snoRNA|piRNA|tRNA)\\.size\\.count\\.txt$", SAM), d = "Category reads x read-size counts"),
  ## step 6 (miRNA strategy comparison)
  list(re = "^Table2m_mature_miRNA_expression\\.csv$", d = "Top expressed mature miRNAs (name, n_unique, n_reads) per sample"),
  list(re = "^Table2n_strategy_comparison_read_annotation\\.csv$", d = "primiRNA vs matmiRNA annotation-composition comparison"),
  ## step 7 (three strategies A/B/C)
  list(re = "^Table3a_strategy_comparison_read_annotation\\.csv$", d = "A/B/C read-annotation composition"),
  list(re = "^Table3b_union_added_matmiRNA_reads\\.csv$", d = "Reads added to the matmiRNA category by the union rule, by size"),
  list(re = "^Table3c_mature_miRNA_expression_B_vs_C\\.csv$", d = "Top mature miRNA expression under strategy B vs C"),
  ## step 8 (three overlap rules)
  list(re = "^Table4a_overlap_rule_composition\\.csv$", d = "A/B/C (any/within/union) overlap-rule composition"),
  list(re = "^Table4b_overlap_rule_category_totals\\.csv$", d = "A/B/C per-category read totals"),
  list(re = "^Table4c_read_movement_B_vs_A\\.csv$", d = "Read movement between within (B) and any (A) strategies"),
  list(re = "^Table4d_read_movement_B_vs_C\\.csv$", d = "Read movement between within (B) and union (C) strategies"),
  ## step 9 (abundance)
  list(re = "^Table5a_category_abundance_summary\\.csv$", d = "Per-category abundance statistics (n loci, total reads, Gini, top-1/5/10 %, top locus, cross-sample tests)"),
  list(re = "^Table5b_per_locus_abundance\\.csv$", d = "Per-locus read abundance, long format"),
  ## step 10 (this script)
  list(re = "^Table6_sanity_checks\\.csv$", d = "Pipeline sanity-check results (step 10)"),
  ## rdata
  list(re = sprintf("^(%s)\\.bam\\.gr\\.RData$", SAM), d = "All aligned reads (GRanges)"),
  list(re = sprintf("^(%s)\\.bam\\.unique\\.gr\\.RData$", SAM), d = "Non-redundant unique reads with 'count' (expression) column"),
  list(re = sprintf("^(%s)\\.bam\\.annotated\\.gr\\.RData$", SAM), d = "Annotated reads (category 'type', gene context 'region', 'count', 'feature.id', 'n_features')"),
  list(re = "^(snoRNA|tRNA)\\.20bp\\.gr\\.RData$", d = "20 bp sliding windows along each gene body (step 5)"),
  list(re = "^(snoRNA|tRNA)\\.20bp\\.dis\\.all\\.df\\.RData$", d = "Weighted read counts per window position (step 5)"),
  ## figures
  list(re = "^Figure1[ab]\\.read_size_vs_count\\.(pdf|png)$", d = "Read size / count distribution figures (step 1)"),
  list(re = "^Figure3[a-e]\\..*\\.(pdf|png)$", d = "Annotation count / percentage barplots (step 3)"),
  list(re = "^Fig2a\\..*_size_barplot\\.pdf$", d = "Per-class read-size barplots (step 4)"),
  list(re = "^Fig2b\\..*_size_barplot\\.percentage\\.pdf$", d = "Per-class read-size percentage barplots (step 4)"),
  list(re = "^Fig2c\\.(tRNA|snoRNA)_pos_barplot(_AS)?\\.pdf$", d = "Position-distribution barplots for tRNA/snoRNA genes (step 5)"),
  list(re = "^Comparison1_read_annotation_primiRNA_vs_matmiRNA\\.pdf$", d = "primiRNA vs matmiRNA read-annotation composition (step 6)"),
  list(re = "^Comparison2_miRNA_size_distribution\\.pdf$", d = "miRNA read-size distribution, primiRNA vs matmiRNA (step 6)"),
  list(re = "^Comparison1_three_strategies_read_annotation\\.pdf$", d = "A/B/C read-annotation composition (step 7)"),
  list(re = "^Comparison2_three_strategies_miRNA_size\\.pdf$", d = "A/B/C miRNA read-size distributions (step 7)"),
  list(re = "^Comparison1_overlap_rules_composition\\.pdf$", d = "A/B/C overlap-rule composition (step 8)"),
  list(re = "^Comparison2_overlap_rules_category_size\\.pdf$", d = "A/B/C per-category read-size distributions (step 8)"),
  list(re = "^Abundance1_per_locus_distribution\\.pdf$", d = "Per-locus log10 abundance distribution per category (step 9)"),
  list(re = "^Abundance2_rank_abundance\\.pdf$", d = "Rank-abundance (Whittaker) curves per category (step 9)"),
  list(re = "^Abundance3_lorenz_matmiRNA\\.pdf$", d = "Lorenz curves for mature-miRNA loci (step 9)")
)
describe.file = function(f){
  for (r in desc.rules) if (grepl(r$re, f)) return(r$d)
  "See producing step for details"
}

dir.map <- list(
  "output/tables"         = "Steps 1-2: read-preprocessing and annotation count tables (primary-miRNA strategy)",
  "output/rdata"          = "Steps 1-2: read and annotated-read GRanges (primary-miRNA strategy)",
  "output/figures"        = "Steps 1, 3-5: read-size, annotation and position-distribution figures",
  "output_matmiRNA/tables"= "Steps 2b, 6, 9, 10: mature-miRNA strategy tables, comparisons, abundance and this summary",
  "output_matmiRNA/rdata" = "Step 2b: annotated-read GRanges (mature-miRNA strategy)",
  "output_matmiRNA/figures" = "Steps 6, 9: strategy-comparison and abundance figures",
  "output_union/tables"   = "Steps 2c, 7: union-containment strategy tables (optional)",
  "output_union/rdata"    = "Step 2c: annotated-read GRanges (optional)",
  "output_union/figures"  = "Step 7: three-strategy comparison figures (optional)",
  "output_any/tables"     = "Steps 2d, 8: any-overlap strategy tables (optional)",
  "output_any/rdata"      = "Step 2d: annotated-read GRanges (optional)",
  "output_any/figures"    = "Step 8: overlap-rule comparison figures (optional)"
)
n.files = 0L
for (d in names(dir.map)){
  full = file.path(project.root, d)
  if (!dir.exists(full)) { rep("### `%s`\n_(not run in this analysis)_\n", d); next }
  files = list.files(full)
  if (length(files) == 0) { rep("### `%s`\n_(empty)_\n", d); next }
  rep("### `%s`\n\n_%s_\n", d, dir.map[[d]])
  for (f in files){
    sz = format(file.info(file.path(full, f))$size, big.mark = ",")
    rep("- `%s` (%s B) - %s", f, sz, describe.file(f))
    n.files = n.files + 1L
  }
  rep("")
}

rep("### Ancillary inputs\n")
rep("- `config/genome.R` - single source of truth for the reference genome, samples and output directories")
rep("- `DB/rdata_%s/` - feature GRanges DB (%s objects), built by step 0", genome, length(list.files(db.dir, pattern = "\\.RData$")))
rep("- `scripts/logs/step*.log` - console log of each step")
rep("")

## -----------------------------------------------------------------------------
## 2. Results summary
## -----------------------------------------------------------------------------
rep("## 2. Results summary\n")

## ---- 2.1 sample-level -------------------------------------------------------
t1a = read.tab(file.path(out.dir.pri, "tables/Table1a_sample_summary.csv"))
if (!is.null(t1a)){
  rep("### 2.1 Sample-level read statistics (Table1a)\n")
  rep("| sample | total reads | unique reads | singletons | %% singletons | median count | max count | dominant size by reads |")
  rep("|---|---|---|---|---|---|---|---|")
  for (i in seq_len(nrow(t1a)))
    rep("| %s | %s | %s | %s | %.1f | %s | %s | %d nt |",
        t1a$sample[i], fmt(t1a$total_reads[i]), fmt(t1a$unique_reads[i]),
        fmt(t1a$singleton_reads[i]), t1a$pct_singletons[i],
        fmt(t1a$median_count[i]), fmt(t1a$max_count[i]), t1a$dominant_size_by_reads_nt[i])
  rep("")
}

## ---- 2.2 read-size composition ----------------------------------------------
t1b = read.tab(file.path(out.dir.pri, "tables/Table1b_read_size_distribution.csv"))
if (!is.null(t1b)){
  rep("### 2.2 Read-size composition (Table1b)\n")
  rep("| sample | reads 20-22 nt | %% of all reads | dominant size |")
  rep("|---|---|---|---|")
  for (s in samples){
    x = t1b[sample == s]
    if (nrow(x) == 0) next
    f = pct2(sum(x$n_reads[x$width %in% 20:22]) / sum(x$n_reads))
    rep("| %s | %s | %.2f | %d nt |", s, fmt(sum(x$n_reads[x$width %in% 20:22])), f,
        x$width[which.max(x$n_reads)])
  }
  rep("")
}

## ---- 2.3 annotation composition (all reads) ---------------------------------
comp.md = function(t2, tag){
  x = t2[category == "read.annotation"]
  if (nrow(x) == 0) return(invisible())
  tot = x[, .(tot = sum(Freq)), by = sample]
  x = merge(x, tot, by = "sample")
  x[, pct := round(100 * Freq / tot, 2)]
  setorder(x, sample, -Freq)
  rep("### %s\n", tag)
  rep("| sample | category | reads | %% of reads |")
  rep("|---|---|---|---|")
  for (i in seq_len(nrow(x)))
    rep("| %s | %s | %s | %.2f |", x$sample[i], x$item[i], fmt(x$Freq[i]), x$pct[i])
  rep("")
}

t2b = read.tab(file.path(out.dir, "tables/Table2b_annotation_count_all_reads.csv"))
if (!is.null(t2b)) comp.md(t2b, "2.3a Annotation composition, all reads (mature-miRNA strategy, output_matmiRNA/Table2b)")
t2a = read.tab(file.path(out.dir, "tables/Table2a_annotation_count_unique_reads.csv"))
if (!is.null(t2a)) comp.md(t2a, "2.3b Annotation composition, unique reads (mature-miRNA strategy, output_matmiRNA/Table2a)")
t2b.pri = read.tab(file.path(out.dir.pri, "tables/Table2b_annotation_count_all_reads.csv"))
if (!is.null(t2b.pri)) comp.md(t2b.pri, "2.3c Annotation composition, all reads (primary-miRNA strategy, output/Table2b)")

## ---- 2.4 mature-miRNA read sizes (20-22 nt) ---------------------------------
rep("### 2.4 Mature-miRNA read sizes (20-22 nt; step 2b size tables)\n")
rep("| sample | matmiRNA reads | matmiRNA reads 20-22 nt | %% |")
rep("|---|---|---|---|")
for (s in samples){
  mt = read.tab(file.path(out.dir, "tables", paste0(s, ".matmiRNA.size.count.txt")))
  if (is.null(mt)) next
  mt[, w := as.integer(Var1)]
  f = pct2(sum(mt$Freq[mt$w %in% 20:22]) / sum(mt$Freq))
  rep("| %s | %s | %s | %.2f |", s, fmt(sum(mt$Freq)),
      fmt(sum(mt$Freq[mt$w %in% 20:22])), f)
}
rep("")

## ---- 2.5 top mature miRNAs ---------------------------------------------------
t2m = read.tab(file.path(out.dir, "tables/Table2m_mature_miRNA_expression.csv"))
if (!is.null(t2m)){
  rep("### 2.5 Top expressed mature miRNAs (Table2m)\n")
  for (s in samples){
    x = t2m[sample == s][order(-n_reads)][1:min(10, .N)]
    rep("**%s**\n", s)
    rep("| rank | mature miRNA | n_unique | n_reads |")
    rep("|---|---|---|---|")
    for (i in seq_len(nrow(x)))
      rep("| %d | %s | %s | %s |", i, x$name[i], fmt(x$n_unique[i]), fmt(x$n_reads[i]))
    rep("")
  }
}

## ---- 2.6 per-category abundance (step 9, if run) ------------------------------
t5a = read.tab(file.path(out.dir, "tables/Table5a_category_abundance_summary.csv"))
if (!is.null(t5a)){
  rep("### 2.6 Per-category abundance (Table5a, step 9)\n")
  rep("| category | sample | n_loci | total reads | Gini | top1 %% | top5 %% | top10 %% | top locus |")
  rep("|---|---|---|---|---|---|---|---|---|")
  for (i in seq_len(nrow(t5a)))
    rep("| %s | %s | %s | %s | %.3f | %.1f | %.1f | %.1f | %s |",
        t5a$category[i], t5a$sample[i], fmt(t5a$n_loci[i]), fmt(t5a$total_reads[i]),
        t5a$gini[i], t5a$top1_pct[i], t5a$top5_pct[i], t5a$top10_pct[i], t5a$top_locus[i])
  rep("")
}

## ---- 2.7 per-unique-read annotation table ------------------------------------
rep("### 2.7 Per-unique-read annotation + abundance table (step 2b)\n")
rep("| sample | unique reads | total reads | categories | top feature (reads) |")
rep("|---|---|---|---|---|")
for (s in samples){
  csv = file.path(out.dir, "tables", paste0(s, ".unique_reads_annotation.csv"))
  if (!file.exists(csv)) next
  dt = fread(csv)
  top = dt[which.max(count)]
  rep("| %s | %s | %s | %d | %s | %s |",
      s, fmt(nrow(dt)), fmt(sum(dt$count)), uniqueN(dt$category),
      ifelse(is.na(top$feature_id), "other", top$feature_id), fmt(top$count))
}
rep("")

## -----------------------------------------------------------------------------
## 3. Sanity checks
## -----------------------------------------------------------------------------
rep("## 3. Sanity checks\n")
rep("Each check returns PASS, FAIL or WARN (WARN = input not available, e.g. optional step not run).\n")

checks = data.table(id = character(), description = character(), status = character(), detail = character())
check = function(id, desc, ok, detail = ""){
  st = if (is.na(ok)) "WARN" else if (ok) "PASS" else "FAIL"
  checks <<- rbind(checks, data.table(id = id, description = desc, status = st, detail = as.character(detail)))
}

total.reads = if (!is.null(t1a)) setNames(as.numeric(t1a$total_reads), t1a$sample) else NULL
uniq.reads  = if (!is.null(t1a)) setNames(as.numeric(t1a$unique_reads), t1a$sample) else NULL

for (s in samples){
  ## totals
  if (!is.null(t2b)){
    v = t2b[category == "read.annotation" & sample == s, sum(Freq)]
    check(paste0("s1_", s), sprintf("sum of all category counts equals total reads (%s)", s),
          !is.null(total.reads) && isTRUE(all.equal(v, as.numeric(total.reads[s]))),
          sprintf("observed %s, expected %s", fmt(v), fmt(total.reads[s])))
  }
  if (!is.null(t2a)){
    v = t2a[category == "read.annotation" & sample == s, sum(Freq)]
    check(paste0("s2_", s), sprintf("sum of unique-read category counts equals unique reads (%s)", s),
          !is.null(uniq.reads) && isTRUE(all.equal(v, as.numeric(uniq.reads[s]))),
          sprintf("observed %s, expected %s", fmt(v), fmt(uniq.reads[s])))
  }
  ## unique-reads annotation CSV
  csv = file.path(out.dir, "tables", paste0(s, ".unique_reads_annotation.csv"))
  if (file.exists(csv)){
    dt = fread(csv)
    check(paste0("s3_", s), sprintf("unique-reads table row count equals unique reads (%s)", s),
          isTRUE(all.equal(nrow(dt), as.numeric(uniq.reads[s]))),
          sprintf("observed %s, expected %s", fmt(nrow(dt)), fmt(uniq.reads[s])))
    check(paste0("s4_", s), sprintf("unique-reads table sum(count) equals total reads (%s)", s),
          isTRUE(all.equal(sum(dt$count), as.numeric(total.reads[s]))),
          sprintf("observed %s, expected %s", fmt(sum(dt$count)), fmt(total.reads[s])))
    if (!is.null(t2b)){
      v = t2b[category == "read.annotation" & item == "matmiRNA" & sample == s, Freq]
      check(paste0("s5_", s), sprintf("Table2b matmiRNA equals sum of matmiRNA reads in unique-reads table (%s)", s),
            isTRUE(all.equal(v, sum(dt$count[dt$category == "matmiRNA"]))),
            sprintf("observed %s, expected %s", fmt(sum(dt$count[dt$category == "matmiRNA"])), fmt(v)))
    }
    check(paste0("s6_", s), sprintf("unique-reads table sorted by count descending, all counts >= 1 (%s)", s),
          dt$count[1] == max(dt$count) && all(dt$count >= 1), "")
    check(paste0("s7_", s), sprintf("sum(cpm) ~ 1e6 in unique-reads table (%s)", s),
          isTRUE(abs(sum(dt$cpm) - 1e6) / 1e6 < 1e-3),
          sprintf("observed %s", fmt(round(sum(dt$cpm))))
    )
  }
  ## read sizes
  if (!is.null(t1b)){
    x = t1b[sample == s]
    ok = isTRUE(all.equal(sum(x$n_reads), as.numeric(total.reads[s]))) &&
         isTRUE(all.equal(sum(x$n_unique), as.numeric(uniq.reads[s])))
    check(paste0("s8_", s), sprintf("Table1b read sizes sum to total/unique reads (%s)", s), ok, "")
  }
  ## matmiRNA 20-22 nt fraction
  mt = read.tab(file.path(out.dir, "tables", paste0(s, ".matmiRNA.size.count.txt")))
  if (!is.null(mt)){
    mt[, w := as.integer(Var1)]
    f = sum(mt$Freq[mt$w %in% 20:22]) / sum(mt$Freq)
    check(paste0("s9_", s), sprintf("matmiRNA reads are predominantly 20-22 nt (>= 0.90) (%s)", s),
          f >= 0.90, sprintf("observed %.4f", f))
  }
  ## n_unique <= n_reads per item
  if (!is.null(t2a) && !is.null(t2b)){
    a = t2a[category == "read.annotation" & sample == s]
    b = t2b[category == "read.annotation" & sample == s]
    m = merge(a, b, by = "item", suffixes = c(".u", ".r"))
    ok = all(m$Freq.u <= m$Freq.r)
    check(paste0("s10_", s), sprintf("unique-read count <= all-read count for every category (%s)", s), ok, "")
  }
  ## primary-miRNA strategy totals
  if (!is.null(t2b.pri)){
    v = t2b.pri[category == "read.annotation" & sample == s, sum(Freq)]
    check(paste0("s11_", s), sprintf("primary-miRNA strategy: category counts sum to total reads (%s)", s),
          isTRUE(all.equal(v, as.numeric(total.reads[s]))),
          sprintf("observed %s, expected %s", fmt(v), fmt(total.reads[s])))
  }
}

## optional strategies: union and any
for (tag in c("output_union", "output_any")){
  f = file.path(project.root, tag, "tables/Table2b_annotation_count_all_reads.csv")
  if (file.exists(f)){
    tb = fread(f)
    for (s in samples){
      v = tb[category == "read.annotation" & sample == s, sum(Freq)]
      check(paste0("s12_", gsub("output_", "", tag), "_", s),
            sprintf("%s strategy: category counts sum to total reads (%s)", tag, s),
            isTRUE(all.equal(v, as.numeric(total.reads[s]))),
            sprintf("observed %s, expected %s", fmt(v), fmt(total.reads[s])))
    }
  }
}

## Table5a vs Table2b consistency
if (!is.null(t5a) && !is.null(t2b)){
  for (s in samples){
    v2 = t2b[category == "read.annotation" & item == "matmiRNA" & sample == s, Freq]
    v5 = t5a[category == "matmiRNA" & sample == s, total_reads]
    if (length(v5) == 1)
      check(paste0("s13_", s), sprintf("Table5a matmiRNA total equals Table2b matmiRNA (%s)", s),
            isTRUE(all.equal(as.numeric(v5), as.numeric(v2))),
            sprintf("observed %s, expected %s", fmt(v5), fmt(v2)))
  }
}

## feature DB
db.files = list.files(db.dir, pattern = "\\.RData$")
check("s14_db", "feature DB has the expected 24 objects", length(db.files) >= 24,
      sprintf("observed %d objects", length(db.files)))
if (file.exists(file.path(db.dir, "matmiRNA.gr.RData"))){
  e = new.env(); load(file.path(db.dir, "matmiRNA.gr.RData"), envir = e)
  n.mir = length(e$matmiRNA.gr)
  if (genome == "mm39")
    check("s15_db", "mature miRNA loci in DB == 2110 (mm39)", n.mir == 2110, sprintf("observed %d", n.mir))
  else
    check("s15_db", sprintf("mature miRNA loci in DB (reported for %s)", genome), NA, sprintf("observed %d", n.mir))
}

## logs
log.dir = file.path(project.root, "scripts", "logs")
expected.logs = c("step1.log", "step2.log", "step2b.log", "step3.log", "step4.log", "step5.log", "step6.log")
missing.logs = expected.logs[!file.exists(file.path(log.dir, expected.logs))]
if (length(missing.logs) > 0){
  check("s16_logs", "all step logs present", FALSE, paste("missing:", paste(missing.logs, collapse = ", ")))
}else{
  check("s16_logs", "all step logs present", TRUE, "")
  bad = sapply(expected.logs, function(l){
    txt = readLines(file.path(log.dir, l), warn = FALSE)
    any(grepl("^Error|Execution halted", txt))
  })
  check("s17_logs", "no fatal errors in step logs", !any(bad),
        ifelse(any(bad), paste("errors in:", paste(names(bad)[bad], collapse = ", ")), ""))
}

## ---- sanity-check summary table ---------------------------------------------
rep("### Sanity-check results\n")
rep("| id | description | status | detail |")
rep("|---|---|---|---|")
for (i in seq_len(nrow(checks)))
  rep("| %s | %s | **%s** | %s |", checks$id[i], checks$description[i], checks$status[i],
      ifelse(checks$detail[i] == "", "-", checks$detail[i]))
rep("")
rep("**Overall**: %d PASS, %d FAIL, %d WARN",
    sum(checks$status == "PASS"), sum(checks$status == "FAIL"), sum(checks$status == "WARN"))

## -----------------------------------------------------------------------------
## write outputs
## -----------------------------------------------------------------------------
report.file = file.path(out.dir, "pipeline_summary.md")
writeLines(rep.lines, report.file)
fwrite(checks, file.path(out.dir, "tables/Table6_sanity_checks.csv"))

cat(sprintf("Pipeline summary written to %s\n", report.file))
cat(sprintf("Sanity checks: %d PASS, %d FAIL, %d WARN; %d files described\n",
            sum(checks$status == "PASS"), sum(checks$status == "FAIL"),
            sum(checks$status == "WARN"), n.files))
