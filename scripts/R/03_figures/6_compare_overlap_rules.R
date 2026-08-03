##Date: 2026
##Author: Zhipeng
## Step 6 (comparison mode only): three-way comparison of SENSE OVERLAP RULES on
## the SAME mature-miRNA (matmiRNA) feature set. The three strategies differ ONLY
## in the sense overlap rule:
##   fully-contained : read fully contained in feature     -> <out.base>/fully_contained/
##   union            : read within feature OR feature within read -> <out.base>/union/
##   any              : any overlap (>= 1 bp)               -> <out.base>/any/
## Contents (covers ALL genome-feature categories, not only miRNA):
##   1) full read.annotation composition (unique + all reads) across strategies
##   2) per-category read-size distributions (matmiRNA, snoRNA, piRNA, tRNA)
##   3) per-category read totals and how reads move between categories
##   4) mature-miRNA-focused: top identities, 5p/3p, sizes, union/any-added reads
##   5) strand specificity for the "any" strategy
## Outputs are written to <out.base>/tables and <out.base>/figures, where out.base
## is output/comparison (config/genome.R).

library(ggplot2)
library(scales)
library(data.table)
library(GenomicRanges)

out = function(...) cat(sprintf(...), "\n")

source("../../../config/genome.R")
stopifnot(is.comparison)

dir.tab = file.path(out.base, "tables")
dir.fig = file.path(out.base, "figures")
dir.create(dir.tab, recursive = TRUE, showWarnings = FALSE)
dir.create(dir.fig, recursive = TRUE, showWarnings = FALSE)

strategy.list = c("fully-contained" = "fully_contained",
                  "union" = "union",
                  "any" = "any")
cols = c("fully-contained" = "#33a02c", "union" = "#1f78b4", "any" = "#e31a1c")
for(st in strategy.list){
  if(!dir.exists(file.path(out.base, st)))
    stop("strategy dir missing: ", file.path(out.base, st), " (run steps 2/3/4/5/9 for each strategy first)")
}

strategy.dir.of = function(st) file.path(out.base, strategy.list[[st]])

## ---- 1) full read.annotation composition ------------------------------------------
d.all = data.table()
for(st in names(strategy.list)){
  dir = strategy.dir.of(st)
  t2a = read.csv(file.path(dir, "tables", "Table2a_annotation_count_unique_reads.csv"))
  t2b = read.csv(file.path(dir, "tables", "Table2b_annotation_count_all_reads.csv"))
  ra.u = as.data.table(t2a)[category == "read.annotation", .(sample, item, n_unique = Freq)]
  ra.a = as.data.table(t2b)[category == "read.annotation", .(sample, item, n_reads = Freq)]
  m = merge(ra.u, ra.a, by = c("sample", "item"))
  m[, pct_unique := 100*n_unique/base::sum(n_unique), by = sample]
  m[, pct_reads := 100*n_reads/base::sum(n_reads), by = sample]
  m[, strategy := st]
  d.all = rbind(d.all, m)
}
d.all[, item := gsub("^AS\\.", "AS_", item)]
d.all[, item := gsub("refGene\\.NM\\.", "gene.", item)]
d.all[, strategy := factor(strategy, levels = names(strategy.list))]
d.all[, item := factor(item, levels = c("matmiRNA", "snoRNA", "piRNA", "tRNA", "RM",
                                        "gene.exon", "gene.intron", "lincRNA.exon",
                                        "AS_matmiRNA", "AS_snoRNA", "AS_piRNA", "AS_tRNA",
                                        "AS_RM", "AS_gene.exon", "AS_gene.intron",
                                        "AS_lincRNA.exon", "other"))]
write.csv(d.all[order(strategy, -pct_reads)], file.path(dir.tab, "Table4a_overlap_rule_composition.csv"), row.names = FALSE)

p.comp = ggplot(d.all, aes(x = item, y = pct_reads, fill = strategy)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  facet_grid(sample ~ ., scales = "free_y") +
  labs(title = "read.annotation composition - all reads (overlap-rule strategies)",
       y = "% of reads", x = NULL) +
  scale_fill_manual(values = cols) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(p.comp, file = file.path(dir.fig, "Comparison1_overlap_rules_composition.pdf"),
       width = 13, height = 6)

## ---- 2) per-category read-size distributions ----------------------------------------
size.all = data.table()
for(st in names(strategy.list)){
  dir = strategy.dir.of(st)
  for(cat in c("matmiRNA", "snoRNA", "piRNA", "tRNA")){
    f = list.files(file.path(dir, "tables"), pattern = paste0("cells.", cat, ".size.count.txt$"))
    for(ff in f){
      x = read.delim(file.path(dir, "tables", ff))
      x$sample = sub(paste0(".", cat, ".size.count.txt"), "", ff)
      x$cat = cat
      x$strategy = st
      size.all = rbind(size.all, as.data.table(x))
    }
  }
}
size.all[, pct := as.numeric(100*Freq/base::sum(Freq)), by = .(sample, cat, strategy)]
size.all[, strategy := factor(strategy, levels = names(strategy.list))]

p.size = ggplot(size.all, aes(x = Var1, y = pct, fill = strategy)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8) +
  facet_grid(cat ~ sample, scales = "free") +
  labs(title = "per-category read-size distribution (all reads)",
       x = "read length (nt)", y = "% within category") +
  scale_fill_manual(values = cols) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(p.size, file = file.path(dir.fig, "Comparison2_overlap_rules_category_size.pdf"),
       width = 12, height = 10)

## ---- 3) per-category totals + read-movement between strategies ----------------------
cat.tot = dcast(d.all, sample + item ~ strategy, value.var = "pct_reads", fill = 0)
cat.tot[, dany := `any` - `fully-contained`]
cat.tot[, dunion := union - `fully-contained`]
write.csv(cat.tot[order(sample, -`fully-contained`)],
          file.path(dir.tab, "Table4b_overlap_rule_category_totals.csv"), row.names = FALSE)

## read-movement cross-tabulation (per read position: type in any vs type in fully-contained)
mov.tab = data.table()
for(s in samples){
  gr.list = list()
  for(st in names(strategy.list)){
    load(file.path(strategy.dir.of(st), "rdata", paste0(s, ".bam.annotated.gr.RData")))
    gr.list[[st]] = data.table(poskey = paste(seqnames(reads.bam.annotated.gr), start(reads.bam.annotated.gr),
                                              end(reads.bam.annotated.gr), strand(reads.bam.annotated.gr)),
                               type = reads.bam.annotated.gr$type,
                               count = as.numeric(reads.bam.annotated.gr$count))
  }
  m = merge(gr.list$any, gr.list$`fully-contained`[, .(poskey, typeB = type)],
            by = "poskey", all = TRUE)
  m = merge(m, gr.list$union[, .(poskey, typeC = type)], by = "poskey", all = TRUE)
  m[is.na(typeB), typeB := "other"]
  m[is.na(typeC), typeC := "other"]
  setnames(m, "type", "typeA")
  mov.tab = rbind(mov.tab, data.table(sample = s, m[, .(typeA, typeB, typeC, count)]))
}
## aggregate read movement (all reads) between fully-contained and any
mov.BA = mov.tab[, .(n_reads = base::sum(count)), by = .(typeB, typeA)]
setorder(mov.BA, -n_reads)
write.csv(mov.BA, file.path(dir.tab, "Table4c_read_movement_contained_vs_any.csv"), row.names = FALSE)
## between fully-contained and union
mov.BC = mov.tab[, .(n_reads = base::sum(count)), by = .(typeB, typeC)]
setorder(mov.BC, -n_reads)
write.csv(mov.BC, file.path(dir.tab, "Table4d_read_movement_contained_vs_union.csv"), row.names = FALSE)

out("== read movement fully-contained -> any [top transitions] ==")
for(k in 1:min(10, nrow(mov.BA)))
  out("  %-14s -> %-14s %12.0f", mov.BA$typeB[k], mov.BA$typeA[k], mov.BA$n_reads[k])
out("== read movement fully-contained -> union [top transitions] ==")
for(k in 1:min(10, nrow(mov.BC)))
  out("  %-14s -> %-14s %12.0f", mov.BC$typeB[k], mov.BC$typeC[k], mov.BC$n_reads[k])

## ---- 4) mature-miRNA-focused -------------------------------------------------------
load(file.path(db.dir, "matmiRNA.gr.RData"))
matmiRNA.gr = matmiRNA.gr[!duplicated(paste(seqnames(matmiRNA.gr),
                                                      start(matmiRNA.gr),
                                                      end(matmiRNA.gr),
                                                      strand(matmiRNA.gr)))]
expr.tab = data.table()
for(st in names(strategy.list)){
  dir = strategy.dir.of(st)
  files = list.files(file.path(dir, "rdata"), pattern = ".bam.annotated.gr.RData$")
  for(f in files){
    load(file.path(dir, "rdata", f))
    s = gsub("\\.bam\\.annotated\\.gr\\.RData", "", f)
    g = reads.bam.annotated.gr[reads.bam.annotated.gr$type == "matmiRNA"]
    if(length(g) == 0) next
    ol = findOverlaps(g, matmiRNA.gr, type = "any")
    hit = data.table(sample = s, strategy = st,
                     name = matmiRNA.gr$Name[subjectHits(ol)],
                     rid = queryHits(ol),
                     count = as.numeric(g$count[queryHits(ol)]))
    hit = hit[, .(n_unique = length(unique(rid)), n_reads = base::sum(count)),
              by = .(sample, strategy, name)]
    expr.tab = rbind(expr.tab, hit)
  }
}
setorder(expr.tab, sample, strategy, -n_reads)
write.csv(expr.tab, file.path(dir.tab, "Table4e_mature_miRNA_expression_strategies.csv"), row.names = FALSE)
for(s in unique(expr.tab$sample)){
  for(st in names(strategy.list)){
    top = expr.tab[sample == s & strategy == st][1:10]
    out("-- %s | %s top mature miRNAs (reads) --", s, st)
    for(k in 1:nrow(top)) out("   %-22s %10.0f", top$name[k], top$n_reads[k])
  }
}
## 5p/3p per strategy
for(s in unique(expr.tab$sample)){
  for(st in names(strategy.list)){
    sub = expr.tab[sample == s & strategy == st]
    arm = ifelse(grepl("-3p$", sub$name), "3p", ifelse(grepl("-5p$", sub$name), "5p", "other"))
    at = data.table(arm, n = sub$n_reads)[, .(n = base::sum(n)), by = arm]
    at[, pct := round(100*n/base::sum(n), 1)]
    out("%s %s 5p/3p: %s", s, st, paste(sprintf("%s=%.1f%%", at$arm, at$pct), collapse = ", "))
  }
}

## ---- 5) strand specificity for the "any" strategy -------------------------------------------
for(s in samples){
  ann = read.delim(file.path(strategy.dir.of("any"), "tables", paste0(s, ".read.annotation.count.txt")))
  dt = as.data.table(ann)
  dt[, type2 := sub("^AS\\.", "", Var1)]
  dt[, str2 := ifelse(grepl("^AS\\.", Var1), "AS", "S")]
  w = dcast(dt, type2 ~ str2, value.var = "Freq", fun.aggregate = base::sum, fill = 0)
  w[, ratio := round(S/pmax(AS, 1), 1)]
  setorder(w, -S)
  out("-- %s strand specificity (strategy: any, top 6) --", s)
  for(k in 1:min(6, nrow(w))) out("   %-18s sense=%.0f antisense=%.0f (ratio %.0f)", w$type2[k], w$S[k], w$AS[k], w$ratio[k])
}

out("Overlap-rule comparison written to %s and %s", dir.tab, dir.fig)
