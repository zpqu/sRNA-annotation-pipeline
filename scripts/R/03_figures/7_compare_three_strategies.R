##Date: 2026
##Author: Zhipeng
## Step 7 (optional): three-way comparison of sense-annotation strategies.
##   A: primiRNA-within (primary miRNA transcripts, read fully contained in feature)  -> output/
##   B: matmiRNA-within (mature miRNA loci, read fully contained in feature)           -> output_matmiRNA/
##   C: matmiRNA-union   (mature miRNA loci, read within feature OR feature within read)-> output_union/
## Contents:
##   1) read.annotation composition (unique vs all reads) across A/B/C, both samples
##   2) miRNA-category read-size distribution across A/B/C
##   3) reads added to the matmiRNA category by the union rule and their read sizes
##   4) top expressed mature miRNAs under B vs C
##   5) strand-specificity check for C
## Outputs are written to output_union/figures and output_union/tables.

library(ggplot2)
library(scales)
library(data.table)
library(GenomicRanges)

out = function(...) cat(sprintf(...), "\n")

## ---- load the three strategies --------------------------------------------------
t2b.list = list(
  primiRNA = read.csv("../../../output/tables/Table2b_annotation_count_all_reads.csv"),
  matmiRNA = read.csv("../../../output_matmiRNA/tables/Table2b_annotation_count_all_reads.csv"),
  matmiRNA_union = read.csv("../../../output_union/tables/Table2b_annotation_count_all_reads.csv")
)
t2a.list = list(
  primiRNA = read.csv("../../../output/tables/Table2a_annotation_count_unique_reads.csv"),
  matmiRNA = read.csv("../../../output_matmiRNA/tables/Table2a_annotation_count_unique_reads.csv"),
  matmiRNA_union = read.csv("../../../output_union/tables/Table2a_annotation_count_unique_reads.csv")
)

## item -> simpler label for the composition comparison
rename.item = function(x){
  x = gsub("^AS\\.", "AS_", x)
  x = gsub("refGene\\.NM\\.", "gene.", x)
  x
}

mk.ra = function(t2a, t2b, strategy){
  ra.u = as.data.table(t2a)[category == "read.annotation"]
  ra.a = as.data.table(t2b)[category == "read.annotation"]
  d.u = ra.u[, .(sample, item, n_unique = Freq)][, pct_unique := 100*n_unique/base::sum(n_unique), by = sample]
  d.a = ra.a[, .(sample, item, n_reads = Freq)][, pct_reads := 100*n_reads/base::sum(n_reads), by = sample]
  d = merge(d.u, d.a, by = c("sample", "item"))
  d[, strategy := strategy]
  d
}
d.all = rbind(mk.ra(t2a.list$primiRNA, t2b.list$primiRNA, "A_primiRNA_within"),
              mk.ra(t2a.list$matmiRNA, t2b.list$matmiRNA, "B_matmiRNA_within"),
              mk.ra(t2a.list$matmiRNA_union, t2b.list$matmiRNA_union, "C_matmiRNA_union"))
d.all[, item := rename.item(item)]
d.sort = d.all[order(-pct_reads)]
write.csv(d.sort, "../../../output_union/tables/Table3a_strategy_comparison_read_annotation.csv", row.names = FALSE)

## ---- 1) composition figure -------------------------------------------------------
lv = c("matmiRNA", "primiRNA", "snoRNA", "piRNA", "tRNA", "RM", "gene.exon",
       "gene.intron", "lincRNA.exon", "AS_matmiRNA", "AS_primiRNA", "AS_snoRNA",
       "AS_piRNA", "AS_tRNA", "AS_RM", "AS_gene.exon", "AS_gene.intron",
       "AS_lincRNA.exon", "other")
d.all[, item := factor(item, levels = lv[!is.na(match(lv, unique(d.all$item)))])]
d.all[, strategy := factor(strategy, levels = c("A_primiRNA_within", "B_matmiRNA_within", "C_matmiRNA_union"))]

p.comp = ggplot(d.all, aes(x = item, y = pct_reads, fill = strategy)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  facet_grid(sample ~ ., scales = "free_y") +
  labs(title = "read.annotation composition - all reads (3 sense-annotation strategies)",
       y = "% of reads", x = NULL) +
  scale_fill_manual(values = c(A_primiRNA_within = "#1f78b4", B_matmiRNA_within = "#33a02c",
                               C_matmiRNA_union = "#e31a1c")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(p.comp, file = "../../../output_union/figures/Comparison1_three_strategies_read_annotation.pdf",
       width = 13, height = 6)

## ---- 2) miRNA-category read-size distribution -------------------------------------
size.all = data.frame()
for(s in unique(t2b.list$matmiRNA$sample)){
  f.a = paste0("../../../output/tables/", s, ".primiRNA.size.count.txt")
  f.b = paste0("../../../output_matmiRNA/tables/", s, ".matmiRNA.size.count.txt")
  f.c = paste0("../../../output_union/tables/", s, ".matmiRNA.size.count.txt")
  a = read.delim(f.a); a$sample = s; a$strategy = "A_primiRNA_within"
  b = read.delim(f.b); b$sample = s; b$strategy = "B_matmiRNA_within"
  c = read.delim(f.c); c$sample = s; c$strategy = "C_matmiRNA_union"
  size.all = rbind(size.all, a, b, c)
}
size.all = as.data.table(size.all)
size.all[, pct := as.numeric(100*Freq/base::sum(Freq)), by = .(sample, strategy)]
size.all[, strategy := factor(strategy, levels = c("A_primiRNA_within", "B_matmiRNA_within", "C_matmiRNA_union"))]

p.size = ggplot(size.all, aes(x = Var1, y = pct, fill = strategy)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8) +
  facet_wrap(~sample, ncol = 1) +
  labs(title = "miRNA-category read-size distribution (all reads)",
       x = "read length (nt)", y = "% of reads in the miRNA class") +
  scale_fill_manual(values = c(A_primiRNA_within = "#1f78b4", B_matmiRNA_within = "#33a02c",
                               C_matmiRNA_union = "#e31a1c")) +
  theme_bw()
ggsave(p.size, file = "../../../output_union/figures/Comparison2_three_strategies_miRNA_size.pdf",
       width = 9, height = 7)

## ---- 3) reads added to the matmiRNA category by the union rule ---------------------
cat.tab = data.table()
for(s in unique(t2b.list$matmiRNA$sample)){
  m.b = read.delim(paste0("../../../output_matmiRNA/tables/", s, ".matmiRNA.size.count.txt"))
  m.c = read.delim(paste0("../../../output_union/tables/", s, ".matmiRNA.size.count.txt"))
  tot = ifelse(grepl("Cumulus", s), 67905011, 72511711)
  sb = sum(m.b$Freq); sc = sum(m.c$Freq)
  mb = as.data.table(m.b)[, .(Var1, Freq)]
  mc = as.data.table(m.c)[, .(Var1, Freq)]
  added = merge(mc, mb, by = "Var1", all = TRUE); added[is.na(added)] = 0
  setnames(added, c("Var1", "Freq.y", "Freq.x"))
  added = added[, added := Freq.y - Freq.x][added > 0]
  setorder(added, -added)
  cat(sprintf("%s: matmiRNA within=%.0f (%.1f%%) | matmiRNA union=%.0f (%.1f%%) | added=%.0f (%.1f%% of within)",
              s, sb, 100*sb/tot, sc, 100*sc/tot, sc-sb, 100*(sc-sb)/max(sb,1)))
  cat("  added-read sizes (top): ")
  cat(paste(sprintf("%snt=%.0f", added$Var1[1:6], added$added[1:6]), collapse = ", "), "\n")
  cat.tab = rbind(cat.tab, data.table(sample = s, within_reads = sb, union_reads = sc,
                                      added_reads = sc - sb))
}
write.csv(cat.tab, "../../../output_union/tables/Table3b_union_added_matmiRNA_reads.csv", row.names = FALSE)

## ---- 4) top mature miRNAs under B vs C ---------------------------------------------
load("../../../DB/rdata/mm39.matmiRNA.gr.RData")
mm39.matmiRNA.gr = mm39.matmiRNA.gr[!duplicated(paste(seqnames(mm39.matmiRNA.gr),
                                                      start(mm39.matmiRNA.gr),
                                                      end(mm39.matmiRNA.gr),
                                                      strand(mm39.matmiRNA.gr)))]
expr.tab = data.table()
for(strategy in c("B_matmiRNA_within", "C_matmiRNA_union")){
  dir = ifelse(strategy == "B_matmiRNA_within", "../../../output_matmiRNA/rdata", "../../../output_union/rdata")
  files = list.files(dir, pattern = ".bam.annotated.gr.RData$")
  for(f in files){
    load(paste0(dir, "/", f))
    s = gsub("\\.bam\\.annotated\\.gr\\.RData", "", f)
    g = reads.bam.annotated.gr[reads.bam.annotated.gr$type == "matmiRNA"]
    if(length(g) == 0) next
    ol = findOverlaps(g, mm39.matmiRNA.gr, type = "any")
    hit = data.table(sample = s, strategy = strategy,
                     name = mm39.matmiRNA.gr$Name[subjectHits(ol)],
                     rid = queryHits(ol),
                     count = as.numeric(g$count[queryHits(ol)]))
    hit = hit[, .(n_unique = length(unique(rid)), n_reads = base::sum(count)), by = .(sample, strategy, name)]
    expr.tab = rbind(expr.tab, hit)
  }
}
setorder(expr.tab, sample, strategy, -n_reads)
write.csv(expr.tab, "../../../output_union/tables/Table3c_mature_miRNA_expression_B_vs_C.csv", row.names = FALSE)
for(s in unique(expr.tab$sample)){
  for(st in c("B_matmiRNA_within", "C_matmiRNA_union")){
    top = expr.tab[sample == s & strategy == st][1:10]
    out("-- %s | %s top mature miRNAs (reads) --", s, st)
    for(k in 1:nrow(top)) out("   %-22s %10.0f", top$name[k], top$n_reads[k])
  }
}

## ---- 5) strand specificity for strategy C -------------------------------------------
for(s in unique(t2b.list$matmiRNA$sample)){
  ann = read.delim(paste0("../../../output_union/tables/", s, ".read.annotation.count.txt"))
  dt = as.data.table(ann)
  dt[, type2 := sub("^AS\\.", "", Var1)]
  dt[, str2 := ifelse(grepl("^AS\\.", Var1), "AS", "S")]
  w = dcast(dt, type2 ~ str2, value.var = "Freq", fun.aggregate = base::sum, fill = 0)
  w[, ratio := round(S/pmax(AS, 1), 1)]
  setorder(w, -S)
  out("-- %s strand specificity (strategy C, top 6) --", s)
  for(k in 1:min(6, nrow(w))) out("   %-18s sense=%.0f antisense=%.0f (ratio %.0f)", w$type2[k], w$S[k], w$AS[k], w$ratio[k])
}

out("Three-way comparison written to output_union/figures/ and output_union/tables/")
