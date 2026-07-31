##Date: 2026
##Author: Zhipeng
## Step 6 (optional): compare the primary-miRNA (primiRNA) annotation strategy
## (output/) with the mature-miRNA (matmiRNA) strategy (output_matmiRNA/).
## 1) side-by-side read.annotation composition (unique vs all reads, both samples)
## 2) miRNA read-size distribution for the two strategies
## 3) top expressed mature miRNAs (identity resolved via mm39.matmiRNA.gr names)
## 4) fraction of primiRNA-annotated reads that are NOT at mature miRNA positions
## Outputs are written to output_matmiRNA/figures and output_matmiRNA/tables.

library(ggplot2)
library(scales)
library(data.table)
library(GenomicRanges)

out = function(...) cat(sprintf(...), "\n")

## ---- load both strategies -----------------------------------------------------
t2a.pri = read.csv("../../../output/tables/Table2a_annotation_count_unique_reads.csv")
t2b.pri = read.csv("../../../output/tables/Table2b_annotation_count_all_reads.csv")
t2a.mat = read.csv("../../../output_matmiRNA/tables/Table2a_annotation_count_unique_reads.csv")
t2b.mat = read.csv("../../../output_matmiRNA/tables/Table2b_annotation_count_all_reads.csv")

## item -> simpler label for the composition comparison
rename.item = function(x){
      x = gsub("^AS\\.", "AS_", x)
      x = gsub("refGene\\.NM\\.", "gene.", x)
      x = gsub("^other$", "other", x)
      x
}

mk.ra = function(t2a, t2b, strategy){
      ra.u = as.data.table(t2a)[category == "read.annotation"]
      ra.a = as.data.table(t2b)[category == "read.annotation"]
      d.u = ra.u[, .(sample, item, n_unique = Freq)][, pct_unique := 100*n_unique/sum(n_unique), by = sample]
      d.a = ra.a[, .(sample, item, n_reads = Freq)][, pct_reads := 100*n_reads/sum(n_reads), by = sample]
      d = merge(d.u, d.a, by = c("sample", "item"))
      d[, strategy := strategy]
      d
}
d.pri = mk.ra(t2a.pri, t2b.pri, "primiRNA")
d.mat = mk.ra(t2a.mat, t2b.mat, "matmiRNA")
d.all = rbind(d.pri, d.mat)
d.all[, item := rename.item(item)]

## ---- 1) composition comparison table + figure ---------------------------------
d.sort = d.all[order(-pct_reads)]
write.csv(d.sort, "../../../output_matmiRNA/tables/Table2n_strategy_comparison_read_annotation.csv", row.names = FALSE)

## pooled item levels (primiRNA-strategy items union matmiRNA-strategy items)
lv = c("matmiRNA", "primiRNA", "snoRNA", "piRNA", "tRNA", "RM", "gene.exon",
       "gene.intron", "lincRNA.exon", "AS_matmiRNA", "AS_primiRNA", "AS_snoRNA",
       "AS_piRNA", "AS_tRNA", "AS_RM", "AS_gene.exon", "AS_gene.intron",
       "AS_lincRNA.exon", "other")
d.all[, item := factor(item, levels = lv[!is.na(match(lv, unique(d.all$item)))])]

p.comp = ggplot(d.all, aes(x = item, y = pct_reads, fill = strategy)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
      facet_grid(sample ~ ., scales = "free_y") +
      labs(title = "read.annotation composition - all reads (primiRNA vs matmiRNA strategy)",
           y = "% of reads", x = NULL) +
      scale_fill_manual(values = c(primiRNA = "#1f78b4", matmiRNA = "#e31a1c")) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(p.comp, file = "../../../output_matmiRNA/figures/Comparison1_read_annotation_primiRNA_vs_matmiRNA.pdf",
       width = 12, height = 6)

## ---- 2) miRNA read-size distribution: primiRNA vs matmiRNA --------------------
size.pri = data.frame()
size.mat = data.frame()
for(s in unique(t2b.pri$sample)){
      f.pri = paste0("../../../output/tables/", s, ".primiRNA.size.count.txt")
      f.mat = paste0("../../../output_matmiRNA/tables/", s, ".matmiRNA.size.count.txt")
      a = read.delim(f.pri); a$sample = s; a$strategy = "primiRNA"
      b = read.delim(f.mat); b$sample = s; b$strategy = "matmiRNA"
      size.pri = rbind(size.pri, a)
      size.mat = rbind(size.mat, b)
}
size.all = rbind(size.pri, size.mat)
size.all = as.data.table(size.all)
size.all[, pct := as.numeric(100*Freq/sum(Freq)), by = .(sample, strategy)]
size.all[, strategy := factor(strategy, levels = c("primiRNA", "matmiRNA"))]

p.size = ggplot(size.all, aes(x = Var1, y = pct, fill = strategy)) +
      geom_bar(stat = "identity", position = "dodge", width = 0.8) +
      facet_wrap(~sample, ncol = 1) +
      labs(title = "miRNA-strategy read-size distribution (all reads)",
           x = "read length (nt)", y = "% of reads in the miRNA class") +
      scale_fill_manual(values = c(primiRNA = "#1f78b4", matmiRNA = "#e31a1c")) +
      theme_bw()
ggsave(p.size, file = "../../../output_matmiRNA/figures/Comparison2_miRNA_size_distribution.pdf",
       width = 8, height = 7)

## ---- 3) top expressed mature miRNAs --------------------------------------------
load("../../../DB/rdata/mm39.matmiRNA.gr.RData")
mm39.matmiRNA.gr = mm39.matmiRNA.gr[!duplicated(paste(seqnames(mm39.matmiRNA.gr),
                                                      start(mm39.matmiRNA.gr),
                                                      end(mm39.matmiRNA.gr),
                                                      strand(mm39.matmiRNA.gr)))]
mat.tab = data.table()
files = list.files("../../../output_matmiRNA/rdata", pattern = ".bam.annotated.gr.RData$")
for(f in files){
      load(paste0("../../../output_matmiRNA/rdata/", f))
      s = gsub("\\.bam\\.annotated\\.gr\\.RData", "", f)
      g = reads.bam.annotated.gr[reads.bam.annotated.gr$type == "matmiRNA"]
      if(length(g) == 0) next
      ol = findOverlaps(g, mm39.matmiRNA.gr, type = "any")
      hit = data.table(sample = s,
                       name = mm39.matmiRNA.gr$Name[subjectHits(ol)],
                       rid = queryHits(ol),
                       count = as.numeric(g$count[queryHits(ol)]))
      hit = hit[, .(n_unique = length(unique(rid)), n_reads = sum(count)), by = .(sample, name)]
      mat.tab = rbind(mat.tab, hit)
}
setorder(mat.tab, sample, -n_reads)
write.csv(mat.tab, "../../../output_matmiRNA/tables/Table2m_mature_miRNA_expression.csv", row.names = FALSE)
for(s in unique(mat.tab$sample)){
      top = mat.tab[sample == s][1:15]
      out("-- %s top mature miRNAs (reads) --", s)
      for(k in 1:nrow(top)) out("   %-22s %10.0f", top$name[k], top$n_reads[k])
}

## ---- 4) primiRNA reads NOT at mature positions ---------------------------------
## From the primiRNA strategy, how many reads fall inside pri-miRNA transcripts
## but outside all mature miRNA loci?
load("../../../DB/rdata/mm39.primiRNA.gr.RData")
load("../../../DB/rdata/mm39.matmiRNA.gr.RData")
for(f in files){
      load(paste0("../../../output_matmiRNA/rdata/", f))
      s = gsub("\\.bam\\.annotated\\.gr\\.RData", "", f)
      g = reads.bam.annotated.gr
      mat.reads = g[g$type == "matmiRNA"]
      ## total reads in primiRNA (strategy A) can be read from strategy-A tables
}
pri.a = as.data.table(t2b.pri)[category == "read.annotation"]
mat.a = as.data.table(t2b.mat)[category == "read.annotation"]
for(s in unique(pri.a$sample)){
      pri.mi = pri.a[sample == s & item == "primiRNA", sum(Freq)]
      mat.mi = mat.a[sample == s & item == "matmiRNA", sum(Freq)]
      out("-- %s: primiRNA-strategy miRNA reads = %.0f | matmiRNA-strategy miRNA reads = %.0f | non-mature reads within pri-miRNA = %.0f (%.1f%% of primiRNA category)",
          s, pri.mi, mat.mi, pri.mi - mat.mi, 100*(pri.mi - mat.mi)/pri.mi)
}

out("Comparison figures written to output_matmiRNA/figures/, tables to output_matmiRNA/tables/")
