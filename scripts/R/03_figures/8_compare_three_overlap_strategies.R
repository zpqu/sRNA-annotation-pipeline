##Date: 2026
##Author: Zhipeng
## Step 8 (optional): three-way comparison of SENSE OVERLAP RULES on the SAME
## mature-miRNA (matmiRNA) feature set. The three strategies differ ONLY in the
## sense overlap rule:
##   A: matmiRNA + any overlap    (read overlaps feature by >= 1 bp)  -> output_any/
##   B: matmiRNA + within         (read fully contained in feature)   -> output_matmiRNA/
##   C: matmiRNA + union          (read within feature OR feature within read) -> output_union/
## Contents (covers ALL genome-feature categories, not only miRNA):
##   1) full read.annotation composition (unique + all reads) across A/B/C
##   2) per-category read-size distributions (matmiRNA, snoRNA, piRNA, tRNA)
##   3) per-category read totals and how reads move between categories under A/B/C
##   4) mature-miRNA-focused: top identities, 5p/3p, sizes, union/any-added reads
##   5) strand specificity for strategy A
## Outputs are written to output_any/tables and output_any/figures.

library(ggplot2)
library(scales)
library(data.table)
library(GenomicRanges)

out = function(...) cat(sprintf(...), "\n")

strategy.list = c("A_matmiRNA_any" = "output_any",
                  "B_matmiRNA_within" = "output_matmiRNA",
                  "C_matmiRNA_union" = "output_union")
cols = c(A_matmiRNA_any = "#e31a1c", B_matmiRNA_within = "#33a02c", C_matmiRNA_union = "#1f78b4")

## ---- 1) full read.annotation composition ------------------------------------------
d.all = data.table()
for(st in names(strategy.list)){
  dir = strategy.list[[st]]
  t2a = read.csv(paste0("../../../", dir, "/tables/Table2a_annotation_count_unique_reads.csv"))
  t2b = read.csv(paste0("../../../", dir, "/tables/Table2b_annotation_count_all_reads.csv"))
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
write.csv(d.all[order(strategy, -pct_reads)], "../../../output_any/tables/Table4a_overlap_rule_composition.csv", row.names = FALSE)

p.comp = ggplot(d.all, aes(x = item, y = pct_reads, fill = strategy)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  facet_grid(sample ~ ., scales = "free_y") +
  labs(title = "read.annotation composition - all reads (overlap-rule strategies)",
       y = "% of reads", x = NULL) +
  scale_fill_manual(values = cols) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(p.comp, file = "../../../output_any/figures/Comparison1_overlap_rules_composition.pdf",
       width = 13, height = 6)

## ---- 2) per-category read-size distributions ----------------------------------------
size.all = data.table()
for(st in names(strategy.list)){
  dir = strategy.list[[st]]
  for(cat in c("matmiRNA", "snoRNA", "piRNA", "tRNA")){
    f = list.files(paste0("../../../", dir, "/tables"), pattern = paste0("cells.", cat, ".size.count.txt$"))
    for(ff in f){
      x = read.delim(paste0("../../../", dir, "/tables/", ff))
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
ggsave(p.size, file = "../../../output_any/figures/Comparison2_overlap_rules_category_size.pdf",
       width = 12, height = 10)

## ---- 3) per-category totals + read-movement between strategies ----------------------
cat.tot = dcast(d.all, sample + item ~ strategy, value.var = "pct_reads", fill = 0)
cat.tot[, dAB := A_matmiRNA_any - B_matmiRNA_within]
cat.tot[, dCB := C_matmiRNA_union - B_matmiRNA_within]
write.csv(cat.tot[order(sample, -B_matmiRNA_within)], "../../../output_any/tables/Table4b_overlap_rule_category_totals.csv", row.names = FALSE)

## read-movement cross-tabulation (per read position: type in A vs type in B)
mov.tab = data.table()
for(s in c("Cumulus-cells", "Granulosa-cells")){
  gr.list = list()
  for(st in names(strategy.list)){
    load(paste0("../../../", strategy.list[[st]], "/rdata/", s, ".bam.annotated.gr.RData"))
    gr.list[[st]] = data.table(poskey = paste(seqnames(reads.bam.annotated.gr), start(reads.bam.annotated.gr),
                                              end(reads.bam.annotated.gr), strand(reads.bam.annotated.gr)),
                               type = reads.bam.annotated.gr$type,
                               count = as.numeric(reads.bam.annotated.gr$count))
  }
  m = merge(gr.list$A_matmiRNA_any, gr.list$B_matmiRNA_within[, .(poskey, typeB = type)],
            by = "poskey", all = TRUE)
  m = merge(m, gr.list$C_matmiRNA_union[, .(poskey, typeC = type)], by = "poskey", all = TRUE)
  m[is.na(typeB), typeB := "other"]
  m[is.na(typeC), typeC := "other"]
  setnames(m, "type", "typeA")
  mov.tab = rbind(mov.tab, data.table(sample = s, m[, .(typeA, typeB, typeC, count)]))
}
## aggregate read movement (all reads) between B and A
mov.BA = mov.tab[, .(n_reads = base::sum(count)), by = .(typeB, typeA)]
setorder(mov.BA, -n_reads)
write.csv(mov.BA, "../../../output_any/tables/Table4c_read_movement_B_vs_A.csv", row.names = FALSE)
## between B and C
mov.BC = mov.tab[, .(n_reads = base::sum(count)), by = .(typeB, typeC)]
setorder(mov.BC, -n_reads)
write.csv(mov.BC, "../../../output_any/tables/Table4d_read_movement_B_vs_C.csv", row.names = FALSE)

out("== read movement B(matmiRNA-within) -> A(matmiRNA-any) [top transitions] ==")
for(k in 1:min(10, nrow(mov.BA)))
  out("  %-14s -> %-14s %12.0f", mov.BA$typeB[k], mov.BA$typeA[k], mov.BA$n_reads[k])
out("== read movement B(matmiRNA-within) -> C(matmiRNA-union) [top transitions] ==")
for(k in 1:min(10, nrow(mov.BC)))
  out("  %-14s -> %-14s %12.0f", mov.BC$typeB[k], mov.BC$typeC[k], mov.BC$n_reads[k])

## ---- 4) mature-miRNA-focused -------------------------------------------------------
load("../../../DB/rdata/mm39.matmiRNA.gr.RData")
mm39.matmiRNA.gr = mm39.matmiRNA.gr[!duplicated(paste(seqnames(mm39.matmiRNA.gr),
                                                      start(mm39.matmiRNA.gr),
                                                      end(mm39.matmiRNA.gr),
                                                      strand(mm39.matmiRNA.gr)))]
expr.tab = data.table()
for(st in names(strategy.list)){
  dir = strategy.list[[st]]
  files = list.files(paste0("../../../", dir, "/rdata"), pattern = ".bam.annotated.gr.RData$")
  for(f in files){
    load(paste0("../../../", dir, "/rdata/", f))
    s = gsub("\\.bam\\.annotated\\.gr\\.RData", "", f)
    g = reads.bam.annotated.gr[reads.bam.annotated.gr$type == "matmiRNA"]
    if(length(g) == 0) next
    ol = findOverlaps(g, mm39.matmiRNA.gr, type = "any")
    hit = data.table(sample = s, strategy = st,
                     name = mm39.matmiRNA.gr$Name[subjectHits(ol)],
                     rid = queryHits(ol),
                     count = as.numeric(g$count[queryHits(ol)]))
    hit = hit[, .(n_unique = length(unique(rid)), n_reads = base::sum(count)),
              by = .(sample, strategy, name)]
    expr.tab = rbind(expr.tab, hit)
  }
}
setorder(expr.tab, sample, strategy, -n_reads)
write.csv(expr.tab, "../../../output_any/tables/Table4e_mature_miRNA_expression_A_B_C.csv", row.names = FALSE)
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

## ---- 5) strand specificity for strategy A -------------------------------------------
for(s in c("Cumulus-cells", "Granulosa-cells")){
  ann = read.delim(paste0("../../../output_any/tables/", s, ".read.annotation.count.txt"))
  dt = as.data.table(ann)
  dt[, type2 := sub("^AS\\.", "", Var1)]
  dt[, str2 := ifelse(grepl("^AS\\.", Var1), "AS", "S")]
  w = dcast(dt, type2 ~ str2, value.var = "Freq", fun.aggregate = base::sum, fill = 0)
  w[, ratio := round(S/pmax(AS, 1), 1)]
  setorder(w, -S)
  out("-- %s strand specificity (strategy A, top 6) --", s)
  for(k in 1:min(6, nrow(w))) out("   %-18s sense=%.0f antisense=%.0f (ratio %.0f)", w$type2[k], w$S[k], w$AS[k], w$ratio[k])
}

out("Overlap-rule comparison written to output_any/tables/ and output_any/figures/")
