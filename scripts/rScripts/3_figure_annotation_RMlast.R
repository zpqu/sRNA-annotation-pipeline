##Date: 21/07/2015
##Author: Zhipeng
## This script is used to make plots for small RNA annotation results (step 3).
## It reads the two consolidated count tables produced by step 2:
##   output/Table2a_annotation_count_unique_reads.csv   (unique reads)
##   output/Table2b_annotation_count_all_reads.csv      (all reads)
## and draws annotation barplots for the "read" and the primiRNA/snoRNA/piRNA/
## tRNA gene-feature classes, showing the unique-reads and all-reads flavors
## side by side (dodged bars). Counts use a log10 y-axis because the annotation
## categories span orders of magnitude; a percentage panel is added beneath.
## Output figures follow the step naming rule: Figure3a, Figure3b, ...

library(ggplot2)
library(scales)
library(data.table)
library(patchwork)

###read consolidated count tables from step 2
unique.tab = read.csv("../../output/Table2a_annotation_count_unique_reads.csv")
all.tab = read.csv("../../output/Table2b_annotation_count_all_reads.csv")
unique.tab$flavor = "unique reads"
all.tab$flavor = "all reads"
all.df = rbind(unique.tab, all.tab)

class.map = c(read = "read.annotation",
              primiRNA = "primiRNA.annotation",
              snoRNA = "snoRNA.annotation",
              piRNA = "piRNA.annotation",
              tRNA = "tRNA.annotation")
class.list = names(class.map)

read.levels = c("primiRNA", "snoRNA", "piRNA", "tRNA", "RM", "refGene.NM.exon",
                "refGene.NM.intron", "lincRNA.exon", "AS.primiRNA", "AS.snoRNA",
                "AS.piRNA", "AS.tRNA", "AS.RM", "AS.refGene.NM.exon",
                "AS.refGene.NM.intron", "AS.lincRNA.exon", "other")
region.levels = c("CDS", "5UTR", "3UTR", "intron", "up1k", "down1k", "RM",
                  "AS.CDS", "AS.5UTR", "AS.3UTR", "AS.intron", "AS.up1k",
                  "AS.down1k", "AS.RM", "intergenic")

rot.theme = theme(axis.text.x = element_text(angle = 45, hjust = 1))
flavor.cols = c("unique reads" = "#1f78b4", "all reads" = "#e31a1c")

for(j in 1:length(class.list)){
      class = class.list[j]
      category = class.map[[class]]
      sub.dt = as.data.table(all.df[all.df$category == category, ])
      if(nrow(sub.dt) == 0) next
      sub.dt$item = factor(sub.dt$item, levels = if(class == "read") read.levels else region.levels)
      sample.num = length(unique(sub.dt$sample))
      fig.width = max(8, sample.num * 4)
      fig.height = 8.5

      ##percentage panel: complete all item x sample x flavor combos (zero-filled)
      all.combo = expand.grid(sample = unique(sub.dt$sample),
                              item = levels(sub.dt$item),
                              flavor = c("unique reads", "all reads"),
                              stringsAsFactors = FALSE)
      sub.full = merge(sub.dt[, .(sample, item, flavor, Freq)], all.combo,
                       by = c("sample", "item", "flavor"), all = TRUE)
      sub.full$Freq[is.na(sub.full$Freq)] = 0
      sub.full$item = factor(sub.full$item, levels = levels(sub.dt$item))
      sub.full = as.data.table(sub.full)
      sub.full[, per := Freq / sum(Freq), by = .(sample, flavor)]

      p.count = ggplot(sub.dt[Freq > 0], aes(x = item, y = Freq, fill = flavor)) +
            geom_bar(stat = "identity", position = position_dodge(preserve = "single"), width = 0.7) +
            scale_y_log10(labels = comma) +
            scale_fill_manual(values = flavor.cols, name = "Reads") +
            facet_wrap(~sample, ncol = 4) +
            labs(title = paste0(class, " annotation - count"),
                 y = "Count (log10 scale)", x = NULL) +
            theme_bw() + rot.theme

      p.per = ggplot(sub.full, aes(x = item, y = per, fill = flavor)) +
            geom_bar(stat = "identity", position = "dodge", width = 0.7) +
            scale_y_continuous(labels = percent) +
            scale_fill_manual(values = flavor.cols, name = "Reads") +
            facet_wrap(~sample, ncol = 4) +
            labs(title = paste0(class, " annotation - percentage"),
                 y = "Percentage", x = NULL) +
            theme_bw() + rot.theme

      fig.base = paste0("../../plots/Figure3", letters[j], ".", class, "_annotation_barplot")
      pdf(paste0(fig.base, ".pdf"), width = fig.width, height = fig.height)
      print(p.count / p.per)
      dev.off()
      png(paste0(fig.base, ".png"), width = fig.width * 300, height = fig.height * 300, res = 300)
      print(p.count / p.per)
      dev.off()
      print(paste("saved:", fig.base, Sys.time()))
}
