##Date: 25/07/2015
##Author: Zhipeng
## This script makes plots for the position distribution of small RNA reads
## with respect to tRNA and snoRNA genes (step 5).
## Features for the reference genome set in config/genome.R are used. Reads are
## taken from the step-2 annotated objects (rdata/*.bam.annotated.gr.RData),
## "AS.tRNA", "snoRNA", "AS.snoRNA", ...) and the non-redundant read 'count'.
## Window counts are weighted by 'count'. For each gene a 21 bp window is slid
## in 1 bp steps along the gene body, and the mean weighted read count per
## window position is plotted per sample.

suppressMessages(library(GenomicRanges))
suppressMessages(library(GenomicAlignments))
suppressMessages(library(ggplot2))
suppressMessages(library(scales))
rm(list = ls())

source("../../../config/genome.R")

## ---- sliding 21 bp windows (1 bp step) along each gene body -----------------
sliding.windows = function(gr){
      win = gr[integer(0)]
      for(i in seq_along(gr)){
            g = gr[i]
            len = width(g)
            if(len <= 20) next
            if(as.character(strand(g)) == "+"){
                  starts = start(g):(end(g) - 20)
            }else{
                  starts = (end(g) - 20):start(g)
            }
            pos = seq_along(starts)
            win = c(win, GRanges(seqnames = seqnames(g),
                                 ranges = IRanges(start = starts, end = starts + 20),
                                 strand = strand(g),
                                 position = pos,
                                 gene_id = values(g)$gene_id))
      }
      win
}

## ---- sum the weighted count of reads overlapping each window -----------------
weighted.counts = function(wins, reads, cnt){
      out = integer(length(wins))
      if(length(reads) > 0){
            ol = findOverlaps(wins, reads, type = "any", ignore.strand = TRUE)
            cs = rowsum(cnt[subjectHits(ol)], queryHits(ol))
            out[as.integer(rownames(cs))] = cs[, 1]
      }
      out
}

## ---- format tRNAs -----------------------------------------------------------
load(file.path(db.dir, "tRNA.gr.RData"))
tRNA.20bp.gr = sliding.windows(tRNA.gr)
save(tRNA.20bp.gr, file = "../../../output/rdata/tRNA.20bp.gr.RData")
print(paste("tRNA windows:", length(tRNA.20bp.gr)))

## ---- format snoRNAs (genes longer than 20 bp) --------------------------------
load(file.path(db.dir, "snoRNA.gr.RData"))
snoRNA.gr = snoRNA.gr[width(snoRNA.gr) > 20]
snoRNA.20bp.gr = sliding.windows(snoRNA.gr)
save(snoRNA.20bp.gr, file = "../../../output/rdata/snoRNA.20bp.gr.RData")
print(paste("snoRNA windows:", length(snoRNA.20bp.gr)))

## ---- count reads per window, per sample --------------------------------------
tRNA.dis.all.df = NULL
snoRNA.dis.all.df = NULL
files = list.files(path = "../../../output/rdata/", pattern = ".bam.annotated.gr.RData$")
if(length(files) == 0) stop("no *.bam.annotated.gr.RData found in ../../../output/rdata/ (run step 2 first)")
for(i in seq(along = files)){
      file.name = paste("../../../output/rdata/", files[i], sep = "")
      sample.name = files[i]
      sample.name = gsub("\\.bam\\.annotated\\.gr\\.RData", "", sample.name)
      sample.name = gsub("\\.(bwa|bowtie2)$", "", sample.name)
      print(paste("Now is processing ...", sample.name, Sys.time()))
      load(file.name)
      cnt = mcols(reads.bam.annotated.gr)$count

      ## tRNA windows
      test.bam.tRNA.gr = reads.bam.annotated.gr[reads.bam.annotated.gr$type == "tRNA"]
      test.bam.tRNA_AS.gr = reads.bam.annotated.gr[reads.bam.annotated.gr$type == "AS.tRNA"]
      tRNA.20bp.gr$sense = weighted.counts(tRNA.20bp.gr, test.bam.tRNA.gr,
                                                cnt[reads.bam.annotated.gr$type == "tRNA"])
      tRNA.20bp.gr$antisense = weighted.counts(tRNA.20bp.gr, test.bam.tRNA_AS.gr,
                                                    cnt[reads.bam.annotated.gr$type == "AS.tRNA"])
      tRNA.20bp.gr$sample = rep(sample.name, length(tRNA.20bp.gr))
      tRNA.dis.all.df = rbind(tRNA.dis.all.df, as.data.frame(tRNA.20bp.gr))

      ## snoRNA windows
      test.bam.snoRNA.gr = reads.bam.annotated.gr[reads.bam.annotated.gr$type == "snoRNA"]
      test.bam.snoRNA_AS.gr = reads.bam.annotated.gr[reads.bam.annotated.gr$type == "AS.snoRNA"]
      snoRNA.20bp.gr$sense = weighted.counts(snoRNA.20bp.gr, test.bam.snoRNA.gr,
                                                  cnt[reads.bam.annotated.gr$type == "snoRNA"])
      snoRNA.20bp.gr$antisense = weighted.counts(snoRNA.20bp.gr, test.bam.snoRNA_AS.gr,
                                                      cnt[reads.bam.annotated.gr$type == "AS.snoRNA"])
      snoRNA.20bp.gr$sample = rep(sample.name, length(snoRNA.20bp.gr))
      snoRNA.dis.all.df = rbind(snoRNA.dis.all.df, as.data.frame(snoRNA.20bp.gr))
}

save(tRNA.dis.all.df, file = "../../../output/rdata/tRNA.20bp.dis.all.df.RData")
save(snoRNA.dis.all.df, file = "../../../output/rdata/snoRNA.20bp.dis.all.df.RData")

## ---- plot tRNA ---------------------------------------------------------------
position.plot = function(dis.all.df, fig.base){
      if(is.null(dis.all.df) || nrow(dis.all.df) == 0) return(NULL)
      sample.num = length(unique(dis.all.df$sample))
      if(sample.num < 4){
            fig.width = sample.num * 4
            fig.height = 4
      }else{
            fig.height = (as.integer((sample.num - 1) / 4) + 1) * 4
            fig.width = 16
      }
      dis.all.df$class = paste(dis.all.df$position, dis.all.df$sample, sep = "-")
      test.sense.df = as.data.frame(tapply(dis.all.df$sense,
                                           as.factor(dis.all.df$class), mean))
      names(test.sense.df) = "freq"
      test.sense.df$position = as.numeric(gsub("\\-.+", "", rownames(test.sense.df)))
      test.sense.df$sample = gsub(".+\\-", "", rownames(test.sense.df))
      test.antisense.df = as.data.frame(tapply(dis.all.df$antisense,
                                               as.factor(dis.all.df$class), mean))
      names(test.antisense.df) = "freq"
      test.antisense.df$position = as.numeric(gsub("\\-.+", "", rownames(test.antisense.df)))
      test.antisense.df$sample = gsub(".+\\-", "", rownames(test.antisense.df))

      p.sense = ggplot(data = test.sense.df, aes(x = position, y = freq)) +
            geom_bar(stat = "identity") +
            xlab("") + ylab("Mean count") +
            scale_y_continuous(labels = comma) +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            facet_wrap(~sample, ncol = 4)
      ggsave(p.sense, file = paste0(fig.base, ".pdf"), width = fig.width, height = fig.height)

      p.antisense = ggplot(data = test.antisense.df, aes(x = position, y = freq)) +
            geom_bar(stat = "identity") +
            xlab("") + ylab("Mean count") +
            scale_y_continuous(labels = comma) +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            facet_wrap(~sample, ncol = 4)
      ggsave(p.antisense, file = paste0(fig.base, "_AS.pdf"), width = fig.width, height = fig.height)
      print(paste("saved:", fig.base, Sys.time()))
}

position.plot(tRNA.dis.all.df, "../../../output/figures/Fig2c.tRNA_pos_barplot")
position.plot(snoRNA.dis.all.df, "../../../output/figures/Fig2c.snoRNA_pos_barplot")
