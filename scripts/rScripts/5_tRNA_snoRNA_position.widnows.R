##Date: 25/07/2015
##Author: Zhipeng
## This script is used to make plot for small RNA position distribution with respect to tRNA

library(GenomicFeatures)
library(GenomicAlignments)
library(ggplot2)
library(scales)
rm(list = ls())

##format tRNAs
load("../../DB/rdata/mm10.tRNA.gr.RData")
mm10.tRNA.20bp.gr = mm10.tRNA.gr[-c(1:length(mm10.tRNA.gr))]
for(i in 1:length(mm10.tRNA.gr)){
      test.gr = mm10.tRNA.gr[i]
      if(as.character(strand(test.gr)) == "+"){
      test.old.start = seq(start(test.gr), end(test.gr), by = 1)
      test.new.end = seq(start(test.gr)+20, end(test.gr), by = 1)
      test.new.start = test.old.start[1:length(test.new.end)]
      }else{
      test.old.end = seq(end(test.gr), start(test.gr), by = -1)
      test.new.start = seq(end(test.gr)-20, start(test.gr), by = -1)
      test.new.end = test.old.end[1:length(test.new.start)]
      }
      test.new.chr = rep(seqnames(test.gr), length(test.new.start))
      test.new.strand = rep(strand(test.gr), length(test.new.start))
      test.new.name = rep(values(test.gr)$gene_id, length(test.new.start))
      test.new.pos = c(1:length(test.new.start))
      test.new.gr = GRanges(seqnames = test.new.chr, 
      		  ranges = IRanges(start = test.new.start, end = test.new.end),
		  strand = test.new.strand, position = test.new.pos, 
		  gene_id = test.new.name)
      mm10.tRNA.20bp.gr = c(mm10.tRNA.20bp.gr, test.new.gr)
}
save(mm10.tRNA.20bp.gr, file = "../../rdata/mm10.tRNA.20bp.gr.RData")

##format snoRNAs
load("../../DB/rdata/mm10.snoRNA.gr.RData")
mm10.snoRNA.gr=mm10.snoRNA.gr[width(mm10.snoRNA.gr) > 20]
mm10.snoRNA.20bp.gr = mm10.snoRNA.gr[-c(1:length(mm10.snoRNA.gr))]
for(i in 1:length(mm10.snoRNA.gr)){
      test.gr = mm10.snoRNA.gr[i]
      if(as.character(strand(test.gr)) == "+"){
      test.old.start = seq(start(test.gr), end(test.gr), by = 1)
      test.new.end = seq(start(test.gr)+20, end(test.gr), by = 1)
      test.new.start = test.old.start[1:length(test.new.end)]
      }else{
      test.old.end = seq(end(test.gr), start(test.gr), by = -1)
      test.new.start = seq(end(test.gr)-20, start(test.gr), by = -1)
      test.new.end = test.old.end[1:length(test.new.start)]
      }
      test.new.chr = rep(seqnames(test.gr), length(test.new.start))
      test.new.strand = rep(strand(test.gr), length(test.new.start))
      test.new.name = rep(values(test.gr)$gene_id, length(test.new.start))
      test.new.pos = c(1:length(test.new.start))
      test.new.gr = GRanges(seqnames = test.new.chr, 
      		  ranges = IRanges(start = test.new.start, end = test.new.end),
		  strand = test.new.strand, position = test.new.pos, 
		  gene_id = test.new.name)
      mm10.snoRNA.20bp.gr = c(mm10.snoRNA.20bp.gr, test.new.gr)
}
save(mm10.snoRNA.20bp.gr, file = "../../rdata/mm10.snoRNA.20bp.gr.RData")

###
tRNA.dis.all.df = NULL
snoRNA.dis.all.df = NULL
files = list.files(path = "../../rdata/", pattern = ".bam.annotated.gr.RData$")
for(i in seq(along = files)){
      file.name = paste("../../rdata/", files[i], sep = "")
      sample.name = files[i]
      sample.name = gsub("\\.bam\\.annotated\\.gr\\.RData", "", sample.name)
      print(paste("Now is processing ...", sample.name, Sys.time()))
      load(file.name)
      
      ##get tRNA coordinates
      mm10.tRNA.20bp.single.df = NULL
      values(mm10.tRNA.20bp.gr)$sense = NULL
      values(mm10.tRNA.20bp.gr)$antisense = NULL 
      values(mm10.tRNA.20bp.gr)$sample = NULL
      
      test.bam.tRNA.gr = reads.bam.annotated.gr[reads.bam.annotated.gr$type == "tRNA"]
      test.bam.tRNA_AS.gr = reads.bam.annotated.gr[reads.bam.annotated.gr$type == "AS.tRNA"]
      
      values(mm10.tRNA.20bp.gr)$sample = rep(sample.name, length(mm10.tRNA.20bp.gr))
      mm10.tRNA.20bp.bam.ol = countOverlaps(mm10.tRNA.20bp.gr, test.bam.tRNA.gr)
      values(mm10.tRNA.20bp.gr)$sense = mm10.tRNA.20bp.bam.ol
      mm10.tRNA.20bp.bam_AS.ol = countOverlaps(mm10.tRNA.20bp.gr, test.bam.tRNA_AS.gr, ignore.strand = T)
      values(mm10.tRNA.20bp.gr)$antisense = mm10.tRNA.20bp.bam_AS.ol
      mm10.tRNA.20bp.single.df = as.data.frame(mm10.tRNA.20bp.gr)
      tRNA.dis.all.df = rbind(tRNA.dis.all.df, mm10.tRNA.20bp.single.df)

      ##get snoRNA coordinates
      mm10.snoRNA.20bp.single.df = NULL
      values(mm10.snoRNA.20bp.gr)$sense = NULL
      values(mm10.snoRNA.20bp.gr)$antisense = NULL 
      values(mm10.snoRNA.20bp.gr)$sample = NULL
      
      test.bam.snoRNA.gr = reads.bam.annotated.gr[reads.bam.annotated.gr$type == "snoRNA"]
      test.bam.snoRNA_AS.gr = reads.bam.annotated.gr[reads.bam.annotated.gr$type == "AS.snoRNA"]
      
      values(mm10.snoRNA.20bp.gr)$sample = rep(sample.name, length(mm10.snoRNA.20bp.gr))
      mm10.snoRNA.20bp.bam.ol = countOverlaps(mm10.snoRNA.20bp.gr, test.bam.snoRNA.gr)
      values(mm10.snoRNA.20bp.gr)$sense = mm10.snoRNA.20bp.bam.ol
      mm10.snoRNA.20bp.bam_AS.ol = countOverlaps(mm10.snoRNA.20bp.gr, test.bam.snoRNA_AS.gr, ignore.strand = T)
      values(mm10.snoRNA.20bp.gr)$antisense = mm10.snoRNA.20bp.bam_AS.ol
      mm10.snoRNA.20bp.single.df = as.data.frame(mm10.snoRNA.20bp.gr)
      snoRNA.dis.all.df = rbind(snoRNA.dis.all.df, mm10.snoRNA.20bp.single.df)
}

save(tRNA.dis.all.df, file = "../../rdata/tRNA.20bp.dis.all.df.RData")
save(snoRNA.dis.all.df, file = "../../rdata/snoRNA.20bp.dis.all.df.RData")

#####plot tRNA
sample.num = length(unique(tRNA.dis.all.df$sample))
if(sample.num < 4){
              fig.width = sample.num*4
              fig.height = 4
}else{
        fig.height = (as.integer((sample.num-1)/4)+1)*4
        fig.width = 16
}
tRNA.dis.all.df$class = paste(tRNA.dis.all.df$position, tRNA.dis.all.df$sample, sep = "-")
test.antisense.df=test.sense.df=NULL
test.sense.df = as.data.frame(tapply(tRNA.dis.all.df$sense, as.factor(tRNA.dis.all.df$class), mean))
names(test.sense.df) = "freq"
test.sense.df$position = as.numeric(gsub("\\-.+", "", rownames(test.sense.df)))
test.sense.df$sample = gsub(".+\\-", "", rownames(test.sense.df))
test.antisense.df = as.data.frame(tapply(tRNA.dis.all.df$antisense, as.factor(tRNA.dis.all.df$class), mean))
names(test.antisense.df) = "freq"
test.antisense.df$position = as.numeric(gsub("\\-.+", "", rownames(test.antisense.df)))
test.antisense.df$sample = gsub(".+\\-", "", rownames(test.antisense.df))
tRNA.pos.freq.df = rbind(test.sense.df, test.antisense.df)
#tRNA.pos.freq.df$position = as.numeric(tRNA.pos.freq.df$position)

p.tRNA.pos.barplot = ggplot(data = test.sense.df, aes(x = position, y = freq)) +
                   geom_bar(stat = "identity") +
                   xlab("") + ylab("Count") +
                   scale_y_continuous(labels = comma) +
                   theme_bw() +
                   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                   facet_wrap(~sample, ncol = 4)
ggsave(p.tRNA.pos.barplot, file = "../../plots/Fig2c.tRNA_pos_barplot.pdf", width = fig.width, height = fig.height)

#p.tRNA.pos.percentage.plot = ggplot(data = test.sense.df, aes(x = position, y = freq, group = sample)) +
#                   geom_histogram(aes(y = ..percentage..), stat = "identity") +
#                   xlab("") + ylab("Percentage") +
#                   scale_y_continuous(labels = comma) +
#                   theme_bw() +
#                   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#                   facet_wrap(~sample, ncol = 4)
#ggsave(p.tRNA.pos.percentage.plot, file = "../../plots/Fig2b.tRNA_pos_percentage_plot.pdf", width = fig.width, height = fig.height)

##antisense
p.tRNA.pos.AS.barplot = ggplot(data = test.antisense.df, aes(x = position, y = freq)) +
                   geom_bar(stat = "identity") +
                   xlab("") + ylab("Count") +
                   scale_y_continuous(labels = comma) +
                   theme_bw() +
                   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                   facet_wrap(~sample, ncol = 4)
ggsave(p.tRNA.pos.AS.barplot, file = "../../plots/Fig2c.tRNA_pos_AS_barplot.pdf", width = fig.width, height = fig.height)

#p.tRNA.pos.AS.percentage.plot = ggplot(data = test.antisense.df, aes(x = position, y = freq, group = sample)) +
#                   geom_histogram(aes(y = ..percentage..), stat = "identity") +
#                   xlab("") + ylab("Percentage") +
#                   scale_y_continuous(labels = comma) +
#                   theme_bw() +
#                   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#                   facet_wrap(~sample, ncol = 4)
#ggsave(p.tRNA.pos.AS.percentage.plot, file = "../../plots/Fig2d.tRNA_pos_AS_percentage_plot.pdf", width = fig.width, height = fig.height)


#####plot snoRNA
sample.num = length(unique(snoRNA.dis.all.df$sample))
if(sample.num < 4){
              fig.width = sample.num*4
              fig.height = 4
}else{
        fig.height = (as.integer((sample.num-1)/4)+1)*4
        fig.width = 16
}
snoRNA.dis.all.df$class = paste(snoRNA.dis.all.df$position, snoRNA.dis.all.df$sample, sep = "-")
test.antisense.df=test.sense.df=NULL
test.sense.df = as.data.frame(tapply(snoRNA.dis.all.df$sense, as.factor(snoRNA.dis.all.df$class), mean))
names(test.sense.df) = "freq"
test.sense.df$position = as.numeric(gsub("\\-.+", "", rownames(test.sense.df)))
test.sense.df$sample = gsub(".+\\-", "", rownames(test.sense.df))
test.antisense.df = as.data.frame(tapply(snoRNA.dis.all.df$antisense, as.factor(snoRNA.dis.all.df$class), mean))
names(test.antisense.df) = "freq"
test.antisense.df$position = as.numeric(gsub("\\-.+", "", rownames(test.antisense.df)))
test.antisense.df$sample = gsub(".+\\-", "", rownames(test.antisense.df))
snoRNA.pos.freq.df = rbind(test.sense.df, test.antisense.df)
#snoRNA.pos.freq.df$position = as.numeric(snoRNA.pos.freq.df$position)

p.snoRNA.pos.barplot = ggplot(data = test.sense.df, aes(x = position, y = freq)) +
                   geom_bar(stat = "identity") +
                   xlab("") + ylab("Count") +
                   scale_y_continuous(labels = comma) +
                   theme_bw() +
                   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                   facet_wrap(~sample, ncol = 4)
ggsave(p.snoRNA.pos.barplot, file = "../../plots/Fig2c.snoRNA_pos_barplot.pdf", width = fig.width, height = fig.height)

#p.snoRNA.pos.percentage.plot = ggplot(data = test.sense.df, aes(x = position, y = freq*100)) +
#                   geom_bar(stat = "identity") +
#                   xlab("") + ylab() +
#                   scale_y_continuous(labels = percent) +
#                   theme_bw() +
#                   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#                   facet_wrap(~sample, ncol = 4)
#ggsave(p.snoRNA.pos.percentage.plot, file = "../../plots/Fig2b.snoRNA_pos_percentage_plot.pdf", width = fig.width, height = fig.height)

##antisense
p.snoRNA.pos.AS.barplot = ggplot(data = test.antisense.df, aes(x = position, y = freq)) +
                   geom_bar(stat = "identity") +
                   xlab("") + ylab("Count") +
                   scale_y_continuous(labels = comma) +
                   theme_bw() +
                   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                   facet_wrap(~sample, ncol = 4)
ggsave(p.snoRNA.pos.AS.barplot, file = "../../plots/Fig2c.snoRNA_pos_AS_barplot.pdf", width = fig.width, height = fig.height)

#p.snoRNA.pos.AS.percentage.plot = ggplot(data = test.antisense.df, aes(x = position, y = freq*100)) +
#                   geom_bar(stat = "identity") +
#                   xlab("") + ylab() +
#                   scale_y_continuous(labels = percent) +
#                   theme_bw() +
#                   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#                   facet_wrap(~sample, ncol = 4)
#ggsave(p.snoRNA.pos.AS.percentage.plot, file = "../../plots/Fig2d.snoRNA_pos_AS_percentage_plot.pdf", width = fig.width, height = fig.height)



