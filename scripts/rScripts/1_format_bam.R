##Date: 16/07/2015
##Author: Zhipeng
## This script is used to format small RNA data as rdata

library(rtracklayer)
library(GenomicFeatures)
library(Rsamtools)
library(GenomicAlignments)
rm(list = ls())

files = list.files(path = "../../bams/", pattern = ".bam$")
for(i in seq(along = files)){
  file.name = paste("../../bams/", files[i], sep = "");
  sam = files[i];
  print(paste("Now is processing ...", sam, "at", Sys.time()))
  sam = gsub("\\.bam", "", sam)
  sam = gsub("\\_R1\\.trimmed\\_NEB\\.bwa\\.hg19", "", sam)
  reads.bam = readGAlignments(file.name, format = "BAM", use.names = T)
  reads.strand = gsub("chr", "", seqnames(reads.bam))
  reads.strand = gsub("Chr", "", reads.strand)
  reads.strand = paste("chr", reads.strand, sep = "")
  reads.bam.gr = GRanges(seqnames = reads.strand,
			    ranges = ranges(reads.bam), 
			    strand = strand(reads.bam))
  gr.name = paste("../../rdata/", sam, ".bam.gr.RData", sep = "");
  save(reads.bam.gr, file = gr.name)
  gc()
}