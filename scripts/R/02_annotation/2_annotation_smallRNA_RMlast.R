##Date: 21/07/2015
##Author: Zhipeng
## This script is used to make annotations for small RNA data (reference genome
## set in config/genome.R)
## for all small RNAs, priority: miRNA>snoRNA>piRNA>tRNA>RM>refGene_NM_exon>refGene_NM_intron>lincRNA
## >antisense_tRNA>antisense_RM>antisense_NM_exon>antisense_NM_intron>anitsense_lincRNA
## for specific smallRNAs, priority: CDS>5UTR>3UTR>intron>up1k>down1k>RM
## >antisense_CDS>antisense_5UTR>antisense_3UTR>antisense_intron>antisense_up1k>antisense_down1k>antisense_RM>intergenic
## Overlap strategy: sense annotation uses fully-contained overlap (type="within",
## reads must be entirely contained in the feature); antisense annotation uses the
## default overlap strategy (type="any", minimum overlap of 1 bp).
## Reads are annotated as non-redundant unique reads (chr/start/end/strand); their
## frequencies (the 'count' metadata column from step 1) are summed when counting.
## Two consolidated count tables are saved for step 2:
##   output/Table2a_annotation_count_unique_reads.csv  (each unique read once)
##   output/Table2b_annotation_count_all_reads.csv     (weighted by read frequency)
## plus the per-sample .txt count files used by the figure scripts.

library(GenomicFeatures)
library(data.table)
rm(list = ls())

source("../../../config/genome.R")

###load genomic features
load(file.path(db.dir, "primiRNA.gr.RData"))
load(file.path(db.dir, "snoRNA.gr.RData"))
load(file.path(db.dir, "piRNA.gr.RData"))
load(file.path(db.dir, "tRNA.gr.RData"))
load(file.path(db.dir, "RM.gr.RData"))
load(file.path(db.dir, "refGene.NM.exon.grl.RData"))
load(file.path(db.dir, "refGene.NM.intron.grl.RData"))
load(file.path(db.dir, "lincRNA.exon.grl.RData"))

#load additional gene features
load(file.path(db.dir, "refGene.NM.CDS.grl.RData"))
load(file.path(db.dir, "refGene.NM.5UTR.grl.RData"))
load(file.path(db.dir, "refGene.NM.3UTR.grl.RData"))
load(file.path(db.dir, "refGene.NM.up1k.grl.RData"))
load(file.path(db.dir, "refGene.NM.down1k.grl.RData"))

####make genomic feature list or gene feature list
genomicFeature.id = c("primiRNA.gr", "snoRNA.gr", "piRNA.gr",
		  "tRNA.gr", "RM.gr", "refGene.NM.exon.grl",
		  "refGene.NM.intron.grl",
		  "lincRNA.exon.grl")
geneFeature.list = c("refGene.NM.CDS.grl", "refGene.NM.5UTR.grl",
	       "refGene.NM.3UTR.grl", "refGene.NM.intron.grl",
	       "refGene.NM.up1k.grl", "refGene.NM.down1k.grl", "RM.gr")

###
source("mygeneFeature.R")

### count unique reads weighted by their frequency (sum of the 'count' column)
count.by = function(group, counts){
  dt = data.table(Var1 = as.character(group), Freq = as.numeric(counts))
  dt = dt[, .(Freq = sum(Freq)), by = Var1]
  return(as.data.frame(dt))
}

### count per group: n_unique (each unique read once) and n_reads (weighted by count)
count.by.2 = function(group, counts){
  dt = data.table(Var1 = as.character(group), count = as.numeric(counts))
  dt = dt[, .(n_unique = .N, n_reads = sum(count)), by = Var1]
  return(dt)
}

### append group counts to the two consolidated tables (unique-based and all-reads-based)
unique.all.tab = data.table()
reads.all.tab = data.table()
add.counts = function(category, group, counts){
  dt = count.by.2(group, counts)
  unique.all.tab <<- rbind(unique.all.tab,
                           data.table(sample = sample, category = category,
                                      item = dt$Var1, Freq = dt$n_unique))
  reads.all.tab <<- rbind(reads.all.tab,
                          data.table(sample = sample, category = category,
                                     item = dt$Var1, Freq = dt$n_reads))
}

files = list.files(path = "../../../output/rdata", pattern = ".bam.unique.gr.RData$")
for(i in seq(along = files)){
      file.name = paste("../../../output/rdata/", files[i], sep = "")
      sample = files[i]
      sample = gsub("\\.bam\\.unique\\.gr\\.RData", "", sample)
      sample = gsub("\\.(bwa|bowtie2)$", "", sample)

      print(paste("Now is processing ...", sample, Sys.time()))
      load(file.name)
      test.all.id = c(1:length(reads.bam.unique.gr))
      
      test.bam.gr = reads.bam.unique.gr
      test.bam.gr$type = NULL
      test.bam.gr$region = NULL
      ## position key for removing already-annotated reads. Names are NOT unique:
      ## multi-mapping reads share one name across loci (step 1 collapses by
      ## chr/start/end/strand), so removal must use this key, not read names.
      test.bam.gr$rm.key = paste(seqnames(test.bam.gr), start(test.bam.gr),
                                 end(test.bam.gr), strand(test.bam.gr), sep = ":")
      test.bam.new.gr = test.bam.gr
      test.bam.new.gr = test.bam.new.gr[-(1:length(test.bam.new.gr))]
      test.other.gr = NULL
      ###miRNA[0]>snoRNA[1]>piRNA[2]>tRNA[3]>RM_refGeneNM[4]>RM[5]>refGene_NM_exon[6]>refGene_NM_intron[7]>lincRNA[8]
      ### >antisense_tRNA>RM_refGeneNM>antisense_RM>antisense_NM_exon>antisense_NM_intron>anitsense_lincRNA

      ##sense
      for(j in 1:length(genomicFeature.id)){
      	    feature.name = genomicFeature.id[j]
	    feature.id = gsub("\\.grl", "", feature.name)
	    feature.id = gsub("\\.gr", "", feature.id)
	    
	    test.bam.ol = NULL
	    test.bam.hit.gr = NULL
	    test.bam.ol = findOverlaps(test.bam.gr, get(feature.name), type = "within")
	    test.bam.hit.gr = test.bam.gr[unique(queryHits(test.bam.ol))]
	    test.bam.hit.gr$region = rep("NA", length(test.bam.hit.gr))
	    if(feature.id == "primiRNA"){
	    		  test.bam.hit.gr = mygeneFeature(bam = "test.bam.hit.gr")
			  primiRNA.geneFeature.count.df = count.by(test.bam.hit.gr$region, test.bam.hit.gr$count)
			  write.table(primiRNA.geneFeature.count.df, file = paste("../../../output/tables/", sample, ".primiRNA.annotation.count.txt", sep = ""),
			   sep = "\t", quote = F, col.names = T, row.names = F)
			  primiRNA.size.count.df = count.by(width(test.bam.hit.gr), test.bam.hit.gr$count)
			  write.table(primiRNA.size.count.df, file = paste("../../../output/tables/", sample, ".primiRNA.size.count.txt", sep = ""),
			   sep = "\t", quote = F, col.names = T, row.names = F)
			  add.counts("primiRNA.annotation", test.bam.hit.gr$region, test.bam.hit.gr$count)
			  add.counts("primiRNA.size", width(test.bam.hit.gr), test.bam.hit.gr$count)
	    }
	    if(feature.id == "snoRNA"){
	    		  test.bam.hit.gr = mygeneFeature("test.bam.hit.gr")
			  snoRNA.geneFeature.count.df = count.by(test.bam.hit.gr$region, test.bam.hit.gr$count)
			  write.table(snoRNA.geneFeature.count.df, file = paste("../../../output/tables/", sample, ".snoRNA.annotation.count.txt", sep = ""),
			   sep = "\t", quote = F, col.names = T, row.names = F)
			  snoRNA.size.count.df = count.by(width(test.bam.hit.gr), test.bam.hit.gr$count)
			  write.table(snoRNA.size.count.df, file = paste("../../../output/tables/", sample, ".snoRNA.size.count.txt", sep = ""),
			   sep = "\t", quote = F, col.names = T, row.names = F)
			  add.counts("snoRNA.annotation", test.bam.hit.gr$region, test.bam.hit.gr$count)
			  add.counts("snoRNA.size", width(test.bam.hit.gr), test.bam.hit.gr$count)
	    }
	    if(feature.id == "piRNA"){
	    		  test.bam.hit.gr = mygeneFeature("test.bam.hit.gr")
			  piRNA.geneFeature.count.df = count.by(test.bam.hit.gr$region, test.bam.hit.gr$count)
			  write.table(piRNA.geneFeature.count.df, file = paste("../../../output/tables/", sample, ".piRNA.annotation.count.txt", sep = ""),
			   sep = "\t", quote = F, col.names = T, row.names = F)
			  piRNA.size.count.df = count.by(width(test.bam.hit.gr), test.bam.hit.gr$count)
			  write.table(piRNA.size.count.df, file = paste("../../../output/tables/", sample, ".piRNA.size.count.txt", sep = ""),
			   sep = "\t", quote = F, col.names = T, row.names = F)
			  add.counts("piRNA.annotation", test.bam.hit.gr$region, test.bam.hit.gr$count)
			  add.counts("piRNA.size", width(test.bam.hit.gr), test.bam.hit.gr$count)
	    }
	    if(feature.id == "tRNA"){
	    		  test.bam.hit.gr = mygeneFeature("test.bam.hit.gr")
			  tRNA.geneFeature.count.df = count.by(test.bam.hit.gr$region, test.bam.hit.gr$count)
			  write.table(tRNA.geneFeature.count.df, file = paste("../../../output/tables/", sample, ".tRNA.annotation.count.txt", sep = ""),
			   sep = "\t", quote = F, col.names = T, row.names = F)
			  tRNA.size.count.df = count.by(width(test.bam.hit.gr), test.bam.hit.gr$count)
			  write.table(tRNA.size.count.df, file = paste("../../../output/tables/", sample, ".tRNA.size.count.txt", sep = ""),
			   sep = "\t", quote = F, col.names = T, row.names = F)
			  add.counts("tRNA.annotation", test.bam.hit.gr$region, test.bam.hit.gr$count)
			  add.counts("tRNA.size", width(test.bam.hit.gr), test.bam.hit.gr$count)
	    }
	    test.bam.hit.gr$type = rep(feature.id, length(test.bam.hit.gr))
	    test.bam.new.gr = c(test.bam.new.gr, test.bam.hit.gr)
	    if(length(unique(queryHits(test.bam.ol))) == 0){
	     test.bam.gr = test.bam.gr
	    }else{
	     test.bam.gr = test.bam.gr[!test.bam.gr$rm.key %in% test.bam.gr$rm.key[unique(queryHits(test.bam.ol))]]
	    }
      }

      ##antisense
      for(k in 1:length(genomicFeature.id)){
            feature.name = genomicFeature.id[k]
            feature.id = gsub("\\.grl", "", feature.name)
            feature.id = gsub("\\.gr", "", feature.id)

            test.bam.ol = NULL
            test.bam.hit.gr = NULL
            test.bam.ol = findOverlaps(test.bam.gr, get(feature.name), type = "any", ignore.strand = T)
            test.bam.hit.gr = test.bam.gr[unique(queryHits(test.bam.ol))]
	    test.bam.hit.gr$region = rep("NA", length(test.bam.hit.gr))
            test.bam.hit.gr$type = rep(paste("AS.", feature.id, sep = ""), length(test.bam.hit.gr))
            test.bam.new.gr = c(test.bam.new.gr, test.bam.hit.gr)
            if(length(unique(queryHits(test.bam.ol))) == 0){
	        test.bam.gr = test.bam.gr  
            }else{
	    	test.bam.gr = test.bam.gr[!test.bam.gr$rm.key %in% test.bam.gr$rm.key[unique(queryHits(test.bam.ol))]]
            }
  
      }

      ##other      
      test.bam.other.gr = test.bam.gr
      test.bam.other.gr$region = rep("NA", length(test.bam.other.gr))
      test.bam.other.gr$type = rep("other", length(test.bam.other.gr))
      test.bam.new.gr = c(test.bam.new.gr, test.bam.other.gr)
      read.annotation.count.df = count.by(test.bam.new.gr$type, test.bam.new.gr$count)
      read.size.count.df = count.by(width(test.bam.new.gr), test.bam.new.gr$count)
      add.counts("read.annotation", test.bam.new.gr$type, test.bam.new.gr$count)
      add.counts("read.size", width(test.bam.new.gr), test.bam.new.gr$count)
      reads.bam.annotated.gr = test.bam.new.gr
      reads.bam.annotated.gr$rm.key = NULL
      save(reads.bam.annotated.gr, file = paste("../../../output/rdata/", sample, ".bam.annotated.gr.RData", sep = ""))
      write.table(read.annotation.count.df, file = paste("../../../output/tables/", sample, ".read.annotation.count.txt", sep = ""), quote = F, sep = "\t", col.names = T, row.names = F)
      write.table(read.size.count.df, file = paste("../../../output/tables/", sample, ".read.size.count.txt", sep = ""), quote = F, sep = "\t", col.names = T, row.names = F)
}

## ----- consolidated count tables (unique reads vs all reads) -------------------
write.csv(unique.all.tab, "../../../output/tables/Table2a_annotation_count_unique_reads.csv", row.names = FALSE)
write.csv(reads.all.tab, "../../../output/tables/Table2b_annotation_count_all_reads.csv", row.names = FALSE)
print("Tables saved to ../../../output/tables/: Table2a_annotation_count_unique_reads.csv, Table2b_annotation_count_all_reads.csv")