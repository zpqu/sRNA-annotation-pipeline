##Date: 21/07/2015
##Author: Zhipeng
## This script is used to make annotations for small RNA data (human)
## for all small RNAs, priority: miRNA>snoRNA>piRNA>tRNA>RM>refGene_NM_exon>refGene_NM_intron>lincRNA
## >antisense_tRNA>antisense_RM>antisense_NM_exon>antisense_NM_intron>anitsense_lincRNA
## for specific smallRNAs, priority: CDS>5UTR>3UTR>intron>up1k>down1k>RM
## >antisense_CDS>antisense_5UTR>antisense_3UTR>antisense_intron>antisense_up1k>antisense_down1k>antisense_RM>intergenic

library(GenomicFeatures)
rm(list = ls())

###load genomic features
load("../../DB/rdata/mm10.primiRNA.gr.RData")
load("../../DB/rdata/mm10.snoRNA.gr.RData")
load("../../DB/rdata/mm10.piRNA.gr.RData")
load("../../DB/rdata/mm10.tRNA.gr.RData")
load("../../DB/rdata/mm10.RM.gr.RData")
load("../../DB/rdata/mm10.refGene.NM.exon.grl.RData")
load("../../DB/rdata/mm10.refGene.NM.intron.grl.RData")
load("../../DB/rdata/mm10.lincRNA.exon.grl.RData")

#load additional gene features
load("../../DB/rdata/mm10.refGene.NM.CDS.grl.RData")
load("../../DB/rdata/mm10.refGene.NM.5UTR.grl.RData")
load("../../DB/rdata/mm10.refGene.NM.3UTR.grl.RData")
load("../../DB/rdata/mm10.refGene.NM.up1k.grl.RData")
load("../../DB/rdata/mm10.refGene.NM.down1k.grl.RData")

####make genomic feature list or gene feature list
genomicFeature.id = c("mm10.primiRNA.gr", "mm10.snoRNA.gr", "mm10.piRNA.gr",
		  "mm10.tRNA.gr", "mm10.RM.gr", "mm10.refGene.NM.exon.grl",
		  "mm10.refGene.NM.intron.grl",
		  "mm10.lincRNA.exon.grl")
geneFeature.list = c("mm10.refGene.NM.CDS.grl", "mm10.refGene.NM.5UTR.grl",
	       "mm10.refGene.NM.3UTR.grl", "mm10.refGene.NM.intron.grl",
	       "mm10.refGene.NM.up1k.grl", "mm10.refGene.NM.down1k.grl", "mm10.RM.gr")

###
source("mygeneFeature.R")

files = list.files(path = "../../rdata", pattern = ".bam.gr.RData$")
for(i in seq(along = files)){
      file.name = paste("../../rdata/", files[i], sep = "")
      sample = files[i]
      sample = gsub("\\.bam\\.gr\\.RData", "", sample)

      print(paste("Now is processing ...", sample, Sys.time()))
      load(file.name)
      test.all.id = c(1:length(reads.bam.gr))
      
      test.bam.gr = reads.bam.gr
      test.bam.gr$type = NULL
      test.bam.gr$region = NULL
      test.bam.new.gr = test.bam.gr
      test.bam.new.gr = test.bam.new.gr[-(1:length(test.bam.new.gr))]
      test.other.gr = NULL
      ###miRNA[0]>snoRNA[1]>piRNA[2]>tRNA[3]>RM_refGeneNM[4]>RM[5]>refGene_NM_exon[6]>refGene_NM_intron[7]>lincRNA[8]
      ### >antisense_tRNA>RM_refGeneNM>antisense_RM>antisense_NM_exon>antisense_NM_intron>anitsense_lincRNA

      ##sense
      for(j in 1:length(genomicFeature.id)){
      	    feature.name = genomicFeature.id[j]
	    feature.id = gsub("mm10\\.", "", feature.name)
	    feature.id = gsub("\\.grl", "", feature.id)
	    feature.id = gsub("\\.gr", "", feature.id)
	    
	    test.bam.ol = NULL
	    test.bam.hit.gr = NULL
	    test.bam.ol = findOverlaps(test.bam.gr, get(feature.name))
	    test.bam.hit.gr = test.bam.gr[unique(queryHits(test.bam.ol))]
	    test.bam.hit.gr$region = rep("NA", length(test.bam.hit.gr))
	    if(feature.id == "primiRNA"){
	    		  test.bam.hit.gr = mygeneFeature(bam = "test.bam.hit.gr")
			  primiRNA.geneFeature.count.df = as.data.frame(table(test.bam.hit.gr$region))
			  write.table(primiRNA.geneFeature.count.df, file = paste("../../output/", sample, ".primiRNA.annotation.count.txt", sep = ""),
			   sep = "\t", quote = F, col.names = T, row.names = F)
			  primiRNA.size.count.df = as.data.frame(table(width(test.bam.hit.gr)))
			  write.table(primiRNA.size.count.df, file = paste("../../output/", sample, ".primiRNA.size.count.txt", sep = ""),
			   sep = "\t", quote = F, col.names = T, row.names = F)
	    }
	    if(feature.id == "snoRNA"){
	    		  test.bam.hit.gr = mygeneFeature("test.bam.hit.gr")
			  snoRNA.geneFeature.count.df = as.data.frame(table(test.bam.hit.gr$region))
			  write.table(snoRNA.geneFeature.count.df, file = paste("../../output/", sample, ".snoRNA.annotation.count.txt", sep = ""),
			   sep = "\t", quote = F, col.names = T, row.names = F)
			  snoRNA.size.count.df = as.data.frame(table(width(test.bam.hit.gr)))
			  write.table(snoRNA.size.count.df, file = paste("../../output/", sample, ".snoRNA.size.count.txt", sep = ""),
			   sep = "\t", quote = F, col.names = T, row.names = F)
	    }
	    if(feature.id == "piRNA"){
	    		  test.bam.hit.gr = mygeneFeature("test.bam.hit.gr")
			  piRNA.geneFeature.count.df = as.data.frame(table(test.bam.hit.gr$region))
			  write.table(piRNA.geneFeature.count.df, file = paste("../../output/", sample, ".piRNA.annotation.count.txt", sep = ""),
			   sep = "\t", quote = F, col.names = T, row.names = F)
			  piRNA.size.count.df = as.data.frame(table(width(test.bam.hit.gr)))
			  write.table(piRNA.size.count.df, file = paste("../../output/", sample, ".piRNA.size.count.txt", sep = ""),
			   sep = "\t", quote = F, col.names = T, row.names = F)
	    }
	    if(feature.id == "tRNA"){
	    		  test.bam.hit.gr = mygeneFeature("test.bam.hit.gr")
			  tRNA.geneFeature.count.df = as.data.frame(table(test.bam.hit.gr$region))
			  write.table(tRNA.geneFeature.count.df, file = paste("../../output/", sample, ".tRNA.annotation.count.txt", sep = ""),
			   sep = "\t", quote = F, col.names = T, row.names = F)
			  tRNA.size.count.df = as.data.frame(table(width(test.bam.hit.gr)))
			  write.table(tRNA.size.count.df, file = paste("../../output/", sample, ".tRNA.size.count.txt", sep = ""),
			   sep = "\t", quote = F, col.names = T, row.names = F)
	    }
	    test.bam.hit.gr$type = rep(feature.id, length(test.bam.hit.gr))
	    test.bam.new.gr = c(test.bam.new.gr, test.bam.hit.gr)
	    if(length(unique(queryHits(test.bam.ol))) == 0){
	     test.bam.gr = test.bam.gr
	    }else{
	     test.bam.gr = test.bam.gr[!(names(test.bam.gr) %in% names(test.bam.gr[unique(queryHits(test.bam.ol))]))]
	    }
      }

      ##antisense
      for(k in 1:length(genomicFeature.id)){
            feature.name = genomicFeature.id[k]
            feature.id = gsub("mm10\\.", "", feature.name)
            feature.id = gsub("\\.grl", "", feature.id)
            feature.id = gsub("\\.gr", "", feature.id)

            test.bam.ol = NULL
            test.bam.hit.gr = NULL
            test.bam.ol = findOverlaps(test.bam.gr, get(feature.name), ignore.strand = T)
            test.bam.hit.gr = test.bam.gr[unique(queryHits(test.bam.ol))]
	    test.bam.hit.gr$region = rep("NA", length(test.bam.hit.gr))
            test.bam.hit.gr$type = rep(paste("AS.", feature.id, sep = ""), length(test.bam.hit.gr))
            test.bam.new.gr = c(test.bam.new.gr, test.bam.hit.gr)
            if(length(unique(queryHits(test.bam.ol))) == 0){
	        test.bam.gr = test.bam.gr  
            }else{
	    	test.bam.gr = test.bam.gr[!(names(test.bam.gr) %in% names(test.bam.gr[unique(queryHits(test.bam.ol))]))]
            }
  
      }

      ##other      
      test.bam.other.gr = test.bam.gr
      test.bam.other.gr$region = rep("NA", length(test.bam.other.gr))
      test.bam.other.gr$type = rep("other", length(test.bam.other.gr))
      test.bam.new.gr = c(test.bam.new.gr, test.bam.other.gr)
      read.annotation.count.df = as.data.frame(table(test.bam.new.gr$type))
      read.size.count.df = as.data.frame(table(width(test.bam.new.gr)))
      reads.bam.annotated.gr = test.bam.new.gr
      save(reads.bam.annotated.gr, file = paste("../../rdata/", sample, ".bam.annotated.gr.RData", sep = ""))
      write.table(read.annotation.count.df, file = paste("../../output/", sample, ".read.annotation.count.txt", sep = ""), quote = F, sep = "\t", col.names = T, row.names = F)
      write.table(read.size.count.df, file = paste("../../output/", sample, ".read.size.count.txt", sep = ""), quote = F, sep = "\t", col.names = T, row.names = F)
}