##Date: 21/07/2015
##Author: Zhipeng
## Modified: mature-miRNA strategy with UNION containment sense annotation (2c)
## Same pipeline as 2b_annotation_smallRNA_matmiRNA_RMlast.R (miRNA category =
## MATURE miRNA loci, mm39.matmiRNA.gr, miRBase v22), but the SENSE annotation
## uses union containment: a read is assigned to a feature when EITHER the read
## is fully contained in the feature (read within feature) OR the feature is
## fully contained in the read (feature within read). This captures reads that
## span an entire small feature (e.g. a long read covering a 22 nt mature miRNA
## locus). Antisense annotation is unchanged (type="any", ignore.strand=T).
## Priority: matmiRNA>snoRNA>piRNA>tRNA>RM>refGene_NM_exon>refGene_NM_intron>lincRNA
## >antisense_tRNA>antisense_RM>antisense_NM_exon>antisense_NM_intron>antisense_lincRNA
## Reads are annotated as non-redundant unique reads (chr/start/end/strand).
## All outputs are written to output_union/ so that results of the primiRNA
## strategy (output/) and matmiRNA-within strategy (output_matmiRNA/) are kept
## intact for the three-way comparison.

library(GenomicFeatures)
library(data.table)
rm(list = ls())

###load genomic features (mature miRNA replaces primary miRNA)
load("../../../DB/rdata/mm39.matmiRNA.gr.RData")
load("../../../DB/rdata/mm39.snoRNA.gr.RData")
load("../../../DB/rdata/mm39.piRNA.gr.RData")
load("../../../DB/rdata/mm39.tRNA.gr.RData")
load("../../../DB/rdata/mm39.RM.gr.RData")
load("../../../DB/rdata/mm39.refGene.NM.exon.grl.RData")
load("../../../DB/rdata/mm39.refGene.NM.intron.grl.RData")
load("../../../DB/rdata/mm39.lincRNA.exon.grl.RData")

#load additional gene features
load("../../../DB/rdata/mm39.refGene.NM.CDS.grl.RData")
load("../../../DB/rdata/mm39.refGene.NM.5UTR.grl.RData")
load("../../../DB/rdata/mm39.refGene.NM.3UTR.grl.RData")
load("../../../DB/rdata/mm39.refGene.NM.up1k.grl.RData")
load("../../../DB/rdata/mm39.refGene.NM.down1k.grl.RData")

####make genomic feature list or gene feature list
genomicFeature.id = c("mm39.matmiRNA.gr", "mm39.snoRNA.gr", "mm39.piRNA.gr",
		  "mm39.tRNA.gr", "mm39.RM.gr", "mm39.refGene.NM.exon.grl",
		  "mm39.refGene.NM.intron.grl",
		  "mm39.lincRNA.exon.grl")
geneFeature.list = c("mm39.refGene.NM.CDS.grl", "mm39.refGene.NM.5UTR.grl",
	       "mm39.refGene.NM.3UTR.grl", "mm39.refGene.NM.intron.grl",
	       "mm39.refGene.NM.up1k.grl", "mm39.refGene.NM.down1k.grl", "mm39.RM.gr")

###
source("mygeneFeature.R")

### union containment: return indices of reads that are either fully contained in
### a feature (read within feature) or fully contain a feature (feature within read)
containment.union.ids = function(reads.gr, feature.gr){
  ol1 = findOverlaps(reads.gr, feature.gr, type = "within")
  ol2 = findOverlaps(feature.gr, reads.gr, type = "within")
  return(unique(c(queryHits(ol1), subjectHits(ol2))))
}

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

      ##sense (union containment: read within feature OR feature within read)
      for(j in 1:length(genomicFeature.id)){
      	    feature.name = genomicFeature.id[j]
	    feature.id = gsub("mm39\\.", "", feature.name)
	    feature.id = gsub("\\.grl", "", feature.id)
	    feature.id = gsub("\\.gr", "", feature.id)
	    
	    test.bam.ol = NULL
	    test.bam.hit.gr = NULL
	    hit.ids = containment.union.ids(test.bam.gr, get(feature.name))
	    test.bam.hit.gr = test.bam.gr[hit.ids]
	    test.bam.hit.gr$region = rep("NA", length(test.bam.hit.gr))
	    if(feature.id == "matmiRNA"){
	    		  test.bam.hit.gr = mygeneFeature(bam = "test.bam.hit.gr", type = "union")
			  matmiRNA.geneFeature.count.df = count.by(test.bam.hit.gr$region, test.bam.hit.gr$count)
			  write.table(matmiRNA.geneFeature.count.df, file = paste("../../../output_union/tables/", sample, ".matmiRNA.annotation.count.txt", sep = ""),
			   sep = "\t", quote = F, col.names = T, row.names = F)
			  matmiRNA.size.count.df = count.by(width(test.bam.hit.gr), test.bam.hit.gr$count)
			  write.table(matmiRNA.size.count.df, file = paste("../../../output_union/tables/", sample, ".matmiRNA.size.count.txt", sep = ""),
			   sep = "\t", quote = F, col.names = T, row.names = F)
			  add.counts("matmiRNA.annotation", test.bam.hit.gr$region, test.bam.hit.gr$count)
			  add.counts("matmiRNA.size", width(test.bam.hit.gr), test.bam.hit.gr$count)
	    }
	    if(feature.id == "snoRNA"){
	    		  test.bam.hit.gr = mygeneFeature("test.bam.hit.gr", type = "union")
			  snoRNA.geneFeature.count.df = count.by(test.bam.hit.gr$region, test.bam.hit.gr$count)
			  write.table(snoRNA.geneFeature.count.df, file = paste("../../../output_union/tables/", sample, ".snoRNA.annotation.count.txt", sep = ""),
			   sep = "\t", quote = F, col.names = T, row.names = F)
			  snoRNA.size.count.df = count.by(width(test.bam.hit.gr), test.bam.hit.gr$count)
			  write.table(snoRNA.size.count.df, file = paste("../../../output_union/tables/", sample, ".snoRNA.size.count.txt", sep = ""),
			   sep = "\t", quote = F, col.names = T, row.names = F)
			  add.counts("snoRNA.annotation", test.bam.hit.gr$region, test.bam.hit.gr$count)
			  add.counts("snoRNA.size", width(test.bam.hit.gr), test.bam.hit.gr$count)
	    }
	    if(feature.id == "piRNA"){
	    		  test.bam.hit.gr = mygeneFeature("test.bam.hit.gr", type = "union")
			  piRNA.geneFeature.count.df = count.by(test.bam.hit.gr$region, test.bam.hit.gr$count)
			  write.table(piRNA.geneFeature.count.df, file = paste("../../../output_union/tables/", sample, ".piRNA.annotation.count.txt", sep = ""),
			   sep = "\t", quote = F, col.names = T, row.names = F)
			  piRNA.size.count.df = count.by(width(test.bam.hit.gr), test.bam.hit.gr$count)
			  write.table(piRNA.size.count.df, file = paste("../../../output_union/tables/", sample, ".piRNA.size.count.txt", sep = ""),
			   sep = "\t", quote = F, col.names = T, row.names = F)
			  add.counts("piRNA.annotation", test.bam.hit.gr$region, test.bam.hit.gr$count)
			  add.counts("piRNA.size", width(test.bam.hit.gr), test.bam.hit.gr$count)
	    }
	    if(feature.id == "tRNA"){
	    		  test.bam.hit.gr = mygeneFeature("test.bam.hit.gr", type = "union")
			  tRNA.geneFeature.count.df = count.by(test.bam.hit.gr$region, test.bam.hit.gr$count)
			  write.table(tRNA.geneFeature.count.df, file = paste("../../../output_union/tables/", sample, ".tRNA.annotation.count.txt", sep = ""),
			   sep = "\t", quote = F, col.names = T, row.names = F)
			  tRNA.size.count.df = count.by(width(test.bam.hit.gr), test.bam.hit.gr$count)
			  write.table(tRNA.size.count.df, file = paste("../../../output_union/tables/", sample, ".tRNA.size.count.txt", sep = ""),
			   sep = "\t", quote = F, col.names = T, row.names = F)
			  add.counts("tRNA.annotation", test.bam.hit.gr$region, test.bam.hit.gr$count)
			  add.counts("tRNA.size", width(test.bam.hit.gr), test.bam.hit.gr$count)
	    }
	    test.bam.hit.gr$type = rep(feature.id, length(test.bam.hit.gr))
	    test.bam.new.gr = c(test.bam.new.gr, test.bam.hit.gr)
	    if(length(hit.ids) == 0){
	     test.bam.gr = test.bam.gr
	    }else{
	     test.bam.gr = test.bam.gr[!test.bam.gr$rm.key %in% test.bam.gr$rm.key[hit.ids]]
	    }
      }

      ##antisense
      for(k in 1:length(genomicFeature.id)){
            feature.name = genomicFeature.id[k]
            feature.id = gsub("mm39\\.", "", feature.name)
            feature.id = gsub("\\.grl", "", feature.id)
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
      save(reads.bam.annotated.gr, file = paste("../../../output_union/rdata/", sample, ".bam.annotated.gr.RData", sep = ""))
      write.table(read.annotation.count.df, file = paste("../../../output_union/tables/", sample, ".read.annotation.count.txt", sep = ""), quote = F, sep = "\t", col.names = T, row.names = F)
      write.table(read.size.count.df, file = paste("../../../output_union/tables/", sample, ".read.size.count.txt", sep = ""), quote = F, sep = "\t", col.names = T, row.names = F)
}

## ----- consolidated count tables (unique reads vs all reads) -------------------
write.csv(unique.all.tab, "../../../output_union/tables/Table2a_annotation_count_unique_reads.csv", row.names = FALSE)
write.csv(reads.all.tab, "../../../output_union/tables/Table2b_annotation_count_all_reads.csv", row.names = FALSE)
print("Tables saved to ../../../output_union/tables/: Table2a_annotation_count_unique_reads.csv, Table2b_annotation_count_all_reads.csv")
