##Date: 21/07/2015
##Author: Zhipeng
## Modified: mature-miRNA strategy (2b)
## Same annotation logic as 2_annotation_smallRNA_RMlast.R, but the miRNA
## category is annotated with MATURE miRNA loci (matmiRNA.gr, miRBase
## v22, ~22 nt) instead of primary miRNA transcripts (primiRNA).
## Priority: matmiRNA>snoRNA>piRNA>tRNA>RM>refGene_NM_exon>refGene_NM_intron>lincRNA
## >antisense_tRNA>antisense_RM>antisense_NM_exon>antisense_NM_intron>antisense_lincRNA
## Reads are annotated as non-redundant unique reads (chr/start/end/strand).
## All outputs are written to output_matmiRNA/ so the primary-miRNA strategy
## results in output/ are kept intact for comparison.
## The reference genome (and the feature DB directory) is set in config/genome.R.
##
## Extra output: per-read annotation + abundance tables for unique reads
## (output_matmiRNA/tables/<sample>.unique_reads_annotation.csv and
## Table_unique_reads_annotation.csv), one row per unique read, sorted by
## expression (count) descending.

library(GenomicFeatures)
library(data.table)
rm(list = ls())

source("../../../config/genome.R")

###load genomic features (mature miRNA replaces primary miRNA)
load(file.path(db.dir, "matmiRNA.gr.RData"))
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
genomicFeature.id = c("matmiRNA.gr", "snoRNA.gr", "piRNA.gr",
		  "tRNA.gr", "RM.gr", "refGene.NM.exon.grl",
		  "refGene.NM.intron.grl",
		  "lincRNA.exon.grl")
geneFeature.list = c("refGene.NM.CDS.grl", "refGene.NM.5UTR.grl",
	       "refGene.NM.3UTR.grl", "refGene.NM.intron.grl",
	       "refGene.NM.up1k.grl", "refGene.NM.down1k.grl", "RM.gr")

### helper: id vector for a genomic feature, aligned to the subject index used
### by findOverlaps on that feature (element index for GRangesList features).
feature.id.vector = function(feature.name, tx2gn){
  f = get(feature.name)
  if(is(f, "GRangesList")){
    u = unlist(f)
    if(!is.null(u$transcript_id) && !is.null(u$gene_name)){
      gn = as.character(u$gene_name); gi = as.character(u$gene_id)
      gn[is.na(gn) | gn == ""] = gi[is.na(gn) | gi == ""]
      dt = unique(data.table(tx = as.character(u$transcript_id), id = gn))
      return(dt$id[match(names(f), dt$tx)])
    }
    return(tx2gn$gene_name[match(names(f), tx2gn$transcript_id)])
  }
  id.col = switch(feature.name,
    matmiRNA.gr = "Name", snoRNA.gr = "gene_id", piRNA.gr = "name",
    tRNA.gr = "gene_id", RM.gr = "gene_id")
  as.character(mcols(f)[[id.col]])
}
## transcript -> gene-name map (for intron grl, whose unlisted ranges carry no mcols)
tx2gn = unique(data.table(transcript_id = as.character(unlist(refGene.NM.exon.grl)$transcript_id),
                          gene_name = as.character(unlist(refGene.NM.exon.grl)$gene_name)))
fids.all = lapply(genomicFeature.id, feature.id.vector, tx2gn = tx2gn)
names(fids.all) = genomicFeature.id

### helper: for a findOverlaps result, return the specific feature id (and the
### number of overlapping features) for each hit read, aligned to
### unique(queryHits(ol)). piRNA uses the highest-scoring locus among the
### overlapping piRNAs; other features use an arbitrary (first) hit.
feature.hit.info = function(ol, feature.name, fids, type, ignore.strand = FALSE){
  qh = queryHits(ol); sh = subjectHits(ol)
  rows = unique(qh)
  if(length(rows) == 0) return(list(feature.id = character(0), n_features = integer(0)))
  nf = tabulate(qh, nbins = length(test.bam.gr))[rows]
  if(feature.name == "piRNA.gr"){
    o = order(qh, -as.numeric(mcols(piRNA.gr)$score[sh]))
    uq = qh[o]; us = sh[o]
    keep = !duplicated(uq)
    sub.map = us[keep]; names(sub.map) = uq[keep]
    fid = fids[sub.map[as.character(rows)]]
  }else{
    fid = fids[findOverlaps(test.bam.gr[rows], get(feature.name), type = type,
                            ignore.strand = ignore.strand, select = "arbitrary")]
  }
  list(feature.id = fid, n_features = nf)
}

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

sample.list = list()
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
	    hit.info = feature.hit.info(test.bam.ol, feature.name, fids.all[[feature.name]], "within")
	    test.bam.hit.gr$feature.id = hit.info$feature.id
	    test.bam.hit.gr$n_features = hit.info$n_features
	    if(feature.id == "matmiRNA"){
	    		  test.bam.hit.gr = mygeneFeature(bam = "test.bam.hit.gr")
			  matmiRNA.geneFeature.count.df = count.by(test.bam.hit.gr$region, test.bam.hit.gr$count)
			  write.table(matmiRNA.geneFeature.count.df, file = paste("../../../output_matmiRNA/tables/", sample, ".matmiRNA.annotation.count.txt", sep = ""),
			   sep = "\t", quote = F, col.names = T, row.names = F)
			  matmiRNA.size.count.df = count.by(width(test.bam.hit.gr), test.bam.hit.gr$count)
			  write.table(matmiRNA.size.count.df, file = paste("../../../output_matmiRNA/tables/", sample, ".matmiRNA.size.count.txt", sep = ""),
			   sep = "\t", quote = F, col.names = T, row.names = F)
			  add.counts("matmiRNA.annotation", test.bam.hit.gr$region, test.bam.hit.gr$count)
			  add.counts("matmiRNA.size", width(test.bam.hit.gr), test.bam.hit.gr$count)
	    }
	    if(feature.id == "snoRNA"){
	    		  test.bam.hit.gr = mygeneFeature("test.bam.hit.gr")
			  snoRNA.geneFeature.count.df = count.by(test.bam.hit.gr$region, test.bam.hit.gr$count)
			  write.table(snoRNA.geneFeature.count.df, file = paste("../../../output_matmiRNA/tables/", sample, ".snoRNA.annotation.count.txt", sep = ""),
			   sep = "\t", quote = F, col.names = T, row.names = F)
			  snoRNA.size.count.df = count.by(width(test.bam.hit.gr), test.bam.hit.gr$count)
			  write.table(snoRNA.size.count.df, file = paste("../../../output_matmiRNA/tables/", sample, ".snoRNA.size.count.txt", sep = ""),
			   sep = "\t", quote = F, col.names = T, row.names = F)
			  add.counts("snoRNA.annotation", test.bam.hit.gr$region, test.bam.hit.gr$count)
			  add.counts("snoRNA.size", width(test.bam.hit.gr), test.bam.hit.gr$count)
	    }
	    if(feature.id == "piRNA"){
	    		  test.bam.hit.gr = mygeneFeature("test.bam.hit.gr")
			  piRNA.geneFeature.count.df = count.by(test.bam.hit.gr$region, test.bam.hit.gr$count)
			  write.table(piRNA.geneFeature.count.df, file = paste("../../../output_matmiRNA/tables/", sample, ".piRNA.annotation.count.txt", sep = ""),
			   sep = "\t", quote = F, col.names = T, row.names = F)
			  piRNA.size.count.df = count.by(width(test.bam.hit.gr), test.bam.hit.gr$count)
			  write.table(piRNA.size.count.df, file = paste("../../../output_matmiRNA/tables/", sample, ".piRNA.size.count.txt", sep = ""),
			   sep = "\t", quote = F, col.names = T, row.names = F)
			  add.counts("piRNA.annotation", test.bam.hit.gr$region, test.bam.hit.gr$count)
			  add.counts("piRNA.size", width(test.bam.hit.gr), test.bam.hit.gr$count)
	    }
	    if(feature.id == "tRNA"){
	    		  test.bam.hit.gr = mygeneFeature("test.bam.hit.gr")
			  tRNA.geneFeature.count.df = count.by(test.bam.hit.gr$region, test.bam.hit.gr$count)
			  write.table(tRNA.geneFeature.count.df, file = paste("../../../output_matmiRNA/tables/", sample, ".tRNA.annotation.count.txt", sep = ""),
			   sep = "\t", quote = F, col.names = T, row.names = F)
			  tRNA.size.count.df = count.by(width(test.bam.hit.gr), test.bam.hit.gr$count)
			  write.table(tRNA.size.count.df, file = paste("../../../output_matmiRNA/tables/", sample, ".tRNA.size.count.txt", sep = ""),
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
	    hit.info = feature.hit.info(test.bam.ol, feature.name, fids.all[[feature.name]], "any", ignore.strand = T)
	    test.bam.hit.gr$feature.id = hit.info$feature.id
	    test.bam.hit.gr$n_features = hit.info$n_features
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
      test.bam.other.gr$feature.id = rep(NA_character_, length(test.bam.other.gr))
      test.bam.other.gr$n_features = rep(NA_integer_, length(test.bam.other.gr))
      test.bam.new.gr = c(test.bam.new.gr, test.bam.other.gr)
      read.annotation.count.df = count.by(test.bam.new.gr$type, test.bam.new.gr$count)
      read.size.count.df = count.by(width(test.bam.new.gr), test.bam.new.gr$count)
      add.counts("read.annotation", test.bam.new.gr$type, test.bam.new.gr$count)
      add.counts("read.size", width(test.bam.new.gr), test.bam.new.gr$count)
      reads.bam.annotated.gr = test.bam.new.gr
      reads.bam.annotated.gr$rm.key = NULL
      save(reads.bam.annotated.gr, file = paste("../../../output_matmiRNA/rdata/", sample, ".bam.annotated.gr.RData", sep = ""))
      write.table(read.annotation.count.df, file = paste("../../../output_matmiRNA/tables/", sample, ".read.annotation.count.txt", sep = ""), quote = F, sep = "\t", col.names = T, row.names = F)
      write.table(read.size.count.df, file = paste("../../../output_matmiRNA/tables/", sample, ".read.size.count.txt", sep = ""), quote = F, sep = "\t", col.names = T, row.names = F)

      ### per-read annotation + abundance table (unique reads, sorted by expression)
      sample.dt = data.table(
        read_id = seq_along(reads.bam.annotated.gr),
        chr = as.character(seqnames(reads.bam.annotated.gr)),
        start = start(reads.bam.annotated.gr),
        end = end(reads.bam.annotated.gr),
        strand = as.character(strand(reads.bam.annotated.gr)),
        size = width(reads.bam.annotated.gr),
        count = as.numeric(reads.bam.annotated.gr$count),
        cpm = round(as.numeric(reads.bam.annotated.gr$count) /
                      sum(as.numeric(reads.bam.annotated.gr$count)) * 1e6, 3),
        category = as.character(reads.bam.annotated.gr$type),
        feature_id = as.character(reads.bam.annotated.gr$feature.id),
        gene_context = as.character(reads.bam.annotated.gr$region),
        n_features = as.integer(reads.bam.annotated.gr$n_features))
      setorder(sample.dt, -count, chr, start, end)
      fwrite(sample.dt, file = paste("../../../output_matmiRNA/tables/", sample, ".unique_reads_annotation.csv", sep = ""))
      sample.list[[i]] = copy(sample.dt)[, sample := sample]
}

## ----- per-read annotation + abundance table, all samples combined ------------
all.samples.dt = rbindlist(sample.list)
fwrite(all.samples.dt, "../../../output_matmiRNA/tables/Table_unique_reads_annotation.csv")

## ----- consolidated count tables (unique reads vs all reads) -------------------
write.csv(unique.all.tab, "../../../output_matmiRNA/tables/Table2a_annotation_count_unique_reads.csv", row.names = FALSE)
write.csv(reads.all.tab, "../../../output_matmiRNA/tables/Table2b_annotation_count_all_reads.csv", row.names = FALSE)
print("Tables saved to ../../../output_matmiRNA/tables/: Table2a_annotation_count_unique_reads.csv, Table2b_annotation_count_all_reads.csv, Table_unique_reads_annotation.csv, <sample>.unique_reads_annotation.csv")
