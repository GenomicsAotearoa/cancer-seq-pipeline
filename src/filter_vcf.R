#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(VariantAnnotation)
library(rtracklayer)


if (length(args)==0) {
  stop("VCF file need to be supplied\n", call.=FALSE)
} else if (length(args)==1) {
  
  file.gz <- args[1]
  stopifnot(file.exists(file.gz))
  file.gz.tbi <- paste(file.gz, ".tbi", sep="")
  if(!(file.exists(file.gz.tbi)))
    indexTabix(file.gz, format="vcf")
  
  ## Start define filters
  
  allelicDepth <- function(x) {
    ad <- geno(x)$AD
    tDepth <- geno(x)$AD[, 2] + geno(x)$RD[, 2]
    nDepth <- geno(x)$AD[, 1] + geno(x)$RD[, 1]
    
    tVarDepth <- geno(x)$AD[, 2]
    tVarFreq <- geno(x)$AD[, 2]/(geno(x)$AD[, 2]+geno(x)$RD[, 2])
    
    test <- (tDepth >= 50) & (nDepth >= 50) & (tVarDepth >= 10) & (tVarFreq >= 0.25)
    as.vector(!is.na(test) & test)
  }
  
  NotInComplexRegion <- function(x){
    grepl("wgEncodeDukeMapabilityRegionsExcludable=.", x, fixed=TRUE) &
      grepl("wgEncodeDacMapabilityConsensusExcludable=.", x, fixed=TRUE) &
      !grepl("SUPPLEMENTARY", x, fixed=TRUE) & grepl("SS=2", x, fixed = TRUE)
  }
  prefilters <- FilterRules(list(region=NotInComplexRegion))
  filters <- FilterRules(list(AD=allelicDepth))
  destination.file <- tempfile()
  
  tabix.file <- TabixFile(file.gz, yieldSize=10000)
  f1 <- filterVcf(tabix.file, "hg19", destination.file, prefilters = prefilters,filters = filters, verbose=TRUE)
  
  vcf.f1 <- readVcf(f1, "hg19")
  exome.bed <- import("S07604514_Regions.bed", format="bed")
  snp.in.exome.index <-  queryHits(findOverlaps(vcf.f1, exome.bed))
  snp.in.exome <- vcf.f1[snp.in.exome.index]
  
  chrom <- seqnames(rowRanges(snp.in.exome))
  SNP_loc <- start(ranges(rowRanges(snp.in.exome)))
  control_BAF <- as.numeric(gsub("%","",geno(snp.in.exome)$FREQ[,1]))/100
  tumor_BAF <- as.numeric(gsub("%","",geno(snp.in.exome)$FREQ[,2]))/100
  control_doc <- geno(snp.in.exome)$DP[,1]
  tumor_doc <- geno(snp.in.exome)$DP[,2]
  
  baf<-data.frame(list(chrom=chrom,SNP_loc=SNP_loc,control_BAF=control_BAF,tumor_BAF=tumor_BAF,control_doc=control_doc,tumor_doc=tumor_doc))
  write.table(file=gsub(".vcf.gz",".baf",file.gz),baf, row.names = F,sep="\t", quote=FALSE)
  
  }