## re-write of code written by Cris Print to identify mutations in a set of genes in a set of VCF files at low stringency

## Comments by Ben : ##


args <- commandArgs(trailingOnly = TRUE)
snpVCF=args[1]
indelVCF=args[2]
out1=args[3]
out2=args[4]


#source("https://bioconductor.org/biocLite.R")

suppressMessages(library("VariantAnnotation"))
suppressMessages(library("tidyverse"))
suppressMessages(library("here"))
suppressMessages(library("httr"))
suppressMessages(library("jsonlite"))
suppressMessages(library("reticulate"))



## Data location
## This section finds the files to be reported up, checks that there is a snp and an indel file for each sample
## and an associated tbi file in the same directory. If successful, it puts the full file paths
## into a dataframe with one row per samples.  

## data location. Files should be paired bzip and corresponding tbi file with
## sample name at the beginning of the file name.
## trailing / is required.
#snpVCF = "tmp/P1017A_1.fastq.gz.align.removeDuplicates.removeSuppplementary.bam.variantCalling.annotation.vcf.gz"
#indelVCF = "tmp/P1017A_1.fastq.gz.align.removeDuplicates.removeSuppplementary.bam.variantCalling.2.annotation.vcf.gz"


here()

metadata = read_tsv("../data/metadata.txt", col_types = cols())
mgl = metadata$Michelle.gene.list[!is.na(metadata$Michelle.gene.list)] # define Michelle's gene list
geneList = paste(mgl,collapse="|")

# prefilters, that are just greps on each line as they are read in.
regionFilter = function(x) grepl("exonic|splicing|UTR5|upstream", x, perl = TRUE)
geneFilter = function(x) grepl(paste0("\b",geneList,"\b"), x, perl = TRUE)

# filters on the info segments
geneFilter = function(x) {
  unlist(info(x)$Gene.refGene %in% mgl)
}

infoFilter = function(x) {
  blackList = info(x)$wgEncodeDacMapabilityConsensusExcludable == "."  # discard variants in Encode DAC black list regions
  germline = info(x)$SS == 1 # germline calls apparently.
  return(unlist(blackList & germline))
}

## filters on the geno region. 
genoFilter = function(x) {
  Ndepth = 10 # N overall depth 
  N.var.depth = 4 # N ALT supporting depth 
  nd = geno(x)$AD[,"NORMAL"] + geno(x)$RD[, "NORMAL"] >= Ndepth
  nvd = geno(x)$AD[,"NORMAL"] >= N.var.depth
  unlist(nd & nvd)
}

## contstruct the filter lists.
PF = FilterRules(list(regionFilter))
FF = FilterRules(list(geneFilter, infoFilter, genoFilter ))


## Process
## This section goes through the data frame holding the file locations, row by row, 
## filters, then extracts the data required for the report, returning a dataframe. 

report.df = read_csv("../data/reportColumns.csv")
infoVariables = report.df %>%
    filter(location == "info")



snpFilter = filterVcf(snpVCF , "hg19", tempfile(), filters = FF, prefilters = PF )
indelFilter = filterVcf(indelVCF , "hg19" , tempfile(), filters = FF, prefilters = PF )

snps = readVcf(snpFilter)
indels = readVcf(indelFilter) 


snp.info = as(info(snps),"DataFrame")
snp.info=as.tibble(snp.info)
snp.info = as.tibble(snp.info, rownames=NULL) %>%
    dplyr::select(one_of(infoVariables$variable))  %>%
    map(function(x) unlist(lapply(x, '[', 1))) %>%
    as.tibble() 
colnames(snp.info) =  infoVariables$name
snp.info$Mutation_type = "snp"
snp.ranges = rowRanges(snps)
snp.info$num_Germline_ALT_reads = as.character(unlist(snp.ranges$ALT))
snp.info$num_Germline_REF_reads = as.character(snp.ranges$REF)
snp.info$sample=paste0("sample","_",names(rowRanges(snps)))


indels.info = as(info(indels),"DataFrame")
indels.info=as.tibble(indels.info)
indels.info = as.tibble(indels.info, rownames=NULL) %>%
    dplyr::select(one_of(infoVariables$variable))  %>%
    map(function(x) unlist(lapply(x, '[', 1))) %>%
    as.tibble()
colnames(indels.info) =  infoVariables$name
indels.info$Mutation_type = "indel"
indels.ranges = rowRanges(indels)
indels.info$num_Germline_ALT_reads = as.character(unlist(indels.ranges$ALT))
indels.info$num_Germline_REF_reads = as.character(indels.ranges$REF)
indels.info$sample=paste0("sample","_",names(rowRanges(indels)))

mutations.info = rbind(snp.info, indels.info)

snpMutations.df = data.frame(chr = seqnames(rowRanges(snps)), pos=start(rowRanges(snps)), ref = mcols(rowRanges(snps))$REF, alt = unlist(mcols(rowRanges(snps))$ALT), sample=paste0("sample","_",names(rowRanges(snps))) )
indelMutations.df = data.frame(chr = seqnames(rowRanges(indels)), pos=start(rowRanges(indels)), ref = mcols(rowRanges(indels))$REF, alt = unlist(mcols(rowRanges(indels))$ALT), sample=paste0("sample","_",names(rowRanges(indels))) )
mutations.df=rbind(snpMutations.df,indelMutations.df)
write_tsv(mutations.df, "mutations.tsv")

 
### Request information from the Cancer Genome Interpreter. 
source_python("cgi.py")

## merge mutation analysis from CGI with the information data frame. 
CGI = read_tsv("mutation_analysis.tsv") %>%
    select(one_of(c("sample", "gene", "cdna" ,"protein", "consequence", "cadd_phred", "driver_gene_source", "driver_statement"))) %>%
    rename_at(vars(colnames(.)), ~ c("sample", "CGI_gene","CGI_cDNA","CGI_protein", "CGI_consequences", "CGI_CADD_Phed", "CGI_annotation", "CGI_driver")) %>%
    mutate_all(as.character)


merged = merge(as.data.frame(mutations.info), CGI)
#info.df
#geno(vcf)
