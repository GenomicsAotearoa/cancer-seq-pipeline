## re-write of code written by Cris Print to identify mutations in a set of genes in a set of VCF files at low stringency

## Comments by Ben : ##


args <- commandArgs(trailingOnly = TRUE)
snpVCF=args[1]
indelVCF=args[2]
intermediateDirectory=args[3]
out1=args[4]
#out2=args[4]


#source("https://bioconductor.org/biocLite.R")

suppressMessages(library("VariantAnnotation"))
suppressMessages(library("tidyverse"))
suppressMessages(library("here"))
suppressMessages(library("httr"))
suppressMessages(library("jsonlite"))
suppressMessages(library("reticulate"))

#options(error=utils::recover) 
options(warn = -1)

## Data location
## This section finds the files to be reported up, checks that there is a snp and an indel file for each sample
## and an associated tbi file in the same directory. If successful, it puts the full file paths
## into a dataframe with one row per samples.  

## data location. Files should be paired bzip and corresponding tbi file with
## sample name at the beginning of the file name.

tumourName = snpVCF %>%
	str_split("/") %>%
	simplify() %>%
	last() %>%
	str_split("\\.") %>%
	simplify() %>%
	first()
 
#tumourName = str_match(snpVCF, "/(.*?)_")[,2]


metadata = read_tsv("../data/references/metadata.txt", col_types = cols())
mgl = metadata$Michelle.gene.list[!is.na(metadata$Michelle.gene.list)] # define Michelle's gene list
geneList = paste(mgl,collapse="|")

# prefilters, that are just greps on each line as they are read in.
regionFilter = function(x) grepl(".*", x, perl = TRUE)
#regionFilter = function(x) grepl("exonic|splicing|UTR5|upstream", x, perl = TRUE)

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


report.df = read_csv("../data/references/dataColumns.csv")
infoVariables = report.df %>%
    filter(location == "info")

snpFilter = filterVcf(snpVCF , "hg19", tempfile(), filters = FF, prefilters = PF )
indelFilter = filterVcf(indelVCF , "hg19" , tempfile(), filters = FF, prefilters = PF )

snps = readVcf(snpFilter)
indels = readVcf(indelFilter) 


snp.info = as(info(snps),"DataFrame")


#snp.info=as.tibble(snp.info, rownames=NULL)

snp.info = as.tibble(snp.info, rownames=NULL) %>%
    dplyr::select(one_of(infoVariables$variable))  %>%
    map(function(x) unlist(lapply(x, '[', 1))) %>%
    as.tibble() 

colnames(snp.info) =  infoVariables$name

snp.ranges = rowRanges(snps)

# This bit is a bunch of odds and ends, currently neccesary due to some weird/stupid/overly complex data
# structures in the R VCF data structures. i.e. not easily converted to data frame. Needs to be tidied up. 
# i.e. it should be generated from a config file, not hardcoded. 

#snp.info$num_Germline_ALT_reads = as.character(unlist(snp.ranges$ALT))
snp.info$num_germline_REF = geno(snps)$RD[, 1]
snp.info$num_germline_ALT = geno(snps)$AD[, 1]
snp.info$num_somatic_REF = geno(snps)$RD[, 2]
snp.info$num_somatic_ALT = geno(snps)$AD[, 2]
snp.info$Indel.or.SNP="snp"
snp.info$DownstreamGene=unlist(lapply(info(snps)$Gene.refGene, `[`, 2))
snp.info$position=names(rowRanges(snps))

# and some of the feild names need to be revised. Put an alt text over the column header in the report or something.
cosmicGenes = metadata$cosmic_10Nov17_somatic
snp.info$In.Cosmic.Census.Somatic.Nov17 = ifelse(snp.info$Main_gene %in% cosmicGenes | snp.info$DownstreamGene %in% cosmicGenes, "YES", "NO")
cosmicGenes = metadata$cosmic_10Nov17_germline
snp.info$In.Cosmic.Census.Somatic.Nov17 = ifelse(snp.info$Main_gene %in% cosmicGenes | snp.info$DownstreamGene %in% cosmicGenes, "YES", "NO")
acmgGenes = metadata$ACMG.genes
snp.info$Variant.gene.in.ACMG.57.gene.list = ifelse(snp.info$Main_gene %in% acmgGenes | snp.info$DownstreamGene %in% acmgGenes, "YES", "NO") 




indels.info = as(info(indels),"DataFrame")
indels.info=as.tibble(indels.info)
indels.info = as.tibble(indels.info, rownames=NULL) %>%
    dplyr::select(one_of(infoVariables$variable))  %>%
    map(function(x) unlist(lapply(x, '[', 1))) %>%
    as.tibble()
colnames(indels.info) =  infoVariables$name
#indels.info$Mutation_type = "indel"
indels.ranges = rowRanges(indels)


# This bit is a bunch of odds and ends, currently neccesary due to some weird/stupid/overly complex data
# structures in the R VCF data structures. i.e. not easily converted to data frame. Needs to be tidied up.

indels.info$num_germline_REF = geno(indels)$RD[, 1]
indels.info$num_germline_ALT = geno(indels)$AD[, 1]
indels.info$num_somatic_REF = geno(indels)$RD[, 2]
indels.info$num_somatic_ALT = geno(indels)$AD[, 2]
indels.info$Indel.or.SNP="indel"
indels.info$DownstreamGene=unlist(lapply(info(indels)$Gene.refGene, `[`, 2))
indels.info$position=names(rowRanges(indels))

cosmicGenes = metadata$cosmic_10Nov17_somatic
indels.info$In.Cosmic.Census.Somatic.Nov17 = ifelse(indels.info$Main_gene %in% cosmicGenes | indels.info$DownstreamGene %in% cosmicGenes, "YES", "NO")
cosmicGenes = metadata$cosmic_10Nov17_germline
indels.info$In.Cosmic.Census.Somatic.Nov17 = ifelse(indels.info$Main_gene %in% cosmicGenes | indels.info$DownstreamGene %in% cosmicGenes, "YES", "NO")
indels.info$Variant.gene.in.ACMG.57.gene.list = ifelse(indels.info$Main_gene %in% acmgGenes | indels.info$DownstreamGene %in% acmgGenes, "YES", "NO") 


#print(colnames(snp.info))
#print(colnames(indels.info))
 

mutations.info = rbind(snp.info, indels.info)
mutations.info$tumour = rep(tumourName,  times = nrow(mutations.info))
mutations.info$sample = rep(tumourName,  times = nrow(mutations.info))
mutations.info$chr = str_extract(mutations.info$position, "[^:]+")


snpMutations.df = data.frame(chr = seqnames(rowRanges(snps)), pos=start(rowRanges(snps)), ref = mcols(rowRanges(snps))$REF, alt = unlist(mcols(rowRanges(snps))$ALT), sample = names(rowRanges(snps)) )
indelMutations.df = data.frame(chr = seqnames(rowRanges(indels)), pos=start(rowRanges(indels)), ref = mcols(rowRanges(indels))$REF, alt = unlist(mcols(rowRanges(indels))$ALT),  sample = names(rowRanges(indels)) )
mutations.df = rbind(snpMutations.df,indelMutations.df)
write_tsv(mutations.df, "mutations.tsv")


### Request information from the Cancer Genome Interpreter. 
source_python("cgi.py")
getDrugData(tumourName)

source("/blue/project/cancer-seq-pipeline/src/plotSampleCRIS.R")

countFile = paste0("../data/intermediate/", tumourName, ".count")
counts = read_tsv(countFile, col_names = FALSE)

colnames(counts) = c("chr", "tileStart","tileEnd", "tumourName" , "tumourBase")

dnaObj = counts %>%
  select(one_of("chr","tileStart","tumourName", "tumourBase"))

#logratio = log2(dnaObj$P1003C+1) - log2(dnaObj$P1003A+1)
dnaObj = dnaObj %>%
	mutate(tumourName=tumourName +1) %>%
	mutate(tumourBase=tumourBase +1) 
logratio = log2(dnaObj$tumourName) - log2(dnaObj$tumourBase)
print(logratio)

	
	

CNA.object = CNA( genomdat = logratio, chrom = dnaObj$chr, maploc = dnaObj$tileStart, data.type = 'logratio')
smoothed.CNA.object <- smooth.CNA(CNA.object, smooth.region = 100000, trim=0.25)
CNA.object <- segment(smoothed.CNA.object, verbose = 1)
plotSampleCP(CNA.object, "/blue/project/cancer-seq-pipeline/data/intermediate/patient01B.pdf")


## merge mutation analysis from CGI with the information data frame. 


cgiFile=paste0(intermediateDirectory,"/", tumourName, "/mutation_analysis.tsv")
print(cgiFile)
CGI = read_tsv(cgiFile) %>%
    dplyr::select(one_of(c("sample", "gene", "cdna" ,"protein", "consequence", "cadd_phred", "driver_gene_source", "driver_statement"))) %>%
    rename_at(vars(colnames(.)), ~ c("position", "CGI_gene","CGI_cDNA","CGI_protein", "CGI_consequences", "CGI_CADD_Phed", "CGI_annotation", "CGI_driver")) %>%
    mutate_all(as.character)



merged = merge(as.data.frame(mutations.info), CGI)

print(out1)
write_tsv(merged, out1)
#info.df
#geno(vcf)
