## re-write of code written by Cris Print to identify mutations in a set of genes in a set of VCF files at low stringency

## Comments by Ben : ##


args <- commandArgs(trailingOnly = TRUE)
snpsCF=args[1]
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
here()

dataLocation = "data.raw"
#files = list.files(dataLocation, pattern="*.gz$")
#files = list(rawVCFfile)
#sampleNames = unlist(lapply(files, function(x) strsplit(x, "\\.")[[1]][1]))
#sampleNames =(unique(sampleNames))
#sampleFiles = list.files(dataLocation, pattern=paste0(sampleNames[1],".*anno.*.gz$"))
#tbiFile=paste0(sampleFiles[1],".tbi")

#files.df = data.frame(row.names=sampleNames, snps = character(length(sampleNames)),indels = character(length(sampleNames)), stringsAsFactors=FALSE)

#for (sample in sampleNames){
    #snpFile = paste0(dataLocation, list.files(dataLocation, pattern=paste0(sample,".snp.*anno.*.gz$")))
    #indelFile = paste0(dataLocation, list.files(dataLocation, pattern=paste0(sample,".indel.*anno.*.gz$")))
#    snpFile = files[1] 
    #tbiFile=paste0(snpFile,".tbi")
    
    #if (!file_test("-f", snpFile)){
    #    print(paste0("No snp vcf file for ", sample, " cannot be found in the data directory"))
    #}#else if (!file_test("-f", tbiFile)){
     #   print(paste0("The corresponding tbi file for ", sample, " snps, cannot be found in the data directory"))
    #}else{
    #    files.df[sample,]$snps = snpFile
    #}
       
    #tbiFile=paste0(indelFile,".tbi") 
    #if (!file_test("-f", indelFile)){
    #    print(paste0("No indel vcf file for ", sample, " cannot be found in the data directory"))
    #}else if (!file_test("-f", tbiFile)){
    #    print(paste0("The corresponding tbi file for ", sample, " indels, cannot be found in the data directory"))
    #}else{
    #    files.df[sample,]$indels = snpFile
    #}
#}
#print(rownames(files.df))
## Filters
## This section sets up pre filters and filters to isolate the desired variant calls


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

#report.df$variable
#print(infoVariables)
#genoVariables = 
#thing=mcols(rowRanges(vcf))
#thing
#thing=as.tibble(rowRanges(vcf))
#thing$ALT[1]$seq
#class(thing$ALT[1])

#report.df$location

for(i in 1:nrow(files.df)) {
    row = files.df[i,]
    filt2 <- filterVcf(row$snps, "hg19", tempfile(), filters = FF, prefilters = PF )
    vcf = readVcf(filt2)
    info.df = as(info(vcf), "DataFrame")
    info.df=as.tibble(info.df)
    info.df = as.tibble(info.df, rownames=NULL) %>%
         dplyr::select(one_of(infoVariables$variable))  %>%
         map(function(x) unlist(lapply(x, '[', 1))) %>%
         as.tibble() 
    colnames(info.df) =  infoVariables$name
    info.df$Mutation_type = "snp"
    
    ranges = rowRanges(vcf)
    info.df$num_Germline_ALT_reads = as.character(unlist(ranges$ALT))
    info.df$num_Germline_REF_reads = as.character(ranges$REF)

    
    #geno.df = as.DataFrame(geno(vcf), row.names=NULL)
    #geno.df = as.tibble(geno.df) %>%
    #     select(one_of(infoVariables$variable))  %>%
    #     map(function(x) unlist(lapply(x, '[', 1))) %>%
    #     as.tibble()

}

ranges[1,1:10]
mcols(ranges)["sampleNames"]

#info.df$sample=paste0(rownames(files.df[i,]),"_",names(rowRanges(vcf)))
#print(head(info.df))

range=rowRanges(vcf)


mutations.df = data.frame(chr = seqnames(rowRanges(vcf)), pos=start(rowRanges(vcf)), ref = mcols(rowRanges(vcf))$REF, alt = unlist(mcols(rowRanges(vcf))$ALT), sample=paste0(rownames(files.df[i,]),"_",names(rowRanges(vcf))) )
#write_tsv(mutations.df, "mutations.tsv")

 
### Request information from the Cancer Genome Interpreter. 
source_python("/home/ben/workspace/prosper/src/cgi.py")

## merge mutation analysis from CGI with the information data frame. 
CGI = read_tsv("mutation_analysis.tsv") %>%
    select(one_of(c("sample", "gene", "cdna" ,"protein", "consequence", "cadd_phred", "driver_gene_source", "driver_statement"))) %>%
    rename_at(vars(colnames(.)), ~ c("sample", "CGI_gene","CGI_cDNA","CGI_protein", "CGI_consequences", "CGI_CADD_Phed", "CGI_annotation", "CGI_driver")) %>%
    mutate_all(as.character)


merged = merge(as.data.frame(info.df), CGI)
#info.df
#geno(vcf)
