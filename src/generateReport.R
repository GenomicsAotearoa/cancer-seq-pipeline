# change WD to dir containing som and gl mutcp format files
# put into the same Dir as these scripts several files
### CGI "drug_prescription.tsv" output file
### The output of the script "Merge_Auck.data_CGI.data.R" (file named "*.WITH_CGI_NO.SUP.somatic_indel_OR_snv.txt")
### The germline output file 
print("********************************wibble*******************************************")

args <- commandArgs(trailingOnly = TRUE)
SomFilePath=args[1]
GLFilePath=args[2]
baseDirectory=args[3]
out1=args[4]
out2=args[5]

##### Put into the Images directory in same folder as this script and other input fikes the pdf B allele freq images from step 1, as well as rteh default image
#SomFilePath = "intermediateFiles/P1014A_1.fastq.align.removeDuplicates.removeSuppplementary.pileUp.annotation.vcf.count"
#GLFilePath = "intermediateFiles/P1014A_1.fastq.align.removeDuplicates.removeSuppplementary.pileUp.annotation.vcf.2.count"

library(knitr)
library(markdown)
library(rmarkdown)
library("here")
library(tidyverse)
suppressMessages(library(dplyr))
suppressMessages(require(ggplot2))

# generate folder named images in same dir as this script and the input files
# Into this place the B allele images (with names such as "B allele frequency for tumour varscan-P1025B_m.snp.vcf.gz .pdf")
# then run the following bash 1-liners to reduce filenamnes to usuable names, after CD-ing into this dir
# for file in *; do mv "${file}" "${file/.snp.vcf.gz /}"; done
# for file in *; do mv "${file}" "${file/B allele frequency for tumour varscan-/}"; done
# leave default image in this dir
#SomFilePath <- list.files(path = "../data/", pattern = "WITH_CGI_NO.SUP.somatic_indel_OR_snv.txt") 
#SomFilePath="../data/2018_WITH_CGI_NO.SUP.somatic_indel_OR_snv.txt"
#S <- read.table(SomFilePath, skip = 0, sep = "\t", as.is = TRUE); colnames(S)<-S[1,]; S<-S[-1,]
S = read_tsv(SomFilePath)

#GLFilePath <- list.files(path = "../data/", pattern = ".germline_indel_OR_snv.txt") 
#GLFilePath="../data/2018_germline_indel_OR_snv.txt"
#G <- read.table(GLFilePath, skip = 0, sep = "\t", as.is = TRUE); colnames(G)<-G[1,]; G<-G[-1,]
G = read_tsv(GLFilePath)

#SCGIFilePath <- list.files(path = ".", pattern = "drug_prescription.tsv") 
#SCGIFilePath="drug_prescription.tsv"
#SCGI <- read.table(SCGIFilePath, skip = 0, sep = "\t", as.is = TRUE); colnames(SCGI)<-SCGI[1,]; SCGI <-SCGI[-1,]

#e1<-lapply(strsplit(SCGI$SAMPLE,"varscan-"), `[`, 2); e1<-lapply(strsplit(as.character(e1),".snp"), `[`, 1)
#e2<-lapply(strsplit(SCGI$SAMPLE,"chr"), `[`, 2)
#SCGIb<-cbind(SCGI,unlist(e1),unlist(e2))
#colnames(SCGIb)[15:16]<-c("tumour", "position")
#SCGIc<-SCGIb[,c(15,16,6,7,8,13:16)]
#SCGIc<-SCGIc[!SCGIc$EVIDENCE %in% c("Pre-clinical", "Early trials"), ]
#SCGIc<-SCGIc[grep("amplification",SCGIc$BIOMARKER, invert = T),]
#SCGIc<-SCGIc[grep("fusion",SCGIc$BIOMARKER, invert = T),]
#SCGIc<-SCGIc[grep("deletion",SCGIc$BIOMARKER, invert = T),]
#SCGIc<-SCGIc[grep("overexpression",SCGIc$BIOMARKER, invert = T),]

tumour.list <- unique(S$sample)
Text2 <- "Somatic Variants"
Text3 <- "Germline Variants"

##### filter out on GL ##### 

# Define filters
GNdepth <- 10 # N overall depth 
GN.var.depth <- 4 # T ALT supporting depth 
GgnomAD_exome_All <- 0.01
GMutation_type_Exclude <-c("synonymous_SNV", "unknown")
Variant_supported_by_supplementary_reads_Exclude <-"YES" 



# generate table for display of GL filters
st1<-c("Min normal tissue read depth", "Min normal tissue ALT base read depth", "Max var frequency gnomAD exome", "Exclude variants supported by sup reads")
st2<-c(GNdepth, GN.var.depth, GgnomAD_exome_All, Variant_supported_by_supplementary_reads_Exclude)
GLRulesText1<-data.frame(cbind(st1,st2)); colnames(GLRulesText1)<-c("Feature", "Rule")
GLRulesText2<- paste("Excluding", GMutation_type_Exclude, "mutations")

# Modify matrix to facilitate analysis
G$Mutation_type[G$Mutation_type =="."]<-"splicing"
G$Variant_position_in_ClinVar[G$Variant_position_in_ClinVar =="."]<-"No record"
# Generate 5 vectors to use while filtering, which operate only on parent matrix
GTGL<-as.numeric(G$num_germline_REF)+as.numeric(G$num_germline_ALT)
GGLVarDepth<-as.numeric(G$num_germline_ALT)
GGnEx<-G$gnomAD_exome_All; GGnEx[GGnEx=="."]<-"0"; GGnEx<-as.numeric(GGnEx)

# Filter
G.filtered<-G[(GTGL >= GNdepth & GGLVarDepth >= GN.var.depth & GGnEx <= GgnomAD_exome_All),]
G.filtered<-G.filtered[! G.filtered$Mutation_type %in% GMutation_type_Exclude,]
#G.filtered<-G.filtered[! G.filtered $Variant_supported_by_supplementary_reads == Variant_supported_by_supplementary_reads_Exclude,]
#f1<-lapply(strsplit(G.filtered$sample,"sample_"), `[`, 2); f1<-lapply(strsplit(as.character(f1),".snp"), `[`, 1)

#G.filtered = G.filtered %>%
#	separate(sample, c("sample", "IGV"), sep="_")

#G.filtered<-cbind(G.filtered, G.filtered$Position, as.character(f1)); colnames(G.filtered)[33:34]<-c("IGV", "sample")

##all this needs to go bac  in
#G.filtered$IGV<-lapply(strsplit(as.character(G.filtered$IGV),"_", fixed = TRUE), `[`, 1)
#colnames(G.filtered)[24]<-"dbSNP147"


##### filter out on Som ##### 

# Define filters
Tdepth <- 20 # T overall depth 
Ndepth <- 20 # N overall depth 
T.var.depth <- 8 # T ALT supporting depth 
T.var.per <- 10 # T ALT supporting % reads
Somatic_P_value<-0.05
gnomAD_exome_All <- 0.01
Mutation_type_Exclude <-c("synonymous_SNV", "unknown")
Variant_supported_by_supplementary_reads_Exclude <-"YES"  


# generate table for display of som filters
st1<-c("Min normal tissue read depth", "Min tumour read depth", "Min tumour ALT base read depth", "Min tumour alt base % reads", "Max Varscan2 som P val", "Max variant pop frequency gnomAD exome all", "Exclude variants supported by sup reads")
st2<-c(Tdepth, Ndepth, T.var.depth, T.var.per, Somatic_P_value, gnomAD_exome_All, Variant_supported_by_supplementary_reads_Exclude)
SomRulesText1<-data.frame(cbind(st1,st2)); colnames(SomRulesText1)<-c("Feature", "Rule")
SomRulesText2<- paste("Excluding", Mutation_type_Exclude, "mutations")


# Modify matrix to facilitate analysis
S$Mutation_type[S$Mutation_type =="."]<-"splicing"
#S$Indel.or.SNP[S$Indel.or.SNP =="SNP"]<-"SNV"


# Generate 5 vectors to use while filtering, which operate only on parent matrix
TGL<-as.numeric(S$num_germline_REF)+as.numeric(S$num_germline_ALT)
TSom<-as.numeric(S$num_somatic_REF)+as.numeric(S$num_somatic_ALT)

SomVarDepth<-as.numeric(S$num_somatic_ALT)
SomVarDepthPerCent<-as.numeric(S$num_somatic_ALT)/TSom*100
GnEx<-S$gnomAD_exome_All; GnEx[GnEx=="."]<-"0"; GnEx<-as.numeric(GnEx)


# Filter
## The p.value is curretly commented out. It shouldn't be. Currently it's the straw that breaks the camels back and filters everything out. 
## And I may hav mentioned somewhere, this does not fall over in a particularly graceful manner.

S.filtered<-S[(TGL >= Ndepth & TSom >= Tdepth & SomVarDepth >= T.var.depth & SomVarDepthPerCent >= T.var.per & GnEx <= gnomAD_exome_All),]
S.filtered<-S.filtered[! S.filtered$Mutation_type %in% Mutation_type_Exclude,]
#S.filtered<-S.filtered[as.numeric(S.filtered$Somatic_P_value) <= Somatic_P_value,]
#S.filtered<-S.filtered[! S.filtered$Variant_supported_by_supplementary_reads == Variant_supported_by_supplementary_reads_Exclude,]

#for (count in 1:10){ # for debugging
for (count in 1:length(tumour.list)){
TT<-tumour.list[count]

#for aneuploidy images


TT.image<-paste0(baseDirectory,"/data/intermediate/patient01B.pdf")


file.name <- paste0(TT,"_",Sys.Date(),".html")
Text1a <- paste("Report for Tumour",TT)
Text1b <- paste("Prepared on", weekdays(Sys.Date()), Sys.Date())
Text0 <- "Notes: In general, a CADD score >25 indicates likely protein functional disruption. Genome reference used = hg19. CGI annotation key; (i) cgc: gene included in the COSMIC Cancer Gene Census, (ii) biomarker: CGI Cancer Biomarkers database, (iii) intogen: gene whose mutations show positive selection across cancer cohorts, (iv) known: known in CGI database in specific tumour types."

if(nrow(S.filtered[S.filtered$sample==TT,])!=0){

	S.filtered.temp<-S.filtered[S.filtered$sample==TT,]	
	for.IGV<-paste(S.filtered.temp[,36], S.filtered.temp[,37], sep=":") # generate IGV-compatible chr and position

#colnames(S.filtered.table)<-c("IGV","Gene","Consequence","Position","AA change", "Pop freq gnomAD exome","GL REF reads", "GL ALT reads","T REF reads", "T ALT reads", "CGI oncodriveMUT status", "CGI annotation","CADD score", "SIFT, MutationTaster, MetaLR, FATTHM","Variant position in ClinVar","Variant position in Cosmic v82", "Variant position in dbSNP147")



#Check these are the right feilds with Cris, because he's referenced them by position and as a consequence, I have no idea if I'm putting the right numbers in the right places.
	reportFile=paste0(baseDirectory,"/data/references/finalReportColumns.csv")

	report.df = read_csv(reportFile) %>%
		drop_na() 

	reportNames=c(report.df$name,"CGI_gene","CGI_cDNA","CGI_protein", "CGI_consequences", "CGI_CADD_Phed", "CGI_annotation", "CGI_driver")

	S.filtered.table = S.filtered.temp %>%
		#rename_all(tolower) %>%
		#rename_all(gsub, pattern = "\\.", replacement = "_") %>% # hack to get around the fact that half these column names have periods in them, half have underscores, whilst trying to transition everything to _
    		#dplyr::select(one_of(tolower(report.df$variable))) %>%
		#rename_at(vars(colnames(.)), ~ report.df$name) %>%
    		#rename_at(vars(colnames(.)), ~ c("Gene", "Position", "Pop freq gnomAD exome", "GL REF reads", "GL ALT reads" , "CGI oncodriveMut")) %>%
		rename(Gene = Main_gene) 
    		#mutate_all(as.character)

	#S.filtered.table = S.filtered.table[order(S.filtered.table$Gene),]
}
else {
	S.filtered.table<-"No Somatic SNVs or Indels Found for this tumour"
}

#"CGI_gene","CGI_cDNA","CGI_protein", "CGI_consequences", "CGI_CADD_Phed", "CGI_annotation", "CGI_driver"

# based on somatic filters, filter CGI drug report (SCGIc)
#if(nrow(S.filtered[S.filtered$sample==TT,])!=0){
#	SCGId<-SCGIc[as.character(SCGIc$tumour) ==TT,]
#	if(nrow(SCGId) !=0){
#		SCGId.table<-SCGId
#	} 
#	else {SCGId.table <-"No CGI recommendations for this tumour"}
#}

if(nrow(G.filtered[G.filtered$sample==TT,])!=0){
    G.filtered.table = G.filtered %>%
        rename(Gene = Main_gene)
}else {
    G.filtered.table<-"No Germline SNVs or Indels Found for this tumour"
}



markdownFile=paste0(baseDirectory, "/data/generate_report.Rmd")
rmarkdown::render(markdownFile, output_file =out1)
}

