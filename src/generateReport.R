# change WD to dir containing som and gl mutcp format files
# put into the same Dir as these scripts several files
### CGI "drug_prescription.tsv" output file
### The output of the script "Merge_Auck.data_CGI.data.R" (file named "*.WITH_CGI_NO.SUP.somatic_indel_OR_snv.txt")
### The germline output file 

##### Put into the Images directory in same folder as this script and other input fikes the pdf B allele freq images from step 1, as well as rteh default image


library(knitr)
library(markdown)
library(rmarkdown)
library("here")
suppressMessages(library(dplyr))
suppressMessages(require(ggplot2))

# generate folder named images in same dir as this script and the input files
# Into this place the B allele images (with names such as "B allele frequency for tumour varscan-P1025B_m.snp.vcf.gz .pdf")
# then run the following bash 1-liners to reduce filenamnes to usuable names, after CD-ing into this dir
# for file in *; do mv "${file}" "${file/.snp.vcf.gz /}"; done
# for file in *; do mv "${file}" "${file/B allele frequency for tumour varscan-/}"; done
# leave default image in this dir
#SomFilePath <- list.files(path = "../data/", pattern = "WITH_CGI_NO.SUP.somatic_indel_OR_snv.txt") 
SomFilePath="../data/2018_WITH_CGI_NO.SUP.somatic_indel_OR_snv.txt"
print(SomFilePath)
S <- read.table(SomFilePath, skip = 0, sep = "\t", as.is = TRUE); colnames(S)<-S[1,]; S<-S[-1,]

#GLFilePath <- list.files(path = "../data/", pattern = ".germline_indel_OR_snv.txt") 
GLFilePath="../data/2018_germline_indel_OR_snv.txt"
G <- read.table(GLFilePath, skip = 0, sep = "\t", as.is = TRUE); colnames(G)<-G[1,]; G<-G[-1,]
print(GLFilePath)

SCGIFilePath <- list.files(path = "../data/", pattern = "drug_prescription.tsv") 
print(SCGIFilePath)
SCGIFilePath="../data/drug_prescription.tsv"
SCGI <- read.table(SCGIFilePath, skip = 0, sep = "\t", as.is = TRUE); colnames(SCGI)<-SCGI[1,]; SCGI <-SCGI[-1,]

e1<-lapply(strsplit(SCGI$SAMPLE,"varscan-"), `[`, 2); e1<-lapply(strsplit(as.character(e1),".snp"), `[`, 1)
e2<-lapply(strsplit(SCGI$SAMPLE,"chr"), `[`, 2)
SCGIb<-cbind(SCGI,unlist(e1),unlist(e2))
colnames(SCGIb)[15:16]<-c("tumour", "position")
SCGIc<-SCGIb[,c(15,16,6,7,8,13:16)]
SCGIc<-SCGIc[!SCGIc$EVIDENCE %in% c("Pre-clinical", "Early trials"), ]
SCGIc<-SCGIc[grep("amplification",SCGIc$BIOMARKER, invert = T),]
SCGIc<-SCGIc[grep("fusion",SCGIc$BIOMARKER, invert = T),]
SCGIc<-SCGIc[grep("deletion",SCGIc$BIOMARKER, invert = T),]
SCGIc<-SCGIc[grep("overexpression",SCGIc$BIOMARKER, invert = T),]

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
GTGL<-as.numeric(G$num_Germline_REF_reads)+as.numeric(G$num_Germline_ALT_reads)
GGLVarDepth<-as.numeric(G$num_Germline_ALT_reads)
GGnEx<-G$gnomAD_exome_All; GGnEx[GGnEx=="."]<-"0"; GGnEx<-as.numeric(GGnEx)

# Filter
G.filtered<-G[(GTGL >= GNdepth & GGLVarDepth >= GN.var.depth & GGnEx <= GgnomAD_exome_All),]
G.filtered<-G.filtered[! G.filtered$Mutation_type %in% GMutation_type_Exclude,]
G.filtered<-G.filtered[! G.filtered $Variant_supported_by_supplementary_reads == Variant_supported_by_supplementary_reads_Exclude,]
f1<-lapply(strsplit(G.filtered$tumour,"varscan-"), `[`, 2); f1<-lapply(strsplit(as.character(f1),".snp"), `[`, 1)
G.filtered<-cbind(G.filtered, G.filtered$Position, as.character(f1)); colnames(G.filtered)[33:34]<-c("IGV", "sample")
G.filtered$IGV<-lapply(strsplit(as.character(G.filtered$IGV),"_", fixed = TRUE), `[`, 1)
colnames(G.filtered)[24]<-"dbSNP147"

##### filter out on Som ##### 

# Define filters
Tdepth <- 20 # T overall depth 
Ndepth <- 20 # N overall depth 
T.var.depth <- 8 # T ALT supporting depth 
T.var.per <- 10 # T ALT supporting % reads
Somatic_P_value<-0.01
gnomAD_exome_All <- 0.01
Mutation_type_Exclude <-c("synonymous_SNV", "unknown")
Variant_supported_by_supplementary_reads_Exclude <-"YES"  

# generate table for display of som filters
st1<-c("Min normal tissue read depth", "Min tumour read depth", "Min tumour ALT base read depth", "Min tumour alt base % reads", "Max Varscan2 som P val", "Max variant pop frequency gnomAD exome all", "Exclude variants supported by sup reads")
st2<-c(Tdepth, Ndepth, T.var.depth, T.var.per, Somatic_P_value, gnomAD_exome_All, Variant_supported_by_supplementary_reads_Exclude)
SomRulesText1<-data.frame(cbind(st1,st2)); colnames(SomRulesText1)<-c("Feature", "Rule")
SomRulesText2<- paste("Excluding", Mutation_type_Exclude, "mutations")

# Modify matrix to facilitate analysis
S$Mutation.type[S$Mutation.type =="."]<-"splicing"
S$Indel.or.SNP[S$Indel.or.SNP =="SNP"]<-"SNV"

# Generate 5 vectors to use while filtering, which operate only on parent matrix
TGL<-as.numeric(S$num.Germline.REF.reads)+as.numeric(S$num.Germline.ALT.reads)
TSom<-as.numeric(S$num.Somatic.REF.reads)+as.numeric(S$num.Somatic.ALT.reads)
SomVarDepth<-as.numeric(S$num.Somatic.ALT.reads)
SomVarDepthPerCent<-as.numeric(S$num.Somatic.ALT.reads)/TSom*100
GnEx<-S$gnomAD.exome.All; GnEx[GnEx=="."]<-"0"; GnEx<-as.numeric(GnEx)

# Filter
S.filtered<-S[(TGL >= Ndepth & TSom >= Tdepth & SomVarDepth >= T.var.depth & SomVarDepthPerCent >= T.var.per & GnEx <= gnomAD_exome_All),]
S.filtered<-S.filtered[! S.filtered$Mutation.type %in% Mutation_type_Exclude,]
S.filtered<-S.filtered[as.numeric(S.filtered$Somatic_P_value) <= Somatic_P_value,]
S.filtered<-S.filtered[! S.filtered$Variant.supported.by.supplementary.reads == Variant_supported_by_supplementary_reads_Exclude,]

#for (count in 1:10){ # for debugging
for (count in 1:length(tumour.list)){
TT<-tumour.list[count]

#for aneuploidy images
TT.image<-paste0("images/",TT,".pdf")

TTT.image<-paste0(TT,".pdf")
if(TTT.image %in% list.files(path = "images")){
TT.image<-paste0("images/",TT,".pdf")
} else {TT.image<-"images/default.pdf"}

file.name <- paste0(TT,"_",Sys.Date(),".html")
Text1a <- paste("Report for Tumour",TT)
Text1b <- paste("Prepared on", weekdays(Sys.Date()), Sys.Date())
Text0 <- "Notes: In general, a CADD score >25 indicates likely protein functional disruption. Genome reference used = hg19. CGI annotation key; (i) cgc: gene included in the COSMIC Cancer Gene Census, (ii) biomarker: CGI Cancer Biomarkers database, (iii) intogen: gene whose mutations show positive selection across cancer cohorts, (iv) known: known in CGI database in specific tumour types."

if(nrow(S.filtered[S.filtered$sample==TT,])!=0){
S.filtered.temp<-S.filtered[S.filtered$sample==TT,]	
for.IGV<-paste(S.filtered.temp[,36], S.filtered.temp[,37], sep=":") # generate IGV-compatible chr and position
S.filtered.table<-cbind(for.IGV,S.filtered.temp[S.filtered.temp$sample==TT,c(5,44,8,43,16,20:23,47,46,45,40, 9,10,25)])
colnames(S.filtered.table)<-c("IGV","Gene","Consequence","Position","AA change", "Pop freq gnomAD exome","GL REF reads", "GL ALT reads","T REF reads", "T ALT reads", "CGI oncodriveMUT status", "CGI annotation","CADD score", "SIFT, MutationTaster, MetaLR, FATTHM","Variant position in ClinVar","Variant position in Cosmic v82", "Variant position in dbSNP147")
S.filtered.table<-S.filtered.table[order(S.filtered.table$Gene),]
} else {S.filtered.table<-"No Somatic SNVs or Indels Found for this tumour"}

# based on somatic filters, filter CGI drug report (SCGIc)
if(nrow(S.filtered[S.filtered$sample==TT,])!=0){
	SCGId<-SCGIc[as.character(SCGIc$tumour) ==TT,]
if(nrow(SCGId) !=0){
	
SCGId.table<-SCGId
} else {SCGId.table <-"No CGI recommendations for this tumour"}
}

if(nrow(G.filtered[G.filtered$sample==TT,])!=0){
G.filtered.table<-G.filtered[G.filtered$sample==TT,c(33,4,17,7,15,19:22,8, 11, 12, 29:32,24)]
colnames(G.filtered.table)<-c("ToIGV","Gene","Mutation type","Position","Pop freq gnomAD exome","GL REF reads", "GL ALT reads","T REF reads", "T ALT reads", "Clinvar", "Gene in Cosmic census germline database", "Gene in ACMG list?", "SIFT", "Mutation Taster", "MetaLR", "FATTHM", "dbSNP147")
} else {G.filtered.table<-"No Germline SNVs or Indels Found for this tumour"}

print(getwd())
rmarkdown::render('generate_report.Rmd', output_file = file.name, output_dir = getwd())
}

