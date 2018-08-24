# script to merge in CGI outpuyt to Auckland variant data
library(here)

#print("CHOOSE Variant data file (of form *.NO.SUP.somatic_indel_OR_snv.txt)")
#FileChoice<-file.choose() # choose input file
file="../data/2018_NO.SUP.somatic_indel_OR_snv.txt"
Auck <- read.delim(file,header=T,sep="\t") # read in data but not header row


#print("CHOOSE CGI mutation output file (named mutation_analysis.tsv)")
#FileChoice<-file.choose() # choose input file
file1="../data/mutation_analysis.tsv"
CGI <- read.delim(file1,header=T,sep="\t") # read in data but not header row
CGI<-cbind(as.character(CGI$sample), as.character(CGI$gene), as.character(CGI$cdna), as.character(CGI$protein), as.character(CGI$consequence), as.character(CGI$cadd_phred), as.character(CGI$driver_gene_source), as.character(CGI$driver_statement))

var.CGI<- merge(Auck, CGI, all = TRUE, by.x=35, by.y=1)
colnames(var.CGI)[41:47]<-c("CGI_gene","CGI_cDNA","CGI_protein", "CGI_consequences", "CGI_CADD_Phed", "CGI_annotation", "CGI_driver")


# write file 
rownames(var.CGI) <- NULL
file <- paste(date(), "WITH_CGI_NO.SUP.somatic_indel_OR_snv.txt",sep = "_")
write.table(t(colnames(var.CGI)), file = file, sep = "\t", row.names = F, col.names = FALSE, append = TRUE)
write.table(var.CGI, file = file, sep = "\t", row.names = F, col.names = FALSE, append = TRUE)
	
	#CGI explanations
	#The gene is reported to be driver according to:
	#•	cgc: gene included in the COSMIC Cancer Gene Census
	#•	biomarkers: gene which shape the response to a given anti-cancer drug upon alterations (see the CGI Cancer Biomarkers database)
	#•	intogen: gene whose mutations show signals of positive selection across cancer cohorts (see IntoGen)
