# 7 July 2018
# Cris Print
# script to identify mutations in a set of genes in a set of VCF files at low stringency
# set working directory to directory containing multiple VCF files (both SNVs and Indels mixed together are OK)
#source("plotSampleCRIS.R")

## Modifed by Ben Curran.
#source("https://bioconductor.org/biocLite.R")

suppressMessages(library(VariantAnnotation))
suppressMessages(library(DNAcopy))
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(VariantAnnotation))
#suppressMessages(library(TxDb.Hsapiens.UCSC.hg19.knownGene))


## read in metadata/list of genes
metadata = read_tsv("data/metadata.txt")
mgl = metadata$Michelle.gene.list[!is.na(metadata$Michelle.gene.list)] # define Michelle's gene list
geneList = paste(mgl,collapse="|")
## location of directory containing files to be preocessed to generate reports. 
## files must be in bgzip format, co-located with a tbi file of the same name. 
## i.e. sample1.vcf.gz, with sample1.vcf.gz.tbi in the same directory.
## This should have a test to verify that some files exist/are found and exit gracefully if not.
raw.snp.vcf.files = list.files(path = "data/toBeProcessed", full.names = TRUE, pattern = ".vcf.gz$") # read in file list

## Haven't figured out what this is for yet. 
CV.tot <- matrix(data = NA, nrow = 0, ncol = 31) # dummy matrix for output
 


## Set the filters
# prefilters,raw greps on each line

regionFilter = function(x) grepl("\bexonic|splicing|UTR5|upstream\b", x, perl = TRUE)
geneFilter = function(x) grepl(paste0("\b",geneList,"\b"), x, perl = TRUE)

geneFilter = function(x) {
  thing=unlist(info(x)$Gene.refGene %in% mgl)
  print(unlist(info(x)$Gene.refGene))
  thing
}
# filters on the info segments
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
FF = FilterRules(list(infoFilter, genoFilter, geneFilter))

## work through files. 
for (file in raw.snp.vcf.files) {
  filt2 <- filterVcf(file, "hg19", tempfile(), index = FALSE, filters = FF, prefilters = PF)
  vcf = readVcf(filt2)
}
unlist(info(vcf)$Gene.refGene)
mgl
vcf
geneList
info(vcf)
testString
testString=unlist(info(vcf)$Gene.refGene)

mgl[mgl %in% unlist(info(vcf)$Gene.refGene)]
# read in and filter contents of the VCF files
for (i in 1:length(raw.snp.vcf.files)) { # iterate through VCFs
    #print(paste(i, "    ", raw.snp.vcf.files[i], "    Size =", prettyNum(file.size(raw.snp.vcf.files[i]),	big.mark = ",", scientific = FALSE), "bytes"))
    fl="/home/ben/workspace/uni/prosper/data/toBeProcessed/varscan-P1003C.snp.vcf.gz"
    o.raw.snp <- readVcf(fl, "hg19") # read in VCF file to start
    raw.snp <- o.raw.snp[(info(o.raw.snp)$SS == 1)] # germline calls
    raw.snp <- raw.snp[unlist(info(raw.snp)$Func.refGene == "exonic") | unlist(info(raw.snp)$Func.refGene == "splicing") | unlist(info(raw.snp)$Func.refGene == "UTR5") | unlist(info(raw.snp)$Func.refGene == "upstream")]  # gene position filtering
    raw.snp <- raw.snp[unlist(info(raw.snp)$wgEncodeDacMapabilityConsensusExcludable == ".")] # discard variants in Encode DAC black list regions
    raw.snp <- raw.snp[(geno(raw.snp)$AD[,1] + geno(raw.snp)$RD[, 1]) >= Ndepth]
    raw.snp <- raw.snp[(geno(raw.snp)$AD[,1] >= N.var.depth)]

    ?stopifnot
  
#	# filter by gene names
    raw.snp <- raw.snp[lapply(info(raw.snp)$Gene.refGene, `[`, 1) %in% mgl | lapply(info(raw.snp)$Gene.refGene, `[`, 2) %in% mgl] 
    print(class(raw.snp))
##writeVcf(raw.snp, paste0("trimmed_",raw.snp.vcf.files[i]))
#	
#
#	# extract data for each variant
	num.germline.REF <- geno(raw.snp)$RD[, 1]
	num.germline.ALT <- geno(raw.snp)$AD[, 1]
	num.somatic.REF <- geno(raw.snp)$RD[, 2]
	num.somatic.ALT <- geno(raw.snp)$AD[, 2]
	fract.somatic.ALT <- round(num.somatic.ALT/(num.somatic.ALT + num.somatic.REF),2)
#
	pre.SPV <- info(raw.snp)$SPV
	SPV <- round(pre.SPV, 5)
#
	supl <- info(raw.snp)$SUPPLEMENTARY
	if (length(supl) == 1) {
		suppl <- "TRUE"
	} else {
		suppl <- "FALSE"
	}
	supl <- gsub(FALSE, ".", supl)
	supl <- gsub(TRUE, "YES", supl)
	positions <- unlist(rowRanges(raw.snp)@ranges@NAMES)
	snp <- lapply(info(raw.snp)$avsnp147, `[`, 1)
	# D: Deleterious; T: Tolerated, see http://annovar.openbioinformatics.org/en/latest/user-guide/filter/#overview
	MetaLR <- info(raw.snp)$MetaLR_pred == "D"
	MetaLR <- gsub(FALSE, ".", MetaLR)
	MetaLR <- gsub(TRUE, "YES", MetaLR)
	# D: Deleterious (sift<=0.05); T: tolerated (sift>0.05)
	sift <- info(raw.snp)$SIFT_pred == "D"
	sift <- gsub(FALSE, ".", sift)
	sift <- gsub(TRUE, "YES", sift)
	# D: Deleterious; T: Tolerated, very similar to SIFT and Polyphen
	FATHMM <- info(raw.snp)$FATHMM_pred == 	"D"
	FATHMM <- gsub(FALSE, ".", FATHMM)
	FATHMM <- gsub(TRUE, "YES", FATHMM)
	#very similar to SIFT and PolyPhen.
	MutationTaster <- info(raw.snp)$MutationTaster_pred == "A" | info(raw.snp)$MutationTaster_pred == "D"
	MutationTaster <- gsub(FALSE, ".", 	MutationTaster)
	MutationTaster <- gsub(TRUE, "YES", MutationTaster)
	repeat_mask <- lapply(info(raw.snp)$rmsk, `[`, 1)
	ENCODE.Dac.unmappable <- lapply(info(raw.snp)$wgEncodeDacMapabilityConsensusExcludable, `[`, 1)
	segdup <- lapply(info(raw.snp)$genomicSuperDups, `[`, 1) # Duplications of >1000 Bases of Non-RepeatMasked Sequence
	ACMG <- lapply(info(raw.snp)$Gene.refGene,`[`, 1) %in% metadata$ACMG.genes | lapply(info(raw.snp)$Gene.refGene, 	`[`, 2) %in% metadata$ACMG.genes
	ACMG <- gsub(FALSE, ".", ACMG)
	ACMG <- gsub(TRUE, "YES", ACMG)
	# The popfreq_max database contains the maximum allele frequency from several population frequency databases, including 1000 Genomes Project (ALL+5 ethnicity groups), ESP6500 (ALL+2 ethnicity groups), ExAC (ALL+7 ethnicity groups), CG46.
	popfreqmax <- lapply(info(raw.snp)$PopFreqMax, 	`[`, 1)
	muttype <- lapply(info(raw.snp)$ExonicFunc.refGene,	`[`, 1)
	gene.position <- lapply(info(raw.snp)$Func.refGene,	`[`, 1)
	clinvar <- as.matrix(lapply(info(raw.snp)$CLINSIG, 	`[`, 1))
	cosmic <- lapply(info(raw.snp)$cosmic82, 	`[`, 1)
	ExAc_ALL <- lapply(info(raw.snp)$ExAC_ALL,`[`, 1)
	gnomAD_genome_ALL <- lapply(info(raw.snp)$gnomAD_genome_ALL,`[`, 1)
	gnomAD_exome_ALL <- lapply(info(raw.snp)$gnomAD_exome_ALL, 	`[`, 1)

	CosCensusSomatic <- lapply(info(raw.snp)$Gene.refGene, `[`, 1) %in% metadata$cosmic_10Nov17_somatic | lapply(info(raw.snp)$Gene.refGene, `[`, 2) %in% metadata$cosmic_10Nov17_somatic
	CosCensusSomatic <- gsub(FALSE, ".", 	CosCensusSomatic)
	CosCensusSomatic <- gsub(TRUE, "YES",	CosCensusSomatic)

	CosCensusGermline <- lapply(info(raw.snp)$Gene.refGene,	`[`, 1) %in% metadata$cosmic_10Nov17_germline | lapply(info(raw.snp)$Gene.refGene, `[`, 2) %in% metadata$cosmic_10Nov17_germline
	CosCensusGermline <- gsub(FALSE, ".", CosCensusGermline)
	CosCensusGermline <- gsub(TRUE, "YES",CosCensusGermline)

	if (length(grep("indel", raw.snp.vcf.files[i])) ==	1) {
		ii <- "INDEL"
	} else {
		ii <- "SNP"
	} # are the variants in this VCF indels or snps?

	if (length(lapply(info(raw.snp)$Gene.refGene, `[`, 2) == 1)) {
        ij <- lapply(info(raw.snp)$Gene.refGene, `[`, 2)
	} else {
		ij <- "NA"
	} # 2nd gene

	# for this VCF file, combine extratced information and name columns
	if (nrow(raw.snp) > 0) {

		CV.out <- cbind(rep(raw.snp.vcf.files[i], nrow(raw.snp)), rep(ii, nrow(raw.snp)), lapply(info(raw.snp)$Gene.refGene, `[`, 1), ij, SPV, positions, 
			clinvar, cosmic, CosCensusSomatic, 
			CosCensusGermline, ACMG, 
			popfreqmax, ExAc_ALL, gnomAD_exome_ALL, 
			gnomAD_genome_ALL, muttype, 
			gene.position, num.germline.REF, 
			num.germline.ALT, num.somatic.REF, 
			num.somatic.ALT, fract.somatic.ALT, 
			snp, repeat_mask, ENCODE.Dac.unmappable, 
			segdup, supl, sift, MutationTaster, 
			MetaLR, FATHMM)

		colnames(CV.out) <- c("tumour", 
			"Indel_or_SNP", "Main_gene", 
			"Downstream_gene", "Somatic_P_value", 
			"Position", "Variant_position_in_ClinVar", 
			"Variant_position_in_Cosmic_v82", 
			"In_Cosmic_Census_Somatic_Nov17", 
			"In_Cosmic_Census_Germline_Nov17", 
			"Variant_gene_in_ACMG_57_gene_list?", 
			"Maximum_Population_Frequency", 
			"ExAc_All", "gnomAD_exome_All", 
    		"gnomAD_genome All", "Mutation_type", 
			"ANNOVAR_gene_position", 
			"num_Germline_REF_reads", 
			"num_Germline_ALT_reads", 
			"num_Somatic_REF_reads", 
			"num_Somatic_ALT_reads", 
			"Fraction_of_somatic_reads_that_are_ALT", 
			"Variant_position_in_dbSNP147 ?", 
			"Variant_position_in_repeat_masked_region?", 
			"Variant_position_in_ENCODE_Dac_unmappable_region?", 
			"Variant_position_in_regions_of_segmental_duplication?", 
			"Variant_supported_by_supplementary_reads?", 
			"Variant_effect_SIFT", 
			"Variant_effect_Mutation_Taster", 
			"Variant_effect_MetaLR", 
			"Variant_effect_FATTHM")

		CV.tot <- rbind(CV.tot, CV.out)
	}
	
}

rownames(CV.tot) <- NULL
timestamp=format(Sys.time(), "%d%m%Y-%H:%M")
file <- paste0("../results/",timestamp, "germline_indel_OR_snv.txt")
write.table(t(colnames(CV.tot)), file = file, 	sep = "\t", row.names = "index", col.names = FALSE, append = TRUE)
write.table(CV.tot, file = file, sep = "\t", row.names = TRUE, col.names = FALSE, 	append = TRUE)
