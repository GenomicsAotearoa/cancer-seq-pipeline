# script to identify mutations in a set of genes in a set of VCF files at low stringency
# set working directory to directory containing multiple VCF files (both SNVs and Indels mixed together are OK)

# atfter this script upload the For CGI output file to CGI, using all tumours.

source("plotSampleCRIS.R")
library(VariantAnnotation)
library(DNAcopy)
library(stringr)


# user querry
print("Do you want to plot B allele freq plots? (y or n)")
BAF.plot<-readLines(con = stdin(), n = 1, ok = TRUE)

#choose filters
Tdepth <- 20 # T overall depth 
Ndepth <- 20 # N overall depth 
T.var.depth <- 8 # T ALT supporting depth 
T.var.per <- 10 # T ALT supporting % reads 

# put Cris' metadata file in with the VCFs (one column) tab delim header row
metadata <- read.table("metadata.txt", 
	skip = 0, sep = "\t", as.is = TRUE)
colnames(metadata) <- metadata[1, ]
metadata <- metadata[-1, ]

raw.snp.vcf.files <- list.files(path = ".", 
	pattern = ".vcf") # read in file list

CV.tot <- matrix(data = NA, nrow = 0, ncol = 31) # dummy matrix for output
time.taken<-system.time(


# read in and filter contents of the VCF files
for (i in 1:length(raw.snp.vcf.files)) { # iterate through VCFs
	print(paste(i, "    ", raw.snp.vcf.files[i], 
		"    Size =", prettyNum(file.size(raw.snp.vcf.files[i]), 
			big.mark = ",", scientific = FALSE), 
		"bytes"))
	o.raw.snp <- readVcf(raw.snp.vcf.files[i], 
		"hg19") # read in VCF file to start
	raw.snp <- o.raw.snp[(info(o.raw.snp)$SS == 2)] # remove all but putative somatic mutations and LoH
	raw.snp <- raw.snp[(geno(raw.snp)$AD[, 
		2] + geno(raw.snp)$RD[, 2]) >= 
		Tdepth]
	raw.snp <- raw.snp[(geno(raw.snp)$AD[, 
		1] + geno(raw.snp)$RD[, 1]) >= 
		Ndepth]
	raw.snp <- raw.snp[(geno(raw.snp)$AD[, 
		2] >= T.var.depth)]
	raw.snp <- raw.snp[(100 * geno(raw.snp)$AD[, 
		2]/(geno(raw.snp)$AD[, 2] + geno(raw.snp)$RD[, 
		2])) >= (T.var.per)]
	raw.snp <- raw.snp[(geno(raw.snp)$AD[, 
		1] == 0)] # no ALT reads in germline
	raw.snp <- raw.snp[unlist(info(raw.snp)$wgEncodeDacMapabilityConsensusExcludable == 
		".")]
	# then filter by gene names
	#raw.snp <- raw.snp[lapply(info(raw.snp)$Gene.refGene, `[`, 1) %in% c("MEN1","MUTYH","DAXX") | lapply(info(raw.snp)$Gene.refGene, `[`, 2) %in% c("MEN1","MUTYH","DAXX")] # don't use gene-specific filtering here

raw.snp <- raw.snp[unlist(info(raw.snp)$Func.refGene == "exonic") | unlist(info(raw.snp)$Func.refGene == "splicing") | unlist(info(raw.snp)$Func.refGene == "UTR5") | unlist(info(raw.snp)$Func.refGene == "upstream")] 
	
#raw.snp <- raw.snp[unlist(info(raw.snp)$ExonicFunc.refGene == "nonsynonymous_SNV") | unlist(info(raw.snp)  $ExonicFunc.refGene == "stopgain")] # gene position filtering
	
# write filtered VCF
writeVcf(raw.snp, paste0("trimmed_",raw.snp.vcf.files[i]))


	# extract data for each variant
	num.germline.REF <- geno(raw.snp)$RD[, 
		1]
	num.germline.ALT <- geno(raw.snp)$AD[, 
		1]
	num.somatic.REF <- geno(raw.snp)$RD[, 
		2]
	num.somatic.ALT <- geno(raw.snp)$AD[, 
		2]
	fract.somatic.ALT <- round(num.somatic.ALT/(num.somatic.ALT + 
		num.somatic.REF),2)

	pre.SPV <- info(raw.snp)$SPV
	SPV <- round(pre.SPV, 5)

	supl <- info(raw.snp)$SUPPLEMENTARY
	if (length(supl) == 1) {
		suppl <- "TRUE"
	} else {
		suppl <- "FALSE"
	}
	supl <- gsub(FALSE, ".", supl)
	supl <- gsub(TRUE, "YES", supl)
	positions <- unlist(rowRanges(raw.snp)@ranges@NAMES)
	snp <- lapply(info(raw.snp)$avsnp147, 
		`[`, 1)
	# D: Deleterious; T: Tolerated, see http://annovar.openbioinformatics.org/en/latest/user-guide/filter/#overview
	MetaLR <- info(raw.snp)$MetaLR_pred == 
		"D"
	MetaLR <- gsub(FALSE, ".", MetaLR)
	MetaLR <- gsub(TRUE, "YES", MetaLR)
	# D: Deleterious (sift<=0.05); T: tolerated (sift>0.05)
	sift <- info(raw.snp)$SIFT_pred == 
		"D"
	sift <- gsub(FALSE, ".", sift)
	sift <- gsub(TRUE, "YES", sift)
	# D: Deleterious; T: Tolerated, very similar to SIFT and Polyphen
	FATHMM <- info(raw.snp)$FATHMM_pred == 
		"D"
	FATHMM <- gsub(FALSE, ".", FATHMM)
	FATHMM <- gsub(TRUE, "YES", FATHMM)
	#very similar to SIFT and PolyPhen.
	MutationTaster <- info(raw.snp)$MutationTaster_pred == 
		"A" | info(raw.snp)$MutationTaster_pred == 
		"D"
	MutationTaster <- gsub(FALSE, ".", 
		MutationTaster)
	MutationTaster <- gsub(TRUE, "YES", 
		MutationTaster)
	repeat_mask <- lapply(info(raw.snp)$rmsk, 
		`[`, 1)
	ENCODE.Dac.unmappable <- lapply(info(raw.snp)$wgEncodeDacMapabilityConsensusExcludable, 
		`[`, 1)
	segdup <- lapply(info(raw.snp)$genomicSuperDups, 
		`[`, 1) # Duplications of >1000 Bases of Non-RepeatMasked Sequence
	ACMG <- lapply(info(raw.snp)$Gene.refGene, 
		`[`, 1) %in% metadata$ACMG.genes | 
		lapply(info(raw.snp)$Gene.refGene, 
			`[`, 2) %in% metadata$ACMG.genes
	ACMG <- gsub(FALSE, ".", ACMG)
	ACMG <- gsub(TRUE, "YES", ACMG)
	# The popfreq_max database contains the maximum allele frequency from several population frequency databases, including 1000 Genomes Project (ALL+5 ethnicity groups), ESP6500 (ALL+2 ethnicity groups), ExAC (ALL+7 ethnicity groups), CG46.
	popfreqmax <- lapply(info(raw.snp)$PopFreqMax, 
		`[`, 1)
	muttype <- lapply(info(raw.snp)$ExonicFunc.refGene, 
		`[`, 1)
	gene.position <- lapply(info(raw.snp)$Func.refGene, 
		`[`, 1)
	clinvar <- as.matrix(lapply(info(raw.snp)$CLINSIG, 
		`[`, 1))
	cosmic <- lapply(info(raw.snp)$cosmic82, 
		`[`, 1)
	ExAc_ALL <- lapply(info(raw.snp)$ExAC_ALL, 
		`[`, 1)
	gnomAD_genome_ALL <- lapply(info(raw.snp)$gnomAD_genome_ALL, 
		`[`, 1)
	gnomAD_exome_ALL <- lapply(info(raw.snp)$gnomAD_exome_ALL, 
		`[`, 1)

	CosCensusSomatic <- lapply(info(raw.snp)$Gene.refGene, 
		`[`, 1) %in% metadata$cosmic_10Nov17_somatic | 
		lapply(info(raw.snp)$Gene.refGene, 
			`[`, 2) %in% metadata$cosmic_10Nov17_somatic
	CosCensusSomatic <- gsub(FALSE, ".", 
		CosCensusSomatic)
	CosCensusSomatic <- gsub(TRUE, "YES", 
		CosCensusSomatic)

	CosCensusGermline <- lapply(info(raw.snp)$Gene.refGene, 
		`[`, 1) %in% metadata$cosmic_10Nov17_germline | 
		lapply(info(raw.snp)$Gene.refGene, 
			`[`, 2) %in% metadata$cosmic_10Nov17_germline
	CosCensusGermline <- gsub(FALSE, 
		".", CosCensusGermline)
	CosCensusGermline <- gsub(TRUE, "YES", 
		CosCensusGermline)

	if (length(grep("indel", raw.snp.vcf.files[i])) == 
		1) {
		ii <- "INDEL"
	} else {
		ii <- "SNP"
	} # are the variants in this VCF indels or snps?

	if (length(lapply(info(raw.snp)$Gene.refGene, 
		`[`, 2) == 1)) {
		ij <- lapply(info(raw.snp)$Gene.refGene, 
			`[`, 2)
	} else {
		ij <- "NA"
	} # 2nd gene

	# for this VCF file, combine extratced information and name columns
	if (nrow(raw.snp) > 0) {

		CV.out <- cbind(rep(raw.snp.vcf.files[i], 
			nrow(raw.snp)), rep(ii, nrow(raw.snp)), 
			lapply(info(raw.snp)$Gene.refGene, 
				`[`, 1), ij, SPV, positions, 
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
			"Indel or SNP", "Main gene", 
			"Downstream gene", "Somatic_P_value", 
			"Position", "Variant position in ClinVar", 
			"Variant position in Cosmic v82", 
			"In Cosmic Census Somatic Nov17", 
			"In Cosmic Census Germline Nov17", 
			"Variant gene in ACMG 57 gene list?", 
			"Maximum Population Frequency", 
			"ExAc All", "gnomAD exome All", 
			"gnomAD_genome All", "Mutation type", 
			"ANNOVAR gene position", 
			"num Germline REF reads", 
			"num Germline ALT reads", 
			"num Somatic REF reads", 
			"num Somatic ALT reads", 
			"Fraction of somatic reads that are ALT", 
			"Variant position in dbSNP147 ?", 
			"Variant position in repeat masked region?", 
			"Variant position in  ENCODE Dac unmappable region?", 
			"Variant position in regions of segmental duplication?", 
			"Variant supported by supplementary reads?", 
			"Variant effect - SIFT", 
			"Variant effect - Mutation_Taster", 
			"Variant effect - MetaLR", 
			"Variant effect - FATTHM")

		CV.tot <- rbind(CV.tot, CV.out)
	}
if (BAF.plot =="y") {	
	if (length(grep("snp", raw.snp.vcf.files[i])) == 
		1) { # a SNP VCF (not an Indel VCF)
		# then plot B allele frequency
		#==============
raw.snp <- o.raw.snp[(geno(o.raw.snp)$AD[, 
			2] + geno(o.raw.snp)$RD[, 
			2]) >= 20]
		raw.snp <- raw.snp[(geno(raw.snp)$AD[, 
			1] + geno(raw.snp)$RD[, 1]) >= 
			20]
		raw.snp <- raw.snp[unlist(info(raw.snp)$Func.refGene == 
			"exonic") | unlist(info(raw.snp)$Func.refGene == 
			"splicing")]
		Graw.snp <- raw.snp[(pbinom(pmin(as.numeric(geno(raw.snp)$RD[, 
			1]), as.numeric(geno(raw.snp)$AD[, 
			1])), as.numeric(geno(raw.snp)$AD[, 
			1] + geno(raw.snp)$RD[, 1]), 
			0.5)) >= 0.05 & (100 * geno(raw.snp)$AD[, 
			1]/(geno(raw.snp)$AD[, 1] + 
			geno(raw.snp)$RD[, 1])) <= 
			70 & (100 * geno(raw.snp)$AD[, 
			1]/(geno(raw.snp)$AD[, 1] + 
			geno(raw.snp)$RD[, 1])) >= 
			30]
		#==============
		a.snp <- strsplit(unlist(rowRanges(Graw.snp)@ranges@NAMES), 
			":")
		chr.snp <- unlist(lapply(a.snp, 
			`[`, 1))
		pos.snp <- unlist(lapply(strsplit(unlist(lapply(a.snp, 
			`[`, 2)), "_"), `[`, 1))
		gn <- unlist(lapply(info(Graw.snp)$Gene.refGene, 
			`[`, 1))
		Tsum <- sum(geno(Graw.snp)$AD[, 
			2] + geno(Graw.snp)$RD[, 
			2])
		Nsum <- sum(geno(Graw.snp)$AD[, 
			1] + geno(Graw.snp)$RD[, 
			1])
		norm.ratt <- Tsum/Nsum
		Treads <- (geno(Graw.snp)$AD[, 
			2] + geno(Graw.snp)$RD[, 
			2])
		NreadsADJ <- (geno(Graw.snp)$AD[, 
			1] + geno(Graw.snp)$RD[, 
			1]) * norm.ratt
		read.ratio <- Treads/NreadsADJ
		read.ratio <- read.ratio * 2 # x2 to make the ratio equivalent to copy number
		ALTratio <- geno(Graw.snp)$AD[, 
			2]/Treads
		forDC <- data.frame(chromosome = chr.snp, 
			start = pos.snp, stop = pos.snp, 
			x = ALTratio, y = read.ratio, 
			name = gn)
		forDC$start <- as.character(forDC$start)
		forDC$stop <- as.character(forDC$stop)
		forDC$chromosome <- as.character(forDC$chromosome)
		forDC$chromosome <- gsub("chr", 
			"", forDC$chromosome, fixed = TRUE)
		cf <- as.numeric(forDC$x)
		cf[cf > 0.5] <- 1 - cf[cf > 0.5]
		CNA.object <- CNA(as.numeric(cf), 
			as.numeric(forDC$chromosome), 
			as.numeric(forDC$start), 
			data.type = "logratio")
		seg.CNA.object <- segment(CNA.object)
		#==============
		f <- paste0(sub(".snp.vcf.gz","",sub(".indel.vcf.gz","",sub("varscan-","", raw.snp.vcf.files[i]))),".pdf")
		pdf(file = f, 12, 4)
		plotSampleCP(seg.CNA.object, 
			ylim = c(-0.05, 0.5), zeroline = FALSE, 
			xmaploc = TRUE, ylab = "Allele Frequency", 
			col = c("blue", "green"), 
			cex.axis = 1.5, cex.lab = 1.5, 
			clabpos = -0.05, pch = ".", 
			main = f, cex.main=0.75, cex.axis=0.75, cex.lab=0.75, cex=1)
		dev.off()
	}}
}
) # end system.time

#write file for excel use and to derive *.mut and CGI input files
# prepare
CV.totd<-data.frame(CV.tot)
tumour1<-sub(".snp.vcf.gz","",sub(".indel.vcf.gz","",sub("varscan-","", CV.totd$tumour))) # get tumour names
P.CGI<-sub(":", " ", as.character(CV.totd$Position)); P.CGI<-sub("_", " ", as.character(P.CGI)); P.CGI<-sub("/", " ", as.character(P.CGI))
S.CGI<-paste(CV.totd$tumour, CV.totd$Position, sep="_")
dfv<-t(data.frame(strsplit(P.CGI, " "))); chr.mut<-dfv[,1]; start.mut<-dfv[,2]; end.mut<-as.numeric(start.mut) - 1 + nchar(dfv[,4])
ccr<-paste0(CV.totd$num.Germline.REF.reads,"_", CV.totd$num.Germline.ALT.reads,"___", CV.totd$num.Somatic.REF.reads,"_", CV.totd$num.Somatic.ALT.reads)
pdr<-paste0(CV.totd$Variant.effect...SIFT,"_", CV.totd$Variant.effect...Mutation_Taster,"_", CV.totd$Variant.effect...MetaLR,"_", CV.totd$Variant.effect...FATTHM)
CV.tote<-data.frame(cbind(CV.totd, tumour1, P.CGI, S.CGI, chr.mut, start.mut, end.mut, ccr, pdr))
colnames(CV.tote)[32:39]<-c("sample","Position_for_CGI", "Sample_for_CGI", "chr", "start", "end", "GL-ref_GL-alt___SOM-ref_SOM-alt", "SIFT_MutationTaster_MetaLR_FATTHM")
CV.totf<-sapply(CV.tote, as.vector)

# then write
rownames(CV.totf) <- NULL
file <- paste(date(), "somatic_indel_OR_snv.txt", 
	sep = "_")
write.table(t(colnames(CV.totf)), file = file, 
	sep = "\t", row.names = "index", 
	col.names = FALSE, append = TRUE)
write.table(CV.totf, file = file, sep = "\t", 
	row.names = TRUE, col.names = FALSE, 
	append = TRUE)

# write file with no variants supported by supp reads
CV.totf.noSupp<-CV.totf[CV.totf[,27]==".",]
rownames(CV.totf.noSupp) <- NULL
file <- paste(date(), "NO.SUP.somatic_indel_OR_snv.txt", 
	sep = "_")
write.table(t(colnames(CV.totf.noSupp)), file = file, 
	sep = "\t", row.names = "index", 
	col.names = FALSE, append = TRUE)
write.table(CV.totf.noSupp, file = file, sep = "\t", 
	row.names = TRUE, col.names = FALSE, 
	append = TRUE)

# write file as CGI input file
CV.totf.noSupp.CGI<-cbind(CV.totf.noSupp[,34],str_split_fixed(as.character(CV.totf.noSupp[,33]), " ", n = Inf))
colnames(CV.totf.noSupp.CGI)<-c("sample", "chr", "pos", "ref", "alt")
rownames(CV.totf.noSupp.CGI) <- NULL
file <- paste(date(), "ForCGI_NO.SUP.somatic_indel_OR_snv.txt", 
	sep = "_")
write.table(t(colnames(CV.totf.noSupp.CGI)), file = file, 
	sep = "\t", row.names = F, 
	col.names = FALSE, append = TRUE)
write.table(CV.totf.noSupp.CGI, file = file, sep = "\t", 
	row.names = F, col.names = FALSE, 
	append = TRUE)