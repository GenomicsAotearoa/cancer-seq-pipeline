#Cris redefined (slightly) the plotSample function from the DNAcopy package
plotSampleCP = function (x, sampleid = NULL, chromlist = NULL, xmaploc = FALSE, 
                        col = c("black", "green"), pch = ".", cex = NULL, altcol = TRUE, 
                        segcol = "red", lwd = 3, zeroline = TRUE, zlcol = "grey", 
                        xlab = NULL, ylab = NULL, main = NULL, clabpos = 0.6, ...) {
    
  if (class(x) != "DNAcopy") 
      stop("First arg must be a DNAcopy object")
  if (missing(sampleid)) {
      sampleid <- 1
  }
  subx <- subset(x, chromlist = chromlist, samplelist = sampleid[1])
  genomdat <- subx$data[, 3]
  ina <- is.finite(genomdat)
  genomdat <- genomdat[ina]
  chrom <- subx$data[ina, 1]
  uchrom <- unique(chrom)
  segres <- subx$output
  if (xmaploc) {
      maploc <- subx$data[ina, 2]
      rmaploc <- sapply(uchrom, function(i, maploc, chrom) range(maploc[chrom == i]), maploc, chrom)
      nc <- length(uchrom)
      if ((nc > 1) && any(rmaploc[1, -1] < rmaploc[2, -nc])) {
          cmaploc <- cumsum(as.numeric(rmaploc[2, ]))
          for (i in 2:nc) {
              maploc[chrom == uchrom[i]] <- cmaploc[i - 1] + maploc[chrom == uchrom[i]]
          }
      }
      xlabel <- "Chromosomal Position"
    }
  else {
    maploc <- 1:sum(ina)
    xlabel <- "Index"
  }
  if (altcol & length(uchrom) > 1) {
    colvec <- rep(1, length(chrom))
    j <- 0
    for (i in uchrom) {
        j <- (j + 1)%%2
        colvec[chrom == i] <- j + 1
    }
  }
  else {
    colvec <- 1
  }
  if (missing(cex)) 
    cex <- ifelse(pch == ".", 3, 1)
  if (missing(main)) 
    main <- ""
  if (missing(xlab)) 
    xlab <- xlabel
  if (missing(ylab)) {
    if (attr(subx$data, "data.type") == "logratio") {
      ylab <- "log(relative CN)"
    }
    else {
      ylab <- "LOH"
    }
  }
  plot(maploc, genomdat, col = col[colvec], pch = pch, cex = cex,  main = main, xlab = xlab, ylab = ylab, ...)
        
        
        
  ii <- cumsum(c(0, segres$num.mark))
  mm <- segres$seg.mean
  kk <- length(ii)
  segments(maploc[ii[-kk] + 1], segres$seg.mean, x1 = maploc[ii[-1]], y1 = segres$seg.mean, col = segcol, lwd = lwd)
  if (zeroline) 
    abline(h = 0, col = zlcol, lwd = lwd)
        
        # add Cris' lines demarking chromosomes
  cbound = c(0,maploc[cumsum(table(chrom))])
  abline(v = cbound, col = "grey")
        
        # add Cris' chromosome names

  ch.na = vector(mode = "numeric",length = length(cbound)-1)
  for (i in 1:(length(cbound) - 1)) {
    ch.na[i] = cbound[i] + (cbound[i + 1] - cbound[i])/2
  }

  text(ch.na,rep(clabpos,22),labels = 1:22, cex = 0.75)       
}

library(DNAcopy)
data(coriell)
library("tidyverse")
counts = read_tsv("src/prosper03_3kb.count", col_names = FALSE)
colnames(counts) = c("chr", "tileStart","tileEnd","P1003A","P1003C","P1008A","P1008B","P1008E","P1008F","P1014A","P1014E","P1015A","P1015B","P1015C")
head(counts)

dnaObj = counts %>% 
  select(one_of("chr","tileStart","P1003A","P1003C"))

logratio = log2(dnaObj$P1003C) - log2(dnaObj$P1003A)

CNA.object = CNA( genomdat = logratio, chrom = dnaObj$chr, maploc = dnaObj$tileStart, data.type = 'logratio')
smoothed.CNA.object <- smooth.CNA(CNA.object, smooth.region = 100000, trim=0.25)
CNA.object <- segment(smoothed.CNA.object, verbose = 1)
plotSampleCP(x = CNA.object )
