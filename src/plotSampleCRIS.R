#Cris redefined (slightly) the plotSample function from the DNAcopy package
library(DNAcopy)

plotSampleCP = function (x, outputFile, col = c("black", "green"), pch = ".", cex = NULL, altcol = TRUE,  segcol = "red", lwd = 3, zeroline = TRUE, zlcol = "grey",  clabpos = 0.6 ) {
    
  if (class(x) != "DNAcopy") 
      stop("First arg must be a DNAcopy object")
  sampleid = 1
  if (missing(outputFile))
      stop("please provide an output filename when callinging this function")
 
  print(outputFile) 	

  subx <- subset(x, samplelist = sampleid)
  genomdat <- subx$data[, 3]
  ina <- is.finite(genomdat)
  genomdat <- genomdat[ina]
  chrom <- subx$data[ina, 1]
  uchrom <- unique(chrom)
  segres <- subx$output
  maploc <- 1:sum(ina)
  xlabel <- "Index"

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

  cex <- ifelse(pch == ".", 3, 1)
  main <- ""
  xlab = xlabel
  ylab <- "log(relative CN)"

  pdf(file = outputFile, width=14, height=7)

  plot(maploc, genomdat, col = col[colvec], pch = pch, cex = cex,  main = main, xlab = xlab, ylab = ylab)
        
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

  dev.off()
}
