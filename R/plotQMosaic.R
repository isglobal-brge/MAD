#' Plot LRR and BAF of a given chromosome
#' 
#' @param x an object of class 'parGADA'
#' @param chr chromosome of interest
#' @param sample sample of interest
#' @param regions regions to be highlighted
#' @param delim same as xlim. Default is NULL
#' @param title should title be added? Default is TRUE
#' @return a plot



plotQMosaic <- function (x, chr, sample, regions, 
                         delim=NULL, title=TRUE,
                         ...)  {
  if (missing(chr)) 
    stop("Please, select a chromosome")
  setwd(x)
  n <- attr( x, "Samples")
  lab <- attr( x, "labels.samples")
  i <- c( 1:n)[ lab == sample]
  load( paste( "SBL/setupGADA", i, sep = ""))
  gen.info <- attr( x, "gen.info")
  o <- gen.info$chr == chr
  pos <- gen.info$pos[o]
  if (!missing(regions))
    region.sel <- regions[regions$chr == chr & regions$sample == sample, ]
  par(mar = c(5, 4, 4, 4) + 0.1)
  plot(pos, temp$LRR[o], ylim = c(-2, 2), las = 1, pch =".", cex = 2, col = 1, ylab = "",xlab = "", main = "", xaxt="n", yaxt="n", ...)
  axis( 2, at = c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2), labels = c("-2.0", -1.5, "-1.0", -0.5, "0.0", 0.5, "1.0", 1.5, "2.0"), las = 1, col = "black", col.axis = "black")
  colBAF <- rep( 2, length(temp$geno[o]))
  colBAF[ temp$geno[o] == "AA" ] <- 2
  colBAF[ temp$geno[o] == "BB" ] <- 2
  if (!missing(regions)) {
    start <- region.sel$IniProbe
    end <- region.sel$EndProbe
    abline(v=c(start, end), lwd=1)
    for (i in 1:nrow(region.sel)){
      colBAF[ temp$geno[o] == "AB" & pos > start[i] & pos < end[i] ] <- 2
    }
  }
  par(new = TRUE)
  plot(pos, temp$B.allele.freq[o], col = colBAF, pch = ".", cex = 2, ylab = "", xlab = "", main = "", axes = F, ...)
  #if (!missing(regions)) {
  #    for (i in 1:nrow(region.sel)){
  #        u <- gen.info$chr == chr & gen.info$pos > region.sel$IniProbe[i] & gen.info$pos < region.sel$EndProbe[i]
  #        upper <- mean(temp$B.allele.freq[ u ][ temp$B.allele.freq[u] > 0.5& temp$geno[u] == "AB"])
  #        bottom <- mean(temp$B.allele.freq[ u ][ temp$B.allele.freq[u] < 0.5& temp$geno[u] == "AB"])
  #        lines(x=c(start[i], end[i]), y=rep(upper,2), col=2, lty=2)
  #        lines(x=c(start[i], end[i]), y=rep(bottom,2), col=2, lty=2)
  #    }
  #}
  abline(h=0.5, col=8)
  abline(h=c(0.33, 0.66), col=8)
  
  mtext("LRR", side = 2, col = "black", line = 2.5, adj = 0.5)
  mtext("BAF", side = 4, col = "black", line = 2.5, adj = 0.5)
  axis( 4, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c("0.0", 0.2, 0.4, 0.6, 0.8, "1.0"), las = 1, col = "black", col.axis = "black")
  xaxis <- seq(min(pos), max(pos), length.out=10)
  if (!is.null(delim)) xaxis <- seq(delim[1], delim[2], length.out=5)
  axis( 1, at = xaxis, labels = as.character(round(xaxis/1000000, 1)), las = 1, col = "black", col.axis = "black")
  mtext("position (Mb)", side = 1, col = "black", line = 2.5)
  if (title){
    title(sample)
    title(paste("Chromosome", chr), line = 0.3)
  }  
}

