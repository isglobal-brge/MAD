#' Plot LRR and BAF of a given chromosome
#' 
#' @param x an object of class 'parGADA'
#' @param chr chromosome of interest
#' @param sample sample of interest
#' @param regions regions to be highlighted
#' @param delim same as xlim. Default is NULL
#' @param title should title be added? Default is TRUE
#' @return a plot highlighting the altered regions



plotQMosaic <- function (x, chr, sample, regions, delim=NULL, title=TRUE,
                         col.dots = c("black", "red"), ...)  {
  if (missing(chr)) 
    stop("Please, select a chromosome")
  n <- attr( x, "Samples" )
  lab <- attr( x, "labels.samples" )
  i <- c(1:n)[lab == sample]
  load(file.path(x, paste0("SBL/setupGADA", i)))
  load(file.path(x, paste0("SBL/gen.info.Rdata")))
  o <- gen.info$chr == chr
  pos <- gen.info$pos[o]
  if (!missing(regions))
    region.sel <- subset(regions, seqnames == paste0("chr", chr)
                          & sample == sample)
  par(mar = c(5, 4, 4, 4) + 0.1)
  plot(pos, temp$LRR[o], ylim = c(-2, 2), las = 1, pch =".", cex = 2, 
       col = col.dots[1], ylab = "",xlab = "", main = "", xaxt="n", yaxt="n", ...)
  axis( 2, at = c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2), 
        labels = c("-2.0", -1.5, "-1.0", -0.5, "0.0", 0.5, "1.0", 1.5, "2.0"), 
        las = 1, col = col.dots[1])
  
  par(new = TRUE)
  plot(pos, temp$B.allele.freq[o], col = col.dots[2], pch = ".", 
       cex = 2, ylab = "", xlab = "", main = "", axes = F, ...)
  abline(h=0.5, col=8)
  abline(h=c(0.33, 0.66), col="gray70", lwd=1)
  
  mtext("LRR", side = 2, col = col.dots[1], line = 2.5, adj = 0.5)
  mtext("BAF", side = 4, col = col.dots[2], line = 2.5, adj = 0.5)
  axis( 4, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels = c("0.0", 0.2, 0.4, 0.6, 0.8, "1.0"), 
        las = 1, col.axis = col.dots[2])
  xaxis <- seq(min(pos), max(pos), length.out=10)
  if (!is.null(delim)) xaxis <- seq(delim[1], delim[2], length.out=5)
  axis( 1, at = xaxis, labels = as.character(round(xaxis/1000000, 1)), las = 1)
  mtext("position (Mb)", side = 1, line = 2)
  
  if (!missing(regions)) {
    start <- GenomicRanges::start(region.sel)
    end <- GenomicRanges::end(region.sel)
    abline(v=c(start, end), lwd=2, col="gray70")
  }
  
  if (title){
    title(sample)
    title(paste("Chromosome", chr), line = 0.3)
  }  
}

