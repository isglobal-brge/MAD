#' Zoom plot of LRR and BAF on regions of interest
#' 
#' @param x an object of class 'parGADA'
#' @param chr chromosome of interest
#' @param sample sample of interest
#' @param regions regions to be highlighted
#' @param delim same as xlim. Default is NULL
#' @return a pannel including two plots


plotZoomQMosaic <- function (x, chr, sample, regions, delim, ...)  {
  par(mfrow=c(2,1))
  ##NOT DETAILED
  plotQMosaic(x, chr, sample, regions, ...)
  ##DETAILED
  region.sel <- subset(regions, seqnames == paste0("chr", chr)
                       & sample == sample)
  start <- GenomicRanges::start(region.sel)
  start <- min(start)
  if (start - 2000000 < 0 ) start <- 0
  end <- GenomicRanges::end(region.sel)
  end <- max(end)
  
  plotQMosaic(x, chr, sample, regions, 
              delim=c(start-2000000, end+2000000), 
              xlim=c(start-2000000, end+2000000),
              title=FALSE, ...)
}
