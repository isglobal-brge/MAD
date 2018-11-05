plotSegments <- function(x, ...){
  xx <- subset(x, State%in%c(1,2,3,4))
  if (length(x)!=length(xx)) 
   warning("Only mosaic gains and loses are depicted")
  xx <- as(xx, "data.frame")[,c(1,2,3,10,11)]
  names(xx) <- c("chromosome", "start", "end",
                 "segmean", "sample") 
  xx$segmean <- xx$segmean
  p <- GenVisR::cnSpec(xx, CNscale = "relative")
  p + scale_fill_gradientn("Mosaicims" , 
                           colours = c("orange", "darkgreen",
                                       "darkblue",
                                       "tomato"), 
                           values = scales::rescale(c(1,4)), 
                           limits = c(1,4), oob = scales::squish) +
    guides(fill=guide_legend(title="Mosaicism"))
}

