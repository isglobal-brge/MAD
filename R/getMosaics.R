getMosaics <- function(x, ...){
  
  load(file.path(x, "SBL/allSegments"))
  
  out <- rbindlist(res)
  out <- out[out$State!=0,]
  
  out.gr <- GenomicRanges::GRanges(seqnames = paste0("chr", out$chr),
                                   ranges = IRanges:::IRanges(
                                     start = out$IniProbe,
                                     end = out$EndProbe),
                                   LenProbe = out$LenProbe,
                                   LRR = out$LRR,
                                   LRR.se = out$`(s.e.)`,
                                   Bdev = out$Bdev,
                                   State = out$State, 
                                   sample = out$sample)
  out.gr
}
