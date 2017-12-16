
plot.Bdev<-function(x, chr, ...)
 {

  if (missing(chr))
   stop("Please, select a chromosome")

  gen.info<-attr(x,"gen.info")

  o<-gen.info$chr==chr
 
  pos<-gen.info$pos[o]


  par(mar=c(5, 4, 4, 4) + 0.1)
  plot(pos, x$LRR[o], ylim=c(-1,1), las=1, pch=1, cex=0.25, col="black", ylab="", xlab="", main="", ...)

    

  par(new=TRUE)
  plot(pos, x$B.allele.freq[o], col="red", pch=1, cex=0.25, ylab="", xlab="", main="", axes=F, ...) 

  mtext("LRR",side=2,col="black",line=2.5,adj=0.5)

  mtext("BAF",side=4,col="black",line=2.5,adj=0.5)
  axis(4,at=c(0,0.25,0.5,0.75,1),labels=c("0.0",0.25,0.5,0.75,"1.0"), las=1,col="black",col.axis="black")


  mtext("Chromosome coordinates",side=1,col="black",line=2.5)


 } 


