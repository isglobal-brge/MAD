
plotChr<-function(x, chr, sample, ...)
 {
   setwd(x)
   n<-attr(x,"Samples") 
   lab<-attr(x,"labels.samples")
   i<-c(1:n)[lab==sample] 
   load(paste("SBL/setupGADA",i,sep=""))
   load("SBL/gen.info.Rdata")
   attr(temp,"gen.info")<-gen.info
   plot.Bdev(temp, chr, ...)
   title(sample)
   title(paste("Chromosome",chr), line=0.3)   
 }


