summary.parBdev<-function(object, Samples, chr=c(1:22), ...)
 {  

  x<-object
  setwd(x)
  
  if (missing(Samples))
    Samples<-attr(x,"Samples")

  if (length(Samples)==1)
    Samples<-c(1,Samples)

  load("SBL/allSegments")
  
 
  ff<-function(x,chr)
   {
    cond<-x[,5]==chr & x$State!=0 
    return(x[cond,])
   } 

  ff2<-function(x,chr)
   {
    cond<-x[,5]==chr & x$State==0 
    return(x[cond,])
   } 


  ans<-list()
  no.mosaic<-list()
  for (i in 1:length(chr))
   {
    ans[[i]] <-lapply(res,FUN=ff,chr=chr[i])
    no.mosaic[[i]] <-lapply(res,FUN=ff2,chr=chr[i])
   }

  attr(ans,"no.mosaic")<-no.mosaic
  attr(ans,"Info")<-x
  attr(ans,"Samples")<-Samples
  attr(ans,"labels.samples")<-labels(x)
  attr(ans,"chr")<-chr
  class(ans)<-"summaryParBdev"
  names(ans)<-paste("chromosome",chr)
  ans
 }

