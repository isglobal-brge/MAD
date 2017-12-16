exportSegments2File<-function(x, file, allSegments=FALSE, ...)
 {  
  setwd(x)
  load("SBL/allSegments")

  if (missing(file))
   file<-paste("segments_", deparse(substitute(x)), ".txt", sep="")

  out<-NULL
  for (i in 1:length(res))
   {
    temp.i<-res[[i]]
    if (!allSegments)
     temp<-temp.i[temp.i$State!=0,]
    out<-rbind(out, temp)
   }

   write.table(out, file=file, row.names=FALSE, quote=FALSE, sep="\t")
 }
