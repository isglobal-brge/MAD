print.summaryParBdev<-function(x, ...)
 {
  if (!inherits(x, "summaryParBdev"))
   stop("object must be of class 'summaryParBdev' ")

  no.mosaic<-attr(x,"no.mosaic")

  nInd<-attr(x,"Samples")[2] 

  out<-list()
  mosaic<-list()

  for (i in 1:22)
  # No cromosoma X Y -> Luego preparar por cromosoma
  {
   tt<-x[[i]]
   UPD <- 0
   LOH <- 0
   n.mosaic <- 0
   for (j in 1:nInd)
    {
     tt2<-tt[[j]]
     n.mosaic.j<-nrow(tt2)
     n.mosaic<-n.mosaic+n.mosaic.j

     UPD<-UPD+sum(tt2$State==1)
     LOH<-LOH+sum(tt2$State==2)
    }

   mosaic[[i]]<-n.mosaic

   out[[i]]<-c(n.mosaic,UPD, LOH)

   n.mosaicT <- sum(unlist(lapply(out, function(x) x[1])))
   UPDT <- sum(unlist(lapply(out, function(x) x[2])))
   LOHT <- sum(unlist(lapply(out, function(x) x[3])))
   per.UPDT<-round((UPDT/(n.mosaicT))*100,1)
   per.LOHT<-round((LOHT/(n.mosaicT))*100,1)

    res<-list(c(n.mosaicT, UPDT, per.UPDT, LOHT, per.LOHT),mosaic,out)
  }
  
  cat("\n")
  cat("------------------------------------------- \n")
  cat("Summary results for", nInd, "individuals \n")
  cat("------------------------------------------- \n")


  total<-data.frame(t(res[[1]]))
  names(total)<-c("# segments","Mosaicism","%","LOH","%")
  dimnames(total)[[1]]<-""

  cat("\n")
  cat("Number of Total Segments:","\n")
  print(total)


  ii<-t(data.frame(res[[3]]))
  dimnames(ii)[[1]]<-paste("Chromosome",1:22)
  dimnames(ii)[[2]]<-c("segments", "Mosaicism", "LOH")

  cat("\n")
  cat("Number of Total Segments by chromosome:","\n")
  print(ii)

  cat("\n")


  invisible(res)
 }

