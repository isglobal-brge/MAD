parBE.B.deviation<-function(x, Samples, T, MinSegLen,
                            verbose=TRUE, mc.cores=1, ...) {

  if (missing(Samples))
   Samples<-attr(x,"Samples") 

  if (verbose)
    cat("Retrieving annotation data ...")

  load(file.path(x, "SBL/gen.info.Rdata"))

  if (verbose)
    cat("done \n")

  if (length(Samples)>2)
   stop(" 'Samples' must be the number of samples or a vector indicating the first and last sample")

  if (length(Samples)==1)
    Samples<-c(1,Samples)
  
  analize.i<-function(i,T, MinSegLen, gen.info, labels, verbose)
    {
      if (verbose)
       cat("   Array #",i,"... ")  
      load(file.path(x, paste0("SBL/sbl",i)))
      attr(step1,"gen.info")<-gen.info
      step2<-BackwardElimination(step1, T=T, MinSegLen=MinSegLen)
      load(file.path(x, paste0("SBL/setupGADA",i)))
      ans<-summary.Bdeviation(step2, temp, print=FALSE, ...)
      save(ans, step2, file=file.path(x, paste0("SBL/segments",i)),
           compress=TRUE)
      ans$sample<-labels[i]
      ans     
    }

   if (verbose)
     cat("Backward elimination procedure for",
         Samples[2]-Samples[1]+1, "samples ... \n")

   res <- mclapply(Samples[1]:Samples[2], 
                                 analize.i, T=T, 
                                 MinSegLen=MinSegLen, 
                                 gen.info=gen.info,
                                 labels=attr(x,"labels.samples"),
                                 verbose=verbose,
                                 mc.cores = mc.cores)

   if (verbose)
     cat("Backward elimination procedure for",
         Samples[2]-Samples[1]+1, "samples ...done \n")

   error<-sum(unlist(lapply(res, function(x) inherits(x, "try-error"))))

   if (error>0)
    {
      cat("WARNING!!! \n")
      cat("  Backward Elimination procedure failed for",sum(error),"samples \n")
      cat("  (type error to see what happened) \n")
      error <<- res
    }

  save(res, file=file.path(x, paste0("SBL/allSegments")),compress=TRUE)
 }

