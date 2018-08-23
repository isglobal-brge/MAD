#' A parallel version of setupGADA.B.deviation function
#' @param folder The folder where data is stored. Not required if the working directory contains a 'rawData' folder
#' @param files The names of the files with pennCNV-format fiels. Not required. By default all files in the 'rawData' folder are analyzed
#' @param verbose Should information about process be printed in the console? The default is TRUE
#' @param sort Should data be sorted by genomic position? Default is TRUE
#' @param MarkerIdCol The column in 'file' containing the name of the marker. Default first column
#' @param ChrNameCol The column in 'file' containing the chromosome. Default second column
#' @param ChrPosCol The column in 'file' containing the genomic position. Default third column
#' @param mc.cores number of cores to be used when using multiple cores (see argument 'mc.cores' from 'mclapply' function of parallel library) 
#' @param ... Other arguments passed through 'setupGADAIllumina'

setupParGADA.B.deviation <- 
  function(folder, files, verbose=TRUE,
           sort=TRUE, MarkerIdCol = 1, 
           ChrNameCol = 2, ChrPosCol = 3,
           mc.cores=1, ...)
 {
 
  if (!"SBL"%in%dir() )
   dir.create("SBL")
  
  if(missing(folder))
    folder<-getwd()
  oldpwd<-getwd()
  setwd(folder)

  if (!"rawData"%in%dir())
   stop("rawData folder with the data files cannot be located")  

  if (missing(files))
    files <- dir("rawData")

  if (verbose)
   {
    cat("\n")
    cat("Creating object with annotation data ... \n")
   }

  gen.info <- fread(file.path("rawData", files[1]), 
                       select=c(MarkerIdCol, 
                                ChrNameCol, 
                                ChrPosCol),
                       data.table = FALSE)
  colnames(gen.info) <- c("marker", "chr", "position")
  

  if (sort)
   {
    
# mitocondrial?
    mito <- is.na(gen.info$chr) 
    gen.info <- gen.info[!mito,]

    gen.info$chr[gen.info$chr=="XY"] <- "X"
    
    gen.info <- gen.info[gen.info$chr%in%
                           c(as.character(1:22), "X", "Y"),]
    
    o<-order(gen.info$chr, gen.info$position)
    gen.info<-gen.info[o,]
    attr(gen.info,"sort") <- TRUE
    attr(gen.info,"orderProbe") <- o
    select <- rownames(gen.info)
   }

  else
   {
    mito <- is.na(gen.info$chr) 
    gen.info <- gen.info[!mito,]
    gen.info$chr[gen.info$chr=="XY"] <- "X"
    
    gen.info <- gen.info[gen.info$chr%in%
                           c(as.character(1:22), "X", "Y"),]
    
    attr(gen.info,"sort") <- FALSE
    attr(gen.info,"orderProbe") <- 1:nrow(gen.info)
    select <- rownames(gen.info)
   } 
  
  save(gen.info, file="SBL/gen.info.Rdata", compress=TRUE)
   

  if (verbose)
   {
    cat("Creating object with annotation data ...done \n")
   }


 if (verbose)
  {
    cat("\n")
    cat("Creating objects of class setupGADA for all input files... \n")
  }


 if (verbose)
  {
  cat("  Applying setupGADA.B.deviation for", length(files) ,"samples ... \n")
  }



  prepare.i<-function(i, files, MarkerIdCol = MarkerIdCol,
                      chrNameCol = chrNameCol, 
                      ChrPosCol = ChrPosCol,
                      ...)
    {
      if (verbose)
       cat("  Importing array: ",files[i],"... ")  

      dd<-paste("rawData/", files[i], sep="")
      temp<-setupGADA.B.deviation(dd, saveGenInfo=FALSE, ...)

      save(temp, file=paste("SBL/setupGADA",i,sep="" ), compress=TRUE)

      if (verbose)
       cat("   Array #",i,"...done \n")  
    }


 res <- mclapply(1:length(files), 
               function(i) try(prepare.i(i, files=files, 
                                         orderProbes=select, ...), TRUE),
               mc.cores=mc.cores)

 if (verbose)
  {
   cat("  Applying setupGADA.B.deviation for", length(files) ,"samples ... done \n")
   cat("Creating objects of class setupGADA for all input files... done \n")
  }


   error<-sum(unlist(lapply(res, function(x) inherits(x, "try-error"))))

   if (error>0)
    {
      cat("WARNING!!! \n")
      cat("  Creating objects procedure failed for",sum(error),"samples \n")
      cat("  (type error to see what happened) \n")
      error <<- res
    }


 ans<-getwd()
 class(ans)<-"parGADA"
 attr(ans,"type")<-"Illumina"
 attr(ans,"labels.samples")<-gsub("sample.","",gsub(".txt","",files))
 attr(ans,"Samples")<-length(files)
 attr(ans,"b.deviation")<-TRUE
 ans

}


