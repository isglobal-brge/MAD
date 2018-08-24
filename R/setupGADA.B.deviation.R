#' Imports array data to MAD
#' @param file A pennCNV-style file 
#' @param NumCols The number of columns in 'file'
#' @param GenoCol The column in 'file' containing genotypes
#' @param log2ratioCol The column in 'file' containing log2ratio intensities
#' @param BAFCol The column in 'file' containing B-allele frequency (BAF)
#' @param name.geno the name of the genotypes (default "AA", "AB" and "BB)
#' @param MarkerIdCol The column in 'file' containing the name of the marker. Default first column
#' @param ChrNameCol The column in 'file' containing the chromosome. Default second column
#' @param ChrPosCol The column in 'file' containing the genomic position. Default third column
#' @param sort should data be sorted by genomic position? The default is TRUE.
#' @param orderProbes used to avoid sorting
#' @param sep how data columns are separated. Default by tab
#' @param saveGenInfo Should annotation data be saved? The default is TRUE. See details

setupGADA.B.deviation <-
function (file, NumCols, GenoCol, log2ratioCol, BAFcol, 
          name.geno = c("AA", "AB", "BB"), MarkerIdCol = 1, 
          ChrNameCol = 2, ChrPosCol = 3, sort = TRUE, 
          orderProbes, sep = "\t", saveGenInfo = TRUE) 
{
    if (missing(GenoCol)) 
        stop("Missing GenoCol. Please, indicate which column of the file contains the genotypes")
    if (missing(BAFcol)) 
        stop("Missing BAFcol. Please, indicate which column of the file contains the BAF")
    if (missing(NumCols)) 
        stop("Missing NumCols argument. Please, indicate the number of columns in the file")
    ans <- list()
    
    x <- fread(file, data.table = FALSE)
    if (!missing(orderProbes))
     x <- x[orderProbes,]
    else {
      if (sort)
        x <- x[order(x[,ChrNameCol], x[, ChrPosCol]), ]
    }
    
    gg <- x[, GenoCol]
    baf <- as.numeric(x[, BAFcol])
    baf[is.na(baf)] <- -999
    nProbes <- length(gg)
    outC <- .C("bDeviation", as.character(gg), as.double(baf), 
        as.integer(nProbes), as.double(rep(0, nProbes)), PACKAGE = "mad")
    b.deviation <- outC[[4]]
    b.deviation[b.deviation %in% c(-9, -999, 999, 999.5, -999.5)] <- NA
    
    ans$log.ratio <- as.numeric(x[, log2ratioCol])
    ans$B.allele.freq <- as.numeric(x[, BAFcol])
    Bdev <- b.deviation
    geno <- x[, GenoCol] 
    
    temp <- Bdev
  
    na.find <- is.na(temp)
    if (any(na.find))
      {
        temp[na.find & geno==name.geno[1]] <- 0
        temp[na.find & geno==name.geno[2]] <- 0.5
        temp[na.find & geno==name.geno[3]] <- 1
      }


    Bdev[temp == 0] <- 0
    qng <- qnorm(Bdev)
    
    m <- median(qng, na.rm = TRUE)
    if (is.infinite(m) & m > 0)
      m <- 100
    if (is.infinite(m) & m < 0)
      m <- -100  
        
    Bdev <- qng - m
    Bdev[temp == 0] <- 0
    Bdev[is.na(Bdev)] <- 0

    mask <- x[,ChrNameCol]%in%c("X", "Y") & geno!="AB"
    Bdev[mask] <- NA
    

    ans$LRR <- ans$log.ratio
    ans$log.ratio <- Bdev
    ans$Bdev.original.scale <- temp
    ans$Bdev.original.scale[geno!=name.geno[2]] <- NA
    
    attr(ans, "type") <- "Illumina"
    attr(ans, "Bdev") <- TRUE
    if (!saveGenInfo) 
      {
        attr(ans, "gen.info") <- TRUE
        ans$geno <- geno
      }
    class(ans) <- "setupGADA"
    ans
}

