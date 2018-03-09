setupGADA.B.deviation <-
function (file, NumCols, GenoCol, log2ratioCol, BAFcol, name.geno = c("AA", 
    "AB", "BB"), MarkerIdCol = 1, ChrNameCol = 2, ChrPosCol = 3, 
    sort = TRUE, orderProbes, sep = "\t", saveGenInfo = TRUE) 
{
    if (missing(GenoCol)) 
        stop("Missing GenoCol. Please, indicate which column of the file contains the genotypes")
    if (missing(BAFcol)) 
        stop("Missing BAFcol. Please, indicate which column of the file contains the BAF")
    if (missing(NumCols)) 
        stop("Missing NumCols argument. Please, indicate the number of columns in the file")
    ans <- list()
    xx <- scan(file, skip = 1, what = c("character"), sep = sep)
    headers <- scan(file, nline = 1, what = c("character"), quiet = TRUE, 
        sep = sep)
    if (headers[MarkerIdCol] != "Name") 
        warning("Expecting 'Name' as the header of MarkerIdCol in an Illumina file")
    if (headers[ChrNameCol] != "Chr") 
        warning("Expecting 'Chr' as the header of ChrNameCol in an Illumina file")
    if (!headers[ChrPosCol] %in% c("position", "Position")) 
        warning("Expecting 'position' as the header of ChrPosCol in an Illumina file")
    if (NROW(grep("GType", headers[GenoCol])) != 1) 
        warning("Expecting 'GType' as the header of  GenoCol in an Illumina file")
    if (NROW(grep("Log.R.Ratio", headers[log2ratioCol])) != 1) 
        warning("Expecting 'Log R Ratio' as the header of  log2ratioCol in an Illumina file")
    if (NROW(grep("B.Allele.Freq", headers[BAFcol])) != 1) 
        warning("Expecting 'B.Allele.Freq' as the header of  BAFcol in an Illumina file")
    x <- matrix(xx, ncol = NumCols, nrow = length(xx)/NumCols, 
        byrow = TRUE)
    gg <- x[, GenoCol]
    baf <- as.numeric(x[, BAFcol])
    baf[is.na(baf)] <- -999
    nProbes <- length(gg)
    outC <- .C("bDeviation", as.character(gg), as.double(baf), 
        as.integer(nProbes), as.double(rep(0, nProbes)), PACKAGE = "mad")
    b.deviation <- outC[[4]]
    b.deviation[b.deviation %in% c(-9, -999, 999, 999.5, -999.5)] <- NA
    if (!sort) {
        x[, ChrNameCol][x[, ChrNameCol] == "XY"] <- "X"
        chr <- factor(x[, ChrNameCol], levels = c(as.character(1:22), 
            "X", "Y"))
        temp <- data.frame(probe = x[, MarkerIdCol], chr = chr, 
            pos = as.numeric(x[, ChrPosCol]), geno=gg, stringsAsFactors = FALSE)
        mito <- is.na(temp$chr)
        temp2 <- temp[(!mito), ]
        ans$log.ratio <- as.numeric(x[!mito, log2ratioCol])
        ans$B.allele.freq <- as.numeric(x[!mito, BAFcol])
        Bdev <- b.deviation[(!mito)]
        geno <- x[!mito, GenoCol] 
        attr(ans, "gen.info") <- temp2
    }
    else {
        x[, ChrNameCol][x[, ChrNameCol] == "XY"] <- "X"
        chr <- factor(x[, ChrNameCol], levels = c(as.character(1:22), 
            "X", "Y"))
        pos <- as.numeric(x[, ChrPosCol])
        if (missing(orderProbes)) 
            o <- order(chr, pos)
        else o <- orderProbes
        temp <- data.frame(probe = x[o, MarkerIdCol], chr = chr[o], 
            pos = pos[o], geno=gg[o], stringsAsFactors = FALSE)
        mito <- is.na(temp$chr)
        temp2 <- temp[(!mito), ]
        aux <- as.numeric(x[o, log2ratioCol])
        ans$log.ratio <- aux[!mito]
        aux2 <- as.numeric(x[o, BAFcol])
        ans$B.allele.freq <- aux2[!mito]
        b.deviation.sorted <- b.deviation[o]
        Bdev <- b.deviation.sorted[(!mito)]
        aux3 <- x[o, GenoCol]
        geno <- aux3[!mito]
        attr(ans, "gen.info") <- temp2
    }
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

    mask <- chr%in%c("X", "Y") & geno!="AB"
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

