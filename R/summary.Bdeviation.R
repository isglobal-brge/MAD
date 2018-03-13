summary.Bdeviation <-function(object, object2, which="both", T, BaseAmp, LRR.threshold=c(-0.09, 0.09), Het.threshold=5,  Hom.threshold=85, print=TRUE, ...)
 {
  
  if(!inherits(object, "BackwardElimination"))
   stop("an object of class 'BackwardElimination' is required")

  if (missing(object2))
   stop(" The argument 'object2' should be given. See the manual")

  which.ok<-match(which,c("Gains","Loses","both"),nomatch=0)
  if (which.ok==0)
   stop ("which should be 'Gains', 'Loses' or 'both' ")
  
  gen.info<-attr(object,"gen.info")

  LRR<-object2$LRR
  Bdev.orig<-object2$Bdev.original.scale
  BAF<-object2$B.allele.freq
  geno.het<-object2$geno.het

  ff<-function(x, y) 
   {
     ini<-x[1]
     end<-x[2] 
     ans<-mean(y[ini:end], na.rm=TRUE)
     ans
   }

  ff2<-function(x, y) 
   {
     ini<-x[1]
     end<-x[2] 
     ans<-sd(y[ini:end], na.rm=TRUE)
     ans
   }

  ff3<-function(x, y, limit=c(0.48,0.52)) 
   {
# % BAF in c(0.48,0.52) %HetWT
     ini<-x[1]
     end<-x[2] 
     temp<-y[ini:end]
     ans<- mean(temp>limit[1] & temp<limit[2], na.rm=TRUE)
     ans
   }


  ff4<-function(x, y, limit=c(0.10,0.90)) 
   {
# % BAF in c(0.10,0.90) %Homocigosity
     ini<-x[1]
     end<-x[2] 
     temp<-y[ini:end]
     ans<- mean(temp<limit[1] | temp>limit[2], na.rm=TRUE)
     ans
   }



#RPR Pulling out the segments, and segments amplitudes
  if (!is.null(gen.info))
   {

#    k<-unique(gen.info$chr)
#    k<-k[!is.na(k)]

# Only autosomes
#    k<-c(1:22)

# All Chromosomes
    k <- c(1:24) 

    Segments<- WextIextToSegments(object[[1]])
    Segments$chr <- attr(object,"chr")[[1]]
    o<-gen.info$chr==1
    LRR.i<-LRR[o]
    Bdev.orig.i<-Bdev.orig[o]
    BAF.i<-BAF[o]

    ax<-round(apply(Segments, 1, ff, y=LRR.i),2)
    ax2<-round(apply(Segments, 1, ff2, y=LRR.i),2)
    ax3<-round(apply(Segments, 1, ff, y=Bdev.orig.i),3)
    ax4<-round(apply(Segments, 1, ff3, y=BAF.i)*100,1)
    ax5<-round(apply(Segments, 1, ff4, y=BAF.i)*100,1)

    Segments$LRR<-ax
    Segments$"(s.e.)"<-ax2
    Segments$"Bdev"<-ax3
    Segments$"%HetWT"<-ax4
    Segments$"%Hom"<-ax5


    for (i in 2:length(k))
     {
       
      temp<-WextIextToSegments(object[[i]]) 
      chr <- attr(object,"chr")[[i]] 
      o<-gen.info$chr==chr
      LRR.i<-LRR[o]
      Bdev.orig.i<-Bdev.orig[o]
      BAF.i<-BAF[o]

      ax<-round(apply(temp, 1, ff, y=LRR.i),2)      
      ax2<-round(apply(temp, 1, ff2, y=LRR.i),2)
      ax3<-round(apply(temp, 1, ff, y=Bdev.orig.i),3)
      ax4<-round(apply(temp, 1, ff3, y=BAF.i)*100,1)
      ax5<-round(apply(temp, 1, ff4, y=BAF.i)*100,1) 

      temp$chr <- chr
      temp[,1]<-temp[,1]+Segments[nrow(Segments),2]  
      temp[,2]<-temp[,2]+Segments[nrow(Segments),2]
      temp$LRR<-ax
      temp$"(s.e.)"<-ax2
      temp$Bdev<-ax3
      temp$"%HetWT"<-ax4
      temp$"%Hom"<-ax5

      Segments<-rbind(Segments,temp)
     }

    }
  else 
   {
    Segments<-WextIextToSegments(object)
    Segments$chr <- 0 #attr(object,"chr")[[i]] #FIXME If no genome.info.... I don't know the chr... maybe we can set 0 as a mark... 
   }

   K <- nrow(Segments)

   sigma2<-object$sigma2



# All chromosomes
       Segments.autosomes<-Segments[Segments$chr%in%c(1:22),]
       K.autosomes<-nrow(Segments.autosomes)
       Segments.X<-Segments[Segments$chr%in%c("X"),]
       K.X<-nrow(Segments.X)
       Segments.Y<-Segments[Segments$chr%in%c("Y"),]
       K.Y<-nrow(Segments.Y)


# Estimation of the reference level as the median across all probes
   if (missing(BaseAmp)) 
   {

# For autosomes 
     aux<-.C("RcallCompAmpMedianMethod",
            SegLen=as.integer(Segments.autosomes$LenProbe),
            SegAmp=as.double(Segments.autosomes$MeanAmp),
            K=as.integer(K.autosomes),
            BaseAmp=double(1),PACKAGE="gada")

     BaseAmp <- ifelse(is.na(aux$BaseAmp),0,aux$BaseAmp)


# For X and Y chromosomes
    if (!is.null(gen.info))
     {
       aux.X<-.C("RcallCompAmpMedianMethod",
              SegLen=as.integer(Segments.X$LenProbe),
              SegAmp=as.double(Segments.X$MeanAmp),
              K=as.integer(K.X),
              BaseAmp=double(1),PACKAGE="mad")

       BaseAmp.X <- ifelse(is.na(aux.X$BaseAmp),0,aux.X$BaseAmp)


       aux.Y<-.C("RcallCompAmpMedianMethod",
              SegLen=as.integer(Segments.Y$LenProbe),
              SegAmp=as.double(Segments.Y$MeanAmp),
              K=as.integer(K.Y),
              BaseAmp=double(1),PACKAGE="mad")

       BaseAmp.Y <- ifelse(is.na(aux.Y$BaseAmp),0,aux.Y$BaseAmp)

       BaseAmp.all<-c(BaseAmp,BaseAmp.X,BaseAmp.Y)
      }
     
     else
      {
       BaseAmp.all<-c(BaseAmp, NA, NA)
      }

   }
  
  else
   {
     BaseAmp.all<-BaseAmp
   }



  if (!is.null(gen.info))
  {
# if T missing use x$T
   if (missing(T))
     T <- object[[1]]$T
   sigma2 <- object[[1]]$sigma2
   MinSegLen <- object[[1]]$MinSegLen
  }

  else
  {
# if T missing use x$T
   if (missing(T))
     T <- object$T
   sigma2 <- object$sigma2
   MinSegLen <- object$MinSegLen
  }


   aux<-.C("RcallClassifySegments",
            SegLen=as.integer(Segments.autosomes$LenProbe),            
            SegAmp=as.double(Segments.autosomes$MeanAmp),
            SegState=double(K.autosomes),
            K=as.integer(K.autosomes),
            BaseAmp=as.double(BaseAmp.all[1]),
            sigma2=as.double(sigma2),
            T=as.double(T),PACKAGE="mad")                   

   Segments.autosomes$State <- aux$SegState



  if (!is.null(gen.info)) # JRG Oct'11
    {

      aux.X<-.C("RcallClassifySegments",  
               SegLen=as.integer(Segments.X$LenProbe),            
               SegAmp=as.double(Segments.X$MeanAmp),
               SegState=double(K.X),
               K=as.integer(K.X),
               BaseAmp=as.double(BaseAmp.all[2]),
               sigma2=as.double(sigma2),
               T=as.double(T),PACKAGE="mad")                   


      aux.Y<-.C("RcallClassifySegments",
               SegLen=as.integer(Segments.Y$LenProbe),            
               SegAmp=as.double(Segments.Y$MeanAmp),
               SegState=double(K.Y),
               K=as.integer(K.Y),
               BaseAmp=as.double(BaseAmp.all[3]),
               sigma2=as.double(sigma2),
               T=as.double(T),PACKAGE="mad")                   
 
      Segments$State <- c(aux$SegState, aux.X$SegState, aux.Y$SegState)
    }

    else
     {
       Segments$State <- aux$SegState
     }





# Get genomic information
#
  if (!is.null(gen.info)) # JRG Oct'09
    {
      ini<-Segments[,1]
      end<-Segments[,2]
      Segments[,1]<-gen.info[ini,3]
      Segments[,2]<-gen.info[end,3]
    }
  else
   {
     ini<-Segments[,1]
     end<-Segments[,2]
   }



#
# Normalize Segments JRG  Oct'11
#

# Autosomes
  if (!is.null(gen.info))
    {
      o<-Segments[,5]%in%c(1:22)
      Segments[o,4]<-Segments[o,4]-BaseAmp.all[1]

if (T) # BaseAmp.all is not used to normalize, only to classify segments.
 {
  # Chr X,Y
   o<-Segments[,5]%in%c("X")
   Segments[o,4]<-Segments[o,4]-BaseAmp.all[2]
   o<-Segments[,5]%in%c("Y")
   Segments[o,4]<-Segments[o,4]-BaseAmp.all[3]
 } 

    }

  else

    {
      Segments[,4]<-Segments[,4]-BaseAmp.all[1]

    }





# rename
  names(Segments)[4]<-"qqBdev"
  Segments$qqBdev<-round(Segments$qqBdev,2)



#
# CLASSIFY SEGMENTS usinng %HetWt %Hom and LRR
#

  stateGADA<-Segments$State

  State2<-Segments$State
  State2[Segments$State!=0 & (Segments$"%HetWT" > Het.threshold | Segments$"Bdev" < 0.01)]<-0 

   
# replace State for the new one
 
  Segments$State<-abs(State2)

# classify segments
# 1: Mosaics 
# 5: LOH


  o<-Segments$State!=0 & (Segments$"%Hom">Hom.threshold & Segments$"Bdev" <= 0.50)
  Segments$State[o]<-5

  o<-Segments$State!=0 & (Segments$"%Hom"<Hom.threshold | Segments$"Bdev" >= 0.50 | (Segments$"%Hom">Hom.threshold & (Segments$LRR< -0.02  |  Segments$LRR>0.3)))
  Segments$State[o]<-1

  o<-Segments$LRR>0.3
  Segments$State[o]<-0



# classify Mosaics
#
#  1: UPD
#  2: deletion
#  3: duplication
#  4: trisomy (deleted - see version 0.9-3)
  

   o<-Segments$State==1 & Segments$LRR<=LRR.threshold[1]
   Segments$State[o]<-2

   o<-Segments$State==1 & Segments$LRR>=LRR.threshold[2]
   Segments$State[o]<-3


  if (print)
   {
    cat("------------------------------------------------------------- \n") 
    cat("Sparse Bayesian Learnig (SBL) algorithm (B-deviation version)\n")
    cat("Backward Elimination procedure with T=",T," and minimun length size=",MinSegLen, "\n",sep="")
    cat(" Number of segments = ",K,"\n")
    cat(" Base Amplitude of quantile normalized B-deviation: chr 1:22:", round(BaseAmp.all[1],4), " chr X:", round(BaseAmp.all[2],4), " chr Y:", round(BaseAmp.all[3],4),"\n", sep="")
    cat("------------------------------------------------------------- \n") 
    cat(" Mosaics (1), LOH (2)\n")
    cat("------------------------------------------------------------- \n") 
    print(Segments[Segments$State!=0,])

   }
  else
   {
    cat("---------------------------------------- \n") 
    cat("Sparse Bayesian Learnig (SBL) algorithm \n")
    cat("Backward Elimination procedure with T=",T," and minimun length size=",MinSegLen, "\n",sep="")
    cat(" Number of segments = ",K,"\n")
    cat(" Base Amplitude of copy number 2: chr 1:22:", round(BaseAmp.all[1],4), ", X=", round(BaseAmp.all[2],4), ", Y=",round(BaseAmp.all[3],4),"\n", sep="")
   }
  
 
 class(Segments)<-c("data.frame","summary.BackwardElimination")
 names(BaseAmp.all)<-c("Autosomes")
 attr(Segments,"BaseAmp")<-BaseAmp.all
 attr(Segments,"index")<-cbind(ini,end)
 if (is.null(gen.info))
  sigma2<-object$sigma2
 else
  sigma2<-object[[1]]$sigma2
 attr(Segments,"sigma2")<-sigma2
 attr(Segments,"stateGADA")<-stateGADA
 attr(Segments,"T")<-T
 attr(Segments,"df")<-MinSegLen-1
 invisible(Segments)

 }



