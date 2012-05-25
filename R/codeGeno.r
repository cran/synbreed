# coding genotypic data

codeGeno <- function(gpData,impute=FALSE,impute.type=c("random","family","beagle","beagleAfterFamily","fix"),replace.value=NULL,
                     maf=NULL,nmiss=NULL,label.heter="AB",keep.identical=TRUE,verbose=FALSE,minFam=5,showBeagleOutput=FALSE,tester=NULL,print.report=FALSE){

  #============================================================
  # read information from arguments
  #============================================================
  
  noHet <- is.null(label.heter)|!is.null(tester) # are there only homozygous genotypes?, we need this for random imputation
  if(is.null(impute.type)) impute.type <- "random"   # default
  if(impute) impute.type <- match.arg(impute.type)
  if (is.character(label.heter)) if(label.heter == "alleleCoding") label.heter <- function(x){substr(x, 1, 1) != substr(x, 3, 3)}

  orgFormat <- class(gpData)
  # check for class 'gpData'
  if(class(gpData)=="gpData"){
    if(is.null(gpData$geno)) stop("no genotypic data available") else res <- gpData$geno
    # family information (population structure) for genotypic data
    # drop unused levels
    if(is.factor(gpData$covar$family)) {
      popStruc <- droplevels(gpData$covar$family[gpData$covar$genotyped]) 
    } else {
      popStruc <- gpData$covar$family[gpData$covar$genotyped] 
    }
    if(gpData$info$codeGeno & !noHet){
       warning("assuming heterozygous genotypes coded as 1. Use 'label.heter' to specify if that is not the case")
       label.heter <- "1"
    } 
  } else { # atm other formats are supported too
    if(impute & impute.type %in% c("beagle","beagleAfterFamily")) stop("using Beagle is only possible for a gpData object")
    res <- gpData
    popStruc <- NULL
    gpData <- list(geno=res)
    gpData$map <- NULL
  }                                
  #  catch errors
  if(class(res)!= "data.frame" & class(res) != "matrix") stop("wrong data format")
  if (any(colMeans(is.na(res))==1)) warning("markers with only missing values in data")
  
  # number of genotypes
  n <- nrow(res)

  # keep names of data object
  cnames <- colnames(res)
  rnames <- rownames(res)
  res <- matrix(unlist(res),nrow=n)
  # tester control
  if(!is.null(tester)){
    if(length(tester)>1) stop("Only one tester is allowed for this function\n")
    if(!tester %in% rnames) stop("Tester has no genotype in your gpData-object\n")
  }

  # elements from control list
  # catch errors 
  if (impute){
    if(!is.logical(impute)) stop("impute has to be logical")
    if(impute.type=="fix" & is.null(replace.value)) stop("'replace.value' must be given for impute.type='fix'")
    # imputing with family information
    if(impute.type=="family" & is.null(popStruc)) stop(paste("family information needed, but '",substitute(gpData),"$covar$family' is empty",sep=""))
    if(impute.type=="family" & !is.null(popStruc)){
      if(length(popStruc)!=n) stop("population structure must have equal length as obsersvations in genotypic data")
      if(any(is.na(popStruc))) warning("missing values in family information, imputation is likely to be incomplete")
    }
  }

  
 
  #============================================================
  # step 1  - remove markers with more than nmiss fraction of missing values (optional, argument nmiss>0)
  #============================================================
 
  if(!is.null(nmiss)){
    if(nmiss<0 | nmiss>1) stop("'nmiss' must be in [0,1]")
    which.miss <- apply(is.na(res),2,mean,na.rm=TRUE)<=nmiss 
    res <- res[,which.miss]
    if (verbose) cat("step 1  :",sum(!which.miss),"marker(s) removed with >",nmiss*100,"% missing values \n")
    cnames <- cnames[which.miss]
    # update map
    if(!is.null(gpData$map)) gpData$map <- gpData$map[which.miss,]
  } else {
    res <- res
    if (verbose) cat("step 1  : No markers removed due to fraction of missing values \n")
  }
    
  #============================================================
  # step 2  - coding alleles
  #============================================================

  if (verbose) cat("step 2  : Recoding alleles \n")
    # identify heterozygous genotypes
    if(!is.null(label.heter)){
      if (is.character(label.heter)) {
        label.heter <- label.heter # 1 label for heterozygous
      } else {
        if (is.function(label.heter)){                          # multiple labels for heterozygous values
          is.heter <- label.heter
          label.heter <- unique(res[which(is.heter(res),arr.ind=TRUE)])
      } else stop("label.heter must be a character string or a function")
    } 
    # make sure that NA is not in label.heter
    # otherwise missing values would be masked
    label.heter <- label.heter[!is.na(label.heter)]
  }

  # inititialize report list
  if(print.report){
    alleles <- apply(res,2,table,useNA="no")
    major.allele <- function(x) names(which.max(x[!names(x) %in% label.heter]))
    minor.allele <- function(x) names(which.min(x[!names(x) %in% label.heter]))
   
    major <- unlist(apply(alleles,2,major.allele))
    minor <- unlist(apply(alleles,2,minor.allele))
    

    names(major) <- names(minor) <- cnames
  }
   
  #============================================================
  # step 2a  - Discarding markers for which the tester is not homozygous or values missing (optional, argument tester = "xxx")
  #============================================================

  if(!is.null(tester)){
    which.miss <- res[rnames==tester,]!=label.heter&!is.na(res[rnames==tester,])
    res <- res[,which.miss]
    cnames <- cnames[which.miss]
    if(sum(!which.miss) > 0){
      if (verbose) cat("step 2a :",sum(!which.miss),"marker(s) discarded because heterozygousity at tester locus or \n          missing values of the tester\n")
    } else {
      if (verbose) cat("step 2a : No marker(s) discarded because heterozygousity at tester locus or \n          missing values of the tester\n")
    }
    # update map
    if(!is.null(gpData$map)) gpData$map <- gpData$map[which.miss,]
  } 

  # function to recode alleles within one locus : 0 = major, 2 = minor 
  codeNumeric <- function(x){
    # names of alleles ordered by allele frequency
    alleles <-  names(table(x)[order(table(x),decreasing=TRUE)])
    # do not use heterozygous values
    alleles <- alleles[!alleles %in% label.heter]
    if (length(alleles)>2) stop("more than 2 marker genotypes found but no 'label.heter' declared")
    x[x %in% alleles] <- (as.numeric(factor(x[x %in% alleles],levels=alleles))-1)*2
    return(x)
  }

  # apply function on whole genotypic data
  res <- apply(as.matrix(res),2,codeNumeric)
 
  # set heterozygous genotypes as 1
  res[res %in% label.heter] <- 1
  res <- matrix(as.numeric(res),nrow=n)
               
  #============================================================
  # step 2b  - Discarding markers for which the tester has the minor allele
  #============================================================

  if(!is.null(tester)){
    which.miss <- res[rnames==tester,] != 2
    res <- res[,which.miss]
    cnames <- cnames[which.miss]
    if(sum(!which.miss) > 0){
      if (verbose) cat("step 2b :",sum(!which.miss),"marker(s) discarded for which the tester has the minor allele\n")
    } else{
      if (verbose) cat("step 2b : No marker(s) discarded for which the tester has the minor allele\n")
    }
    # update map
    if(!is.null(gpData$map)) gpData$map <- gpData$map[which.miss,]
  } 

  #============================================================
  # step 2c  - Discarding homozygout values of the minor allele and markers with more than nmiss values    
  #============================================================

  if(!is.null(tester)){
    res[res == 2] <- NA
    res <- matrix(as.numeric(res), nrow = n)
    if(!is.null(nmiss)){
      which.miss <- apply(is.na(res),2,mean,na.rm=TRUE) <= nmiss
      res <- res[,which.miss]
      cnames <- cnames[which.miss]
      if (verbose) cat("step 2c :",sum(!which.miss),"marker(s) discarded with >",nmiss*100,"% false genotyping values \n")
      # update map
      if(!is.null(gpData$map)) gpData$map <- gpData$map[which.miss,]
    } else{
      if (verbose) cat("step 2c : No markers discarded due to fraction of missing values \n")
    }
  }
  #============================================================
  # step 2.1 - discard duplicated markers   (optional, argument keep.identical=FALSE)
  #============================================================
  
  if(!keep.identical){
    nv <- 0
    if(!is.null(maf)){
      if(maf<0 | maf>1) stop("'maf' must be in [0,1]")
      nv <- ncol(res)
      part <- round(seq(0,nv, length.out=ceiling(nv/50000)+1))
      colnames(res) <- cnames
      for(i in length(part):2){
         res <- res[, part[i-1]+((part[i-1]+1):part[i])[colMeans(res, na.rm=TRUE)!=0]]
      }
      cnames <- colnames(res)
      if(!is.null(gpData$map)) gpData$map <- gpData$map[cnames, ]
      nv <- nv-ncol(res)
    }
    which.duplicated <- duplicated(res,MARGIN=2)
    res <- res[,!which.duplicated]
    cnames <- cnames[!which.duplicated]
    if(!is.null(gpData$map)) gpData$map <- gpData$map[!which.duplicated,]
    step1<- sum(which.duplicated)
    which.duplicated <- rep(FALSE, ncol(res))  
    which.miss <- apply(is.na(res),2,sum)>0
    which.miss <- (1:length(which.miss))[which.miss]
    if(length(which.miss[which.miss]) == ncol(res)) which.miss <- which.miss[1:(length(which.miss)-1)]
    for(i in which.miss){
      if(which.duplicated[i]) next  
      for(j in ((i+1):ncol(res))[!which.duplicated[(i+1):ncol(res)]])
        if(all(res[, i] == res[, j], na.rm = TRUE)){
          if(sum(is.na(res[, i])) >= sum(is.na(res[, j]))){
            which.duplicated[i] <- TRUE
            break
          } else {
            which.duplicated[j] <- TRUE
          }
        }
    }
    res <- res[,!which.duplicated]
    cnames <- cnames[!which.duplicated]
    if (verbose) cat("step 2.1:",sum(which.duplicated)+step1,"duplicated marker(s) removed \n")
    if (verbose) if(nv>0) cat("         and", nv, "monomorphic marker(s)\n")
    # update map 
    if(!is.null(gpData$map)) gpData$map <- gpData$map[!which.duplicated,]
  } else {
    if (verbose) cat("step 2.1: No duplicated markers discarded \n")
  }

  # coding of SNPs finished

  #============================================================
  # step 3  - imputing missing genotypes  (optional, argument impute=TRUE)
  #============================================================
  
  # initialize counter  
    cnt1 <- rep(0,ncol(res))   # for nr. of imputations with family structure
    cnt2 <- rep(0,ncol(res))    # for nr. of beagle imputations
    cnt3 <- rep(0,ncol(res))    # for nr. of random imputations
    names(cnt1) <- names(cnt2) <- names(cnt3) <- cnames
  
  # start of imputing
  if(impute){
    # number of markers
    M <- ncol(res)
    if(M==0) stop("no markers remained after step 1 (to many missing values)")          
  
    # number of missing values
    nmv <- sum(is.na(res))
 
    

    ###########################################################################
    # if impute.type="fix", replace missing values according to specified value
    ###########################################################################
    if(impute.type=="fix"){  
      res[is.na(res)] <- replace.value
      if (verbose) cat("step 3a : Replace missing values by",replace.value," \n")
    } 
    #########################################################
    # impute missing values according to population structure
    #########################################################
    if(impute.type %in% c("family" ,"beagleAfterFamily")){
      if (verbose) cat("step 3b : Imputing of missing values by family information \n")
      # initialize counter (- number of heterozygous values) 
      # loop over all markers
      probList <- list(c(1), c(.5,.5), c(.25,.5,.25))
      for (j in 1:M){
        if(sum(!is.na(res[,j]))>0){
          if(j==1) ptm <- proc.time()[3]
          try({# compute population structure  as counts
               poptab <- table(popStruc,res[,j])
               nFam <- table(popStruc)
               rS <- rowSums(poptab)     
         
               # continue only if there are missing values
               if(sum(is.na(res[,j]))>0 ){
                 # compute otherstatistics
                 major.allele <- unlist(attr(poptab,"dimnames")[[2]][apply(poptab,1,which.max)])
          
                 # look if SNP is segregating  for this population
                 polymorph <- apply(poptab,1,length) >1 & ( apply(poptab,1,min) != 0) 
                 polymorph2 <- apply(poptab,1,min) ==0  | apply(poptab,1,max) < minFam 
                 polymorph[!polymorph2] <- TRUE

                 # count missing vlalues
                 nmissfam <- tapply(is.na(res[,j]),popStruc,sum)
          
                 # must be a named list
                 names(major.allele) <- names(polymorph)
          
                 # loop over all families          
                 for ( i in rownames(poptab)[nmissfam>0] ){                          
                   # impute values for impute.type="family" : all missing genotypes
                   allTab <- table(res[popStruc == i, j])
                   if(length(allTab) == 0 & noHet) {allTab <- table(c(0,2))
                   } else if(all(names(allTab) == c(0, 2)) & !noHet)  allTab <- table(c(0,1,1,2))
                    if (impute.type=="family"){
                      res[is.na(res[,j]) & popStruc == i ,j] <- ifelse(length(allTab)>1,sample(as.numeric(names(allTab)),size=nmissfam[i],prob=probList[[length(allTab)]],replace=TRUE),as.numeric(names(allTab)))
                     # update counter
                     ifelse(polymorph[i],cnt3[j] <- cnt3[j] + nmissfam[i],cnt1[j] <- cnt1[j] + nmissfam[i])  
                   }
                 if (impute.type=="beagleAfterFamily"){
                   if (is.na(gpData$map$pos[j])){     # if no position is available use family algorithm
                     res[is.na(res[,j]) & popStruc == i ,j] <- ifelse(length(allTab)>1,sample(as.numeric(names(allTab)),size=nmissfam[i],prob=probList[[length(allTab)]],replace=TRUE),as.numeric(names(allTab)))
                     # update counter
                     ifelse(polymorph[i],cnt3[j] <- cnt3[j] +  nmissfam[i],cnt1[j] <- cnt1[j] +  nmissfam[i])  
                   } else{ # use Beagle and impute NA for polymorphic families
                     # impute values for impute.type="beagleAfterfamily"  : only monomorph markers
                     res[is.na(res[,j]) & popStruc == i ,j] <- as.numeric(ifelse(polymorph[i],NA,rep(major.allele[i],nmissfam[i])))
                     # update counter
                     ifelse(polymorph[i], cnt3[j] <- cnt3[j] +  0, cnt1[j] <- cnt1[j] + nmissfam[i]) 
                   } 
                 }  
               }
             }
             if(verbose) if(j==ceiling(M/100)) cat("          approximative run time ",(proc.time()[3] - ptm)*99," seconds ... \n",sep="")
          }) # end try
        }   # end of if(sum(!is.na(res[,j]))>0)
      }  # end of marker loop
    
    ###########################
    # run beagle for imputation
    ###########################
    }
    if(impute.type %in% c("beagle","beagleAfterFamily")){
      if (verbose) cat("step 3c : Imputing of missing values by Beagle \n")
      #if(any(grep(" ",.path.package()[grep("synbreed", .path.package())]))) warning("The package is installed in folder ",.path.package()[grep("synbreed", .path.package())]," which contains a space. To run beagle properly, please install the package to a differnt folder without spaces.")
      # use Beagle and impute NA for polymorphic families
      chr <- unique(gpData$map$chr)
      chr <- chr[!is.na(chr)]
      if(!is.null(tester)) 
        res <- res*2
      rownames(res) <- rownames(gpData$geno)
      colnames(res) <- rownames(gpData$map) 
      cnt2 <- apply(is.na(res),2,sum)
      # loop over chromosomses
      for (lg in seq(along=chr)){
        if(verbose) cat("          chromosome ", as.character(chr)[lg], "\n")
        sel <- rownames(gpData$map[is.na(gpData$map$pos) | gpData$map$chr != chr[lg],])
        if (length(sel)>0) {
           markerTEMPbeagle <- discard.markers(gpData,which=sel)
           markerTEMPbeagle$geno <- res[, colnames(res)[!colnames(res) %in% sel]]
        } else {
          markerTEMPbeagle <- gpData    # this occurs for only 1 chr
          markerTEMPbeagle$geno <- res
        }
        # recode for Beagle
        markerTEMPbeagle$geno[markerTEMPbeagle$geno==0] <- "AA"
        markerTEMPbeagle$geno[markerTEMPbeagle$geno==1] <- "AB"
        markerTEMPbeagle$geno[markerTEMPbeagle$geno==2] <- "BB"
     
        # update counter
        #cnt2 <- cnt2 + sum(is.na(markerTEMPbeagle$geno))
     
        # write input files for beagle
        pre <- paste("chr",chr[lg],sep="")
        # create new directory "beagle" for beagle input and output files
        if(!"beagle" %in% list.files()){
           if(.Platform$OS.type == "unix") system("mkdir beagle")
           if(.Platform$OS.type == "windows") shell("mkdir beagle")
        } 
        write.beagle(markerTEMPbeagle,file.path(getwd(),"beagle"),prefix=pre)
        output <- system(paste("java -Xmx1000m -jar ", shQuote(.path.package()[grep("synbreed", .path.package())][1]),     # caution with more than one pacakge with names synbreed*, assume synbreed to be the first one
                     "/exec/beagle.jar unphased=beagle/",pre,"input.bgl markers=beagle/",pre,"marker.txt missing=NA out=",sep=""),
                     intern=!showBeagleOutput)
        if(.Platform$OS.type == "unix") system(paste("gzip -d -f beagle/",pre,"input.bgl.dose.gz",sep=""))
        if(.Platform$OS.type == "windows") shell(paste("gzip -d -f beagle/",pre,"input.bgl.dose.gz",sep=""))

        # read data from beagle
        resTEMP <- read.table(paste("beagle/",pre,"input.bgl.dose",sep=""),header=TRUE,row.names=1)
        resTEMP <- t(resTEMP[,-c(1:2)])
      
        # convert dose to genotypes
        if(noHet){
          resTEMP[resTEMP<1] <- 0
          resTEMP[resTEMP>=1] <- 2
        } else {
          resTEMP <- round(resTEMP,0) # 0, 1, and 2
        }
      
        if (length(sel)>0) {
          res[,!colnames(res) %in% sel] <- resTEMP
        } else {
          res <- resTEMP
        }
      }  
    } 
    #########################################################################
    # impute missing values with no population structure or missing positions
    #########################################################################
    if(impute.type %in% c("random", "beagle")){
      if (verbose) cat("step 3d : Random imputing of missing values \n")
      # initialize counter (- number of heterozygous values) 
      for (j in 1:M){
        cnt3[j] <-  sum(is.na(res[,j]))
        # estimation of running time after the first iteration
        if(j==1) ptm <- proc.time()[3]
          p <- mean(res[,j],na.rm=TRUE)/2  # minor allele frequency
          if(noHet){        # assuming only 2 homozygous genotypes
            res[is.na(res[,j]),j] <- sample(c(0,2),size=sum(is.na(res[,j])),prob=c(1-p,p),replace=TRUE)
          } else {                            # assuming 3 genotypes
            res[is.na(res[,j]),j] <- sample(c(0,1,2),size=sum(is.na(res[,j])),prob=c((1-p)^2,2*p*(1-p),p^2),replace=TRUE)
          }
          if(j==ceiling(M/100) & verbose) cat("         approximate run time ",(proc.time()[3] - ptm)*99," seconds \n",sep=" ")
        } 
        # update counter for Beagle, remove those counts which where imputed ranomly 
        cnt2 <- cnt2-cnt3  
      }
      if(!is.null(tester) & impute.type %in% c("random","beagle", "beagleAfterFamily")) res <- res/2
    
  
    #============================================================
    # step 4 - recoding
    #============================================================
  
    # recode again if allele frequeny changed to to imputing
    if(any(colMeans(res,na.rm=TRUE)>1)){
      if (verbose) cat("step 4  : Recode alleles due to imputation \n")
      res[,which(colMeans(res,na.rm=TRUE)>1)] <- 2 - res[,which(colMeans(res,na.rm=TRUE)>1)]     
    } else{
      if (verbose) cat("step 4  : No recoding of alleles necessary after imputation \n") 
    }
  }
  # update report list
  if(print.report){
    major[which(colMeans(res,na.rm=TRUE)>1)] <- minor[which(colMeans(res,na.rm=TRUE)>1)]
  }
  
  
  #============================================================
  # step 5 - remove markers with minor allele frequency < maf  (optional, argument maf>0)
  #============================================================
  
  if(!is.null(maf)){
    if(maf<0 | maf>1) stop("'maf' must be in [0,1]")
    if(is.null(tester)) which.maf <- colMeans(res,na.rm=TRUE)>=2*maf else
      which.maf <- colMeans(res,na.rm=TRUE)>=maf & colMeans(res,na.rm=TRUE)<=1-maf
    if (verbose) cat("step 5  :",sum(!which.maf),"marker(s) removed with maf <",maf,"\n")
    res <- res[,which.maf]
    cnames <- cnames[which.maf] 
    # update map
    if(!is.null(gpData$map)) gpData$map <- gpData$map[which.maf,]
     # update report list

    
  } else {
    if (verbose) cat("step 5  : No markers discarded due to minor allele frequency \n")
  }
  
 
  
  #============================================================
  # step 6 - discard duplicated markers   (optional, argument keep.identical=FALSE)
  #============================================================
  
  if(!keep.identical){
       which.duplicated <- duplicated(res,MARGIN=2)
       res <- res[,!which.duplicated]
       cnames <- cnames[!which.duplicated]
       if (verbose) cat("step 6  :",sum(which.duplicated),"duplicated marker(s) removed \n")
       # update map 
       if(!is.null(gpData$map)) gpData$map <- gpData$map[!which.duplicated,]
   # update report list

       
  } else{
    if (verbose) cat("step 6  : No duplicated markers discarded \n")
  }
     
  #============================================================
  # step 6a - discard markers for which only the tester is different
  #============================================================
  
  if(!is.null(tester)){
    which.fixed <- apply(res, 2, sum) == nrow(res)-1
    res <- res[,!which.fixed]
    cnames <- cnames[!which.fixed]
    if(!is.null(gpData$map)) gpData$map <- gpData$map[!which.fixed,]
    if (verbose)
      if(sum(which.fixed) != 0){
        cat("step 6a :",sum(which.fixed),"in crosses fixed marker(s) removed \n")
      } else {
        cat("step 6a : No in crosses fixed marker(s) removed \n")
      }
  }
  
  #============================================================
  # step 7 - restoring original data format
  #============================================================
  
  #if (verbose) cat("step 7 : Restoring original data format \n")
  rownames(res) <- rnames
  colnames(res) <- cnames
  if(orgFormat == "matrix"){
    res <- matrix(res,nrow=n)
  }
  if(orgFormat == "data.frame"){
    res  <- as.data.frame(res)
  }

  if (verbose) cat("End      :",ncol(res),"marker(s) remain after the check\n")

  #============================================================
  # print summary of imputation
  #============================================================
  
  if(impute){
    cat("\n")
    cat("Summary of imputation \n")
    cat(paste("  total number of missing values                :",nmv,"\n"))
    if(impute.type %in% c("family","beagleAfterFamily"))                   cat(paste("  number of imputations by family structure     :",sum(cnt1),"\n"))
    if(impute.type %in% c("beagle","beagleAfterFamily"))                   cat(paste("  number of Beagle imputations                  :",sum(cnt2),"\n"))
    if(impute.type %in% c("beagle","random","family","beagleAfterFamily")) cat(paste("  number of random imputations                  :",sum(cnt3),"\n"))
  }
  
  # overwrite original genotypic data
  if(orgFormat == "gpData") {
    gpData$geno <- res
    gpData$info$codeGeno <- TRUE
  } else gpData <- res
  
  if(print.report){
    if (verbose) cat("  Writing report to file 'SNPreport.txt' \n")
    report.list <- data.frame(SNPname=cnames,major=major[cnames],minor=minor[cnames],MAF=round(colMeans(res)/2,3),impute.fam=cnt1[cnames],impute.beagle=cnt2[cnames],impute.ran=cnt3[cnames])
    write.table(report.list,file="SNPreport.txt",quote=FALSE,row.names=FALSE)
   } 

  # return a gpData object (or a matrix)
  return(gpData)
}
