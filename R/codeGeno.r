# coding genotypic data

codeGeno <- function(gpData,impute=FALSE,impute.type=c("random","family","beagle","beagleAfterFamily","fix"),replace.value=NULL,
                     maf=NULL,nmiss=NULL,label.heter="AB",reference.allele="minor",#keep.list=NULL,
                     keep.identical=TRUE,verbose=FALSE,minFam=5,showBeagleOutput=FALSE,tester=NULL,print.report=FALSE,check=FALSE){

  #============================================================
  # read information from arguments
  ##### rownames(res)[apply(is.na(res), 1, mean)>.5]
  #============================================================

  SEED <- round(runif(2,1,1000000),0)
  noHet <- is.null(label.heter)|!is.null(tester) # are there only homozygous genotypes?, we need this for random imputation

  if(is.null(impute.type)) impute.type <- "random"   # default
  MG <- rownames(gpData$geno)[apply(is.na(gpData$geno),1,all)]
  if(check)
    if(impute) {
      impute.type <- match.arg(impute.type)
      if(impute.type %in% c("beagle", "beagleAfterFamily"))
        if(!"beaglePath" %in% ls()){
          cat("Due to the policy of R, we removed the executable of beagle\n",
              "from our package. If you like to use the beagle imputation,\n",
              "you can put the executable beagle.jar somewhere. Then define\n",
              "with 'beaglePath' a variable with the location of beagle.jar.\n",
              "beagle.jar can be found in synbreed versions up to 0.9 in\n",
              "the synbreed subfolder exec")
          return(NULL)
          beaglePath <- NULL
        }
      if(impute.type %in% c("beagle", "beagleAfterFamily") & !is.null(gpData$map))
        if(grepl("string mismatches", all.equal(rownames(gpData$map), colnames(gpData$geno))))
          stop("Order of markers in geno and map does not fit!")
      if(impute.type %in% c("beagle") & length(MG) > 0) stop(paste("Of genotype(s) ", MG, " all genotypic values are missing!", sep=" "))
      else if(length(MG) > 0) warning(paste("Of genotype(s) ", MG, " all genotypic values are missing! \nImputation may be erroneus.", sep=" "))
    } else if(length(MG) > 0) warning(paste("Of genotype(s) ", MG, " all genotypic values are missing!", sep=" "))

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
  if(check){
    if(class(res)!= "data.frame" & class(res) != "matrix") stop("wrong data format")
    if(any(colMeans(is.na(res))==1)) warning("markers with only missing values in data")
    if(length(reference.allele)>1 & length(reference.allele)!=ncol(res)) stop("'reference allele' should be of length 1 or match the number of markers")
    if(class(reference.allele) != "character") stop("'reference allele' should be of class character")
  }
  # number of genotypes
  n <- nrow(res)

  # keep names of data object
  cnames <- colnames(res)
  rnames <- rownames(res)
  res <- matrix(unlist(res),nrow=n)
  # tester control
  if(!is.null(tester)){
    if(length(tester)>1) stop("Only one tester is allowed for this function\n")
    if(!tester %in% rnames) stop("Tester has no genotype in the gpData-object\n")
  }

  # elements from control list
  # catch errors
  if (impute){
    if(!is.logical(impute)) stop("impute has to be logical")
    if(impute.type=="fix" & is.null(replace.value)) stop("'replace.value' must be given for impute.type='fix'")
    # imputing with family information
    if((impute.type=="family" | impute.type=="beagleAfterFamily") & is.null(popStruc)) stop(paste("family information needed, but '",substitute(gpData),"$covar$family' is empty",sep=""))
    if((impute.type=="family" | impute.type=="beagleAfterFamily") & !is.null(popStruc)){
      if(any(is.na(popStruc))) warning("missing values in family information, imputation is likely to be incomplete")
      if(length(popStruc)!=n) stop("population structure must have equal length as obsersvations in genotypic data")
    }
  }

  # use same reference allele for all markers if not specified differently
  if(reference.allele[1]!="minor" & length(reference.allele)==1)     reference.allele <- rep(reference.allele,ncol(res))


  #============================================================
  # step 1  - remove markers with more than nmiss fraction of missing values (optional, argument nmiss>0)
  #============================================================

  if(!is.null(nmiss)){
    if(nmiss<0 | nmiss>1) stop("'nmiss' have to be in [0,1]")
    which.miss <- apply(is.na(res),2,mean,na.rm=TRUE)<=nmiss
    res <- res[,which.miss]
    if(!(reference.allele[1]=="minor" | reference.allele[1]=="keep"))  reference.allele <- reference.allele[which.miss]
    if (verbose) cat("   step 1  :", sum(!which.miss),"marker(s) removed with >",nmiss*100,"% missing values \n")
    cnames <- cnames[which.miss]
    # update map
    gpData$geno <- gpData$geno[, which.miss]
    if(!is.null(gpData$map)) gpData$map <- gpData$map[which.miss,]
    rm(which.miss)
  } else {
    if (verbose) cat("   step 1  : No markers removed due to fraction of missing values \n")
  }

  #============================================================
  # step 2  - coding alleles
  #============================================================

  if (verbose) cat("   step 2  : Recoding alleles \n")
  if(gpData$info$codeGeno) {
    if(reference.allele[1]=="minor"){
      afCols <- (1:ncol(res))[colMeans(res, na.rm=TRUE)>1]
      res[, afCols] <-  rep(1, nrow(res)) %*% t(rep(2, length(afCols))) - res[, afCols]
      # inititialize report list
      if(print.report){
        alleles <- apply(res,2,table,useNA="no")
        major.allele <- function(x) names(which.max(x[!names(x) %in% label.heter]))
        minor.allele <- function(x) names(which.min(x[!names(x) %in% label.heter]))
        if(class(alleles) != "list"){
          mytable <- function(x,...) as.data.frame(table(x,...),stringsAsFactors=FALSE)
          alleles <- apply(res,2,mytable,useNA="no")
          major.allele <- function(x) x$x[which.max(x$Freq)]
          minor.allele <- function(x) x$x[which.min(x$Freq)]
        }
        major <- unlist(sapply(alleles,major.allele))
        minor <- unlist(sapply(alleles,minor.allele))
        names(major) <- names(minor) <- cnames
      }

    } else {
      if(reference.allele[1]!="keep") gpData$info$codeGeno==FALSE
    }
  }
  if(!gpData$info$codeGeno) {
    if(reference.allele[1]=="minor"){

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
        if(class(alleles) != "list"){
          mytable <- function(x,...) as.data.frame(table(x,...),stringsAsFactors=FALSE)
          alleles <- apply(res,2,mytable,useNA="no")
          major.allele <- function(x) x$x[which.max(x$Freq)]
          minor.allele <- function(x) x$x[which.min(x$Freq)]
        }
        major <- unlist(sapply(alleles,major.allele))
        minor <- unlist(sapply(alleles,minor.allele))
        names(major) <- names(minor) <- cnames
      }
      # function to recode alleles within one locus : 0 = major, 2 = minor
      codeNumeric <- function(x){
        # names of alleles ordered by allele frequency
        alleles <-  names(table(x)[order(table(x),decreasing=TRUE)])
        # do not use heterozygous values
        alleles <- alleles[!alleles %in% label.heter]
        if (length(alleles)>2) stop("more than 2 marker genotypes found, but 'label.heter' is not declared")
        x[x %in% alleles] <- (as.numeric(factor(x[x %in% alleles],levels=alleles))-1)*2
        return(x)
      }

      # apply function on whole genotypic data
      res <- apply(as.matrix(res),2,codeNumeric)

      # set heterozygous genotypes as 1
      res[res %in% label.heter] <- 1
      res <- matrix(as.numeric(res),nrow=n)
    } else { #if reference.alle !="minor"
      count.ref.alleles <- function(x,ref){
   	    sum(x==ref)		
      }
      recode.by.ref.alleles <- function(x){
	    ref <- x[1]
	    x <- x[-1]
	    x2 <- ifelse(!is.na(x),strsplit(x,split=""),NA) # split genotype into both alleles
	    nr.ref.alleles <- unlist(lapply(x2,count.ref.alleles,ref))
	    return(nr.ref.alleles)
      }
      res_ref <- rbind(reference.allele,res)
      get.nonref.allele <- function(x){
        xx <- unlist(strsplit(x[-1],split=""))
        unique(xx)[unique(xx) != x[1] & !is.na(unique(xx))]
      }
      res <- apply(res_ref,2,recode.by.ref.alleles)

      if(print.report){
        major <- reference.allele
        minor <-  apply(res_ref,2,get.nonref.allele)
      }
    }
  }

  #============================================================
  # step 3  - Discarding markers for which the tester is not homozygous or values missing (optional, argument tester = "xxx")
  #============================================================

  if(!is.null(tester)){
    which.miss <- res[rnames==tester,]!=label.heter&!is.na(res[rnames==tester,])
    res <- res[,which.miss]
    cnames <- cnames[which.miss]
    if(sum(!which.miss) > 0){
      if (verbose) cat("   step 3 :",sum(!which.miss),"marker(s) discarded because heterozygousity at tester locus or \n          missing values of the tester\n")
    } else {
      if (verbose) cat("   step 3 : No marker(s) discarded because heterozygousity at tester locus or \n          missing values of the tester\n")
    }
    # update map
    if(!is.null(gpData$map)) gpData$map <- gpData$map[which.miss,]
  }


  #============================================================
  # step 4 - remove markers with minor allele frequency < maf  (optional, argument maf>0)
  #============================================================

  if(!is.null(maf)){
    if(maf<0 | maf>1) stop("'maf' must be in [0,1]")
    if(is.null(tester)){
      which.maf <- colMeans(res,na.rm=TRUE)>=2*maf
    } else {
      which.maf <- colMeans(res,na.rm=TRUE)>=maf & colMeans(res,na.rm=TRUE)<=1-maf
    }
    if (verbose) cat("   step 4  :",sum(!which.maf),"marker(s) removed with maf <",maf,"\n")
    res <- res[,which.maf]
    cnames <- cnames[which.maf]
    # update map
    if(!is.null(gpData$map)) gpData$map <- gpData$map[which.maf,]
     # update report list


  } else {
    if (verbose) cat("   step 4  : No markers discarded due to minor allele frequency \n")
  }
  #============================================================
  # step 5  - Discarding markers for which the tester has the minor allele
  #============================================================

  if(!is.null(tester)){
    which.miss <- res[rnames==tester,] != 2
    res <- res[,which.miss]
    cnames <- cnames[which.miss]
    if(sum(!which.miss) > 0){
      if (verbose) cat("   step 5  :",sum(!which.miss),"marker(s) discarded for which the tester has the minor allele\n")
    } else{
      if (verbose) cat("   step 5  : No marker(s) discarded for which the tester has the minor allele\n")
    }
    # update map
    if(!is.null(gpData$map)) gpData$map <- gpData$map[which.miss,]
  }

  #============================================================
  # step 6  - Discarding homozygout values of the minor allele and markers with more than nmiss values
  #============================================================

  if(!is.null(tester)){
    res[res == 2] <- NA
    res <- matrix(as.numeric(res), nrow = n)
    if(!is.null(nmiss)){
      which.miss <- apply(is.na(res),2,mean,na.rm=TRUE) <= nmiss
      res <- res[,which.miss]
      cnames <- cnames[which.miss]
      if (verbose) cat("   step 6  :",sum(!which.miss),"marker(s) discarded with >",nmiss*100,"% false genotyping values \n")
      # update map
      if(!is.null(gpData$map)) gpData$map <- gpData$map[which.miss,]
    } else{
      if (verbose) cat("   step 6  : No markers discarded due to fraction of missing values \n")
    }
  }

  #============================================================
  # step 7  - imputing missing genotypes  (optional, argument impute=TRUE)
  #============================================================

  # initialize counter
    cnt1 <- rep(0,ncol(res))   # for nr. of imputations with family structure
    cnt2 <- rep(0,ncol(res))    # for nr. of beagle imputations
    cnt3 <- rep(0,ncol(res))    # for nr. of random imputations
    names(cnt1) <- names(cnt2) <- names(cnt3) <- cnames

  # start of imputing
  if(impute){
    set.seed(SEED[1])
    # number of markers
    M <- ncol(res)
    if(M==0) stop(" no markers remained after step 1 (to many missing values)")
    if (verbose) cat("   step 7  : Imputing of missing values \n")

    # number of missing values
    nmv <- sum(is.na(res))

    ###########################################################################
    # if impute.type="fix", replace missing values according to specified value
    ###########################################################################
    if(impute.type=="fix"){
      res[is.na(res)] <- replace.value
      if (verbose) cat("   step 7a : Replace missing values by",replace.value," \n")
    }
    #########################################################
    # impute missing values according to population structure
    #########################################################
    if(impute.type %in% c("family" ,"beagleAfterFamily")){
      if (verbose) cat("   step 7b : Imputing of missing values by family information \n")
      # initialize counter (- number of heterozygous values)
      # loop over all markers
      probList <- list(c(1), c(.5,.5), c(.25,.5,.25))
      vec.cols <- (1:M)[is.na(colSums(res, na.rm = FALSE))]
      for (j in vec.cols){
        if(j==vec.cols[1]) ptm <- proc.time()[3]
        if(sum(!is.na(res[,j]))>0){
          try({# compute population structure  as counts
               poptab <- table(popStruc,res[,j])
               nFam <- table(popStruc)
               rS <- rowSums(poptab)
               # compute otherstatistics
               major.allele <- unlist(attr(poptab,"dimnames")[[2]][apply(poptab,1,which.max)])
               # look if SNP is segregating  for this population
               polymorph <- apply(poptab,1,length) > 1 & (apply(poptab,1,min) != 0)
               polymorph2 <- rS > minFam
               polymorph[!polymorph2] <- TRUE
               # count missing values
               nmissfam <- tapply(is.na(res[,j]),popStruc,sum)
               # must be a named list
               names(major.allele) <- names(polymorph)
               # loop over all families
               for (i in rownames(poptab)[nmissfam > 0]){
                 # impute values for impute.type="family" : all missing genotypes
                 allTab <- table(res[popStruc == i, j])
                 if(length(allTab) == 0 & noHet) {
                   allTab <- table(c(0,2))
                 } else if(all(names(allTab) == c(0, 2)) & !noHet)  allTab <- table(c(0,1,1,2))
                 if (impute.type=="family"){
                   res[is.na(res[,j]) & popStruc == i ,j] <- ifelse(length(allTab)>1, sample(as.numeric(names(allTab)),size=nmissfam[i],prob=probList[[length(allTab)]],replace=TRUE),as.numeric(names(allTab)))
                   # update counter
                   if(polymorph[i]) cnt3[j] <- cnt3[j] + nmissfam[i] else cnt1[j] <- cnt1[j] + nmissfam[i]
                 }
                 if(impute.type=="beagleAfterFamily"){
                   if (is.na(gpData$map$pos[j])){     # if no position is available use family algorithm
                     res[is.na(res[,j]) & popStruc == i ,j] <- ifelse(length(allTab)>1, sample(as.numeric(names(allTab)),size=nmissfam[i],prob=probList[[length(allTab)]],replace=TRUE),as.numeric(names(allTab)))
                     # update counter
                     if(polymorph[i]) cnt3[j] <- cnt3[j] +  nmissfam[i] else cnt1[j] <- cnt1[j] +  nmissfam[i]
                   } else { # use Beagle and impute NA for polymorphic families
                     if(!polymorph[i]){
                       # impute values for impute.type="beagleAfterfamily"  : only monomorph markers
                       res[is.na(res[,j]) & popStruc == i ,j] <- as.numeric(rep(major.allele[i],nmissfam[i]))
                       # update counter
                       cnt1[j] <- cnt1[j] + nmissfam[i]
                     }
                   }
                 }
               }
               if(j==ceiling(length(vec.cols)/100)) if(verbose) cat("         approximative run time for imputation by family information ",(proc.time()[3] - ptm)*99," seconds ... \n",sep="")
          }) # end try
        }   # end of if(sum(!is.na(res[,j]))>0)
      } # end of marker loop
    }

    ###########################
    # run beagle for imputation
    ###########################
    if(impute.type %in% c("beagle","beagleAfterFamily")){
      if (verbose) cat("   step 7c : Imputing of missing values by Beagle \n")
        #if(any(grep(" ",path.package()[grep("synbreed", path.package())]))) warning("The package is installed in folder ",path.package()[grep("synbreed", path.package())]," which contains a space. To run beagle properly, please install the package to a differnt folder without spaces.")
        # use Beagle and impute NA for polymorphic families
        chr <-  unique(gpData$map$chr)
        chr <- chr[!is.na(chr)]
        if(!is.null(tester))
          res <- res*2
        rownames(res) <- rownames(gpData$geno)
        colnames(res) <- cnames
        cnt2 <- apply(is.na(res),2,sum)
        # loop over chromosomses
        for (lg in seq(along=chr)){
          if(lg==1) ptm <- proc.time()[3]
          if(verbose) cat("          chromosome ", as.character(chr)[lg], "\n")
          sel <- rownames(gpData$map[is.na(gpData$map$pos) | gpData$map$chr != chr[lg] | !rownames(gpData$map) %in% colnames(res) ,])
          if (length(sel)>0) {
            markerTEMPbeagle <- discard.markers(gpData,which=sel)
            markerTEMPbeagle$geno <- res[, colnames(res)[!colnames(res) %in% sel]]
          } else {
            markerTEMPbeagle <- gpData    # this occurs for only 1 chr
          }
          # recode for Beagle
          markerTEMPbeagle$geno[markerTEMPbeagle$geno==0] <- "AA"
          markerTEMPbeagle$geno[markerTEMPbeagle$geno==1] <- "AB"
          markerTEMPbeagle$geno[markerTEMPbeagle$geno==2] <- "BB"

          # update counter
          #cnt2 <- cnt2 + sum(is.na(markerTEMPbeagle$geno))

          # write input files for beagle
          pre <- paste("chr",chr[lg],sep="")
          pre <- gsub(" ", "_", pre, fixed=TRUE)
          # create new directory "beagle" for beagle input and output files
          if(!"beagle" %in% list.files()){
            if(.Platform$OS.type == "unix")    system("mkdir beagle")
            if(.Platform$OS.type == "windows") shell("mkdir beagle")
          } else if(lg == 1 & length(list.files("beagle")) > 0 ) {
            if(.Platform$OS.type == "unix")    system("rm beagle/*.*")
            if(.Platform$OS.type == "windows") shell("rm beagle/*.*")
          }
          write.beagle(markerTEMPbeagle,file.path(getwd(),"beagle"),prefix=pre)
          output <- system(paste("java -Xmx1000m -jar ",
          # sort(path.package()[grep("synbreed", path.package())])[1]
                           shQuote(paste(beaglePath, "/beagle.jar", sep="")),
                           # caution with more than one pacakge with names synbreed*, assume synbreed to be the first one
                           " unphased=beagle/", pre, "input.bgl markers=beagle/", pre, "marker.txt missing=NA out=", sep=""),
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
          if(lg==1) if(verbose) cat("\n         approximative run time for beagle imputation ",
                                        (proc.time()[3] - ptm)/ncol(markerTEMPbeagle$geno)*ncol(res)," seconds ... \n\n",sep="")
        }
      }

      #########################################################################
      # impute missing values with no population structure or missing positions
      #########################################################################
      if(impute.type %in% c("random", "beagle", "beagleAfterFamily")){
        if (verbose) cat("   step 7d : Random imputing of missing values \n")
        # initialize counter (- number of heterozygous values)
        for (j in (1:M)[apply(is.na(res), 2, sum)>0]){
          cnt3[j] <-  cnt3[j] + sum(is.na(res[,j]))
          # estimation of running time after the first iteration
          if(j==1) ptm <- proc.time()[3]
            p <- mean(res[,j],na.rm=TRUE)/2  # minor allele frequency
            if(noHet){        # assuming only 2 homozygous genotypes
              res[is.na(res[,j]),j] <- sample(c(0,2),size=sum(is.na(res[,j])),prob=c(1-p,p),replace=TRUE)
            } else {                            # assuming 3 genotypes
              res[is.na(res[,j]),j] <- sample(c(0,1,2),size=sum(is.na(res[,j])),prob=c((1-p)^2,2*p*(1-p),p^2),replace=TRUE)
            }
            if(j==ceiling(M/100)) if(verbose) cat("         approximate run time for random imputation ",(proc.time()[3] - ptm)*99," seconds \n",sep=" ")
          }
          # update counter for Beagle, remove those counts which where imputed ranomly
          if(impute.type == "beagle") cnt2 <- cnt2-cnt3
        }
        if(!is.null(tester) & impute.type %in% c("random","beagle", "beagleAfterFamily")) res <- res/2


    #============================================================
    # step 8 - recoding
    #============================================================

    # recode again if allele frequeny changed to to imputing
    if(any(colMeans(res,na.rm=TRUE)>1)){
      if (verbose) cat("   step 8  : Recode alleles due to imputation \n")
      res[,which(colMeans(res,na.rm=TRUE)>1)] <- 2 - res[,which(colMeans(res,na.rm=TRUE)>1)]
    } else{
      if (verbose) cat("   step 8  : No recoding of alleles necessary after imputation \n")
    }
  }
  # update report list
  if(print.report){
    major[which(colMeans(res,na.rm=TRUE)>1)] <- minor[which(colMeans(res,na.rm=TRUE)>1)]
  }

  #============================================================
  # step 9 - remove markers with minor allele frequency < maf  (optional, argument maf>0)
  #============================================================

  if(!is.null(maf) & impute){
    if(maf<0 | maf>1) stop("'maf' must be in [0,1]")
    if(is.null(tester)) which.maf <- colMeans(res,na.rm=TRUE)>=2*maf else
      which.maf <- colMeans(res,na.rm=TRUE)>=maf & colMeans(res,na.rm=TRUE)<=1-maf
    if (verbose) cat("   step 9  :",sum(!which.maf),"marker(s) removed with maf <",maf,"\n")
    res <- res[,which.maf]
    cnames <- cnames[which.maf]
    # update map
    if(!is.null(gpData$map)) gpData$map <- gpData$map[which.maf,]
     # update report list


  } else {
    if (verbose & impute) cat("   step 9  : No markers discarded due to minor allele frequency \n")
  }

  #============================================================
  # step 10 - discard duplicated markers   (optional, argument keep.identical=FALSE)
  #============================================================

  if(!keep.identical){
    set.seed(SEED[2])
    cnms <- sample(1:ncol(res))
    which.duplicated <- duplicated(res[, cnms],MARGIN=2)
    cnms[!which.duplicated]
    res <- res[, sort(cnms[!which.duplicated])]
    cnames <- cnames[sort(cnms[!which.duplicated])]
    if (verbose) cat("   step 10 :",sum(which.duplicated),"duplicated marker(s) removed \n")
    # update map
    if(!is.null(gpData$map)) gpData$map <- gpData$map[sort(cnms[!which.duplicated]),]
   # update report list


  } else{
    if (verbose) cat("   step 10 : No duplicated markers removed \n")
  }

  #============================================================
  # step 10a - discard markers for which only the tester is different
  #============================================================

  if(!is.null(tester)){
    which.fixed <- apply(res, 2, sum) == nrow(res)-1
    res <- res[,!which.fixed]
    cnames <- cnames[!which.fixed]
    if(!is.null(gpData$map)) gpData$map <- gpData$map[!which.fixed,]
    if (verbose)
      if(sum(which.fixed) != 0){
        cat("   step 10a:",sum(which.fixed),"in crosses fixed marker(s) removed \n")
      } else {
        cat("   step 10a: No in crosses fixed marker(s) removed \n")
      }
  }

  #============================================================
  # step 11 - restoring original data format
  #============================================================

  rownames(res) <- rnames
  colnames(res) <- cnames
  if(orgFormat == "matrix"){
    res <- matrix(res,nrow=n)
  }
  if(orgFormat == "data.frame"){
    res  <- as.data.frame(res)
  }

  if (verbose) cat("   End     :",ncol(res),"marker(s) remain after the check\n")

  #============================================================
  # print summary of imputation
  #============================================================

  if(impute){
    cat("\n")
    cat("     Summary of imputation \n")
    cat(paste("    total number of missing values                :",nmv,"\n"))
    if(impute.type %in% c("family","beagleAfterFamily"))                   cat(paste("    number of imputations by family structure     :",sum(cnt1),"\n"))
    if(impute.type %in% c("beagle","beagleAfterFamily"))                   cat(paste("    number of Beagle imputations                  :",sum(cnt2),"\n"))
    if(impute.type %in% c("beagle","random","family","beagleAfterFamily")) cat(paste("    number of random imputations                  :",sum(cnt3),"\n"))
  }

  # overwrite original genotypic data
  if(orgFormat == "gpData") {
    gpData$geno <- res
    gpData$info$codeGeno <- TRUE
  } else gpData <- res

  if(print.report){
    if (verbose) cat("  Writing report to file 'SNPreport.txt' \n")
    report.list <- data.frame(SNPname=cnames,major=major[cnames],minor=minor[cnames],MAF=round(colMeans(res,na.rm=TRUE)/2,3),impute.fam=cnt1[cnames],impute.beagle=cnt2[cnames],impute.ran=cnt3[cnames])
    write.table(report.list,file="SNPreport.txt",quote=FALSE,row.names=FALSE)
   }

  # return a gpData object (or a matrix)
  return(gpData)
}
