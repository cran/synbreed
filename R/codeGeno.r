# coding genotypic data

codeGeno <- function(gpData,impute=FALSE,impute.type=c("random","family","beagle","beagleAfterFamily","beagleNoRand","beagleAfterFamilyNoRand","fix"),replace.value=NULL,
                     maf=NULL,nmiss=NULL,label.heter="alleleCoding",reference.allele="minor",keep.list=NULL,
                     keep.identical=TRUE,verbose=FALSE,minFam=5,showBeagleOutput=FALSE,tester=NULL,print.report=FALSE,
                     check=FALSE,ploidy=2,cores=1){

  #============================================================
  # read information from arguments
  ##### rownames(res)[apply(is.na(res), 1, mean)>.5]
  #============================================================

  impute.type <- match.arg(impute.type)
  infoCall <- match.call()
  SEED <- round(runif(2,1,1000000),0)
  noHet <- is.null(label.heter)|!is.null(tester) # are there only homozygous genotypes?, we need this for random imputation
  multiLapply <- function(x,y,...,mc.cores=1){
    if(.Platform$OS.type == "windows" & cores>1){
      cl <- makeCluster(min(mc.cores, detectCores()))
      registerDoParallel(cl)
      parLapply(cl,x,y,...)
      stopCluster(cl)
    } else {
      mclapply(x,y,...,mc.cores=mc.cores)
    }
  }
  multiCor <- function(x, y=NULL, use="everything", method = c("pearson", "kendall", "spearman"), cores=1){
    method <- match.arg(method)
    if(cores==1) cor(x, y, use=use, method=method) else if(is.null(y)){
      ncolX <- ncol(x); namesX <- colnames(x)
      mat <- matrix(NA, nrow=ncolX, ncol=ncolX)
      nBlocks <- 1:50; nBlocks <- which.min((nBlocks*(nBlocks+1)*.5)%%cores)
      if(nBlocks==1) {nBlocks <- 1:50; nBlocks <- which.min((nBlocks*(nBlocks+1)*.5)%%(cores-1))}
      blocks <- unique(t(apply(cbind(rep(1:nBlocks, nBlocks), rep(1:nBlocks, each=nBlocks)), 1, sort)))
      groups <- split(1:ncolX,sort(rep(1:nBlocks, ceiling(ncolX/nBlocks))[1:ncolX]))
      cl <- makeCluster(min(cores, detectCores()))
      registerDoParallel(cl)
      cors <- list()
      res <- foreach(i = 1:nrow(blocks)) %dopar% {
		cors[[i]] <- cor(x[, groups[[blocks[i,1]]]], x[, groups[[blocks[i,2]]]], use=use, method=method)
      }
      stopCluster(cl)
      for(i in 1:nrow(blocks)) {
        mat[groups[[blocks[i,1]]], groups[[blocks[i,2]]]] <- res[[i]]
        mat[groups[[blocks[i,2]]], groups[[blocks[i,1]]]] <- t(res[[i]])
      }
      return(mat)
    } else cor(x, y, use=use, method=method)
  }
  MG <- rownames(gpData$geno)[unlist(multiLapply(as.data.frame(is.na(t(gpData$geno))),all, mc.cores=cores))]
  if(check)
    if(impute) {
      if(impute.type %in% c("beagle", "beagleAfterFamily", "beagleNoRand", "beagleAfterFamilyNoRand") & !is.null(gpData$map))
        if(grepl("string mismatches", all.equal(rownames(gpData$map), colnames(gpData$geno))))
          stop("Order of markers in geno and map does not fit!")
      if(impute.type %in% c("beagle") & length(MG) > 0) stop(paste("Of genotype(s) ", MG, " all genotypic values are missing!", sep=" "))
      else if(length(MG) > 0) warning(paste("Of genotype(s) ", MG, " all genotypic values are missing! \nImputation may be erroneus.", sep=" "))
    } else if(length(MG) > 0) warning(paste("Of genotype(s) ", MG, " all genotypic values are missing!", sep=" "))

  orgFormat <- class(gpData)
  # check for class 'gpData'
  if(class(gpData)=="gpData"){
    if(is.null(gpData$geno)) stop("no genotypic data available")
    # family information (population structure) for genotypic data
    # drop unused levels
    if(!is.null(gpData$covar$family)){
      if(is.factor(gpData$covar$family)) {
        popStruc <- as.character(droplevels(gpData$covar$family[gpData$covar$genotyped]))
      } else {
        popStruc <- as.character(gpData$covar$family[gpData$covar$genotyped])
      }
      names(popStruc) <- gpData$covar$id[gpData$covar$genotyped]
    } else {
      popStruc <- NULL
    }
    if(gpData$info$codeGeno & !is.null(label.heter))
      if(is.function(label.heter)){
        warning("assuming heterozygous genotypes coded as 1. Use 'label.heter' to specify if that is not the case")
        label.heter <- "1"
      } else if(is.character(label.heter[1])){
        warning("assuming heterozygous genotypes coded as 1. Use 'label.heter' to specify if that is not the case")
        label.heter <- 1
      }
  } else { # atm other formats are supported too
    if(impute & impute.type %in% c("beagle","beagleAfterFamily","beagleNoRand","beagleAfterFamilyNoRand")) stop("using Beagle is only possible for a gpData object")
    res <- gpData
    popStruc <- NULL
    gpData <- list(geno=res, info=list(map.unit="NA", codeGeno=FALSE))
    rm(res)
  }
  if(impute.type %in% c("beagle","beagleAfterFamily","beagleNoRand","beagleAfterFamilyNoRand")){
    gpData$map <- gpData$map[rownames(gpData$map) %in% colnames(gpData$geno),]
    if(is.null(gpData$map)){
      warning("Beagle imputation makes no sense without map information!")
    } else {
      if(nrow(gpData$map[!is.na(gpData$map$pos),]) != nrow(unique(gpData$map[!is.na(gpData$map$pos),])) &
         !gpData$info$map.unit %in% c("cM", "M"))
        warning("Markers with identical positions will be coded with subsequent basepair numbers!")
    }
  }
  if(is.null(gpData$map)){
    gpData$map <- data.frame(row.names=colnames(gpData$geno),chr=rep(NA, ncol(gpData$geno)),pos=rep(NA, ncol(gpData$geno)))
  } else if(nrow(gpData$map) < ncol(gpData$geno)){
    mnm <- colnames(gpData$geno)[!colnames(gpData$geno) %in% rownames(gpData$map)]
    gpData$map <- rbind(gpData$map, data.frame(row.names=mnm, chr=rep(NA, length(mnm)), pos=rep(NA, length(mnm))))
    gpData$map <- gpData$map[colnames(gpData$geno),]
  }
  if(!is.null(attr(gpData$geno, "identical"))){
    df.ldOld <- attr(gpData$geno, "identical")
    if(!"removed.refer" %in% colnames(df.ldOld)) {
      df.ldOld$removed.refer <- NA
      df.ldOld$removed.alter <- NA
    }
  }  else df.ldOld <- NULL
  #  catch errors
  if(check){
    if(class(gpData$geno)!= "data.frame" & class(gpData$geno) != "matrix") stop("wrong data format")
    if(any(colMeans(is.na(gpData$geno))==1)) warning("markers with only missing values in data")
    if(length(reference.allele)>1 & length(reference.allele)!=ncol(gpData$geno)) stop("'reference allele' should be of length 1 or match the number of markers")
    if(class(reference.allele) != mode(gpData$geno)) stop("'reference allele' should be of class character")
  }
  # number of genotypes
  n <- nrow(gpData$geno)

  # keep names of data object
  cnames <- colnames(gpData$geno)
  rnames <- rownames(gpData$geno)
  popStruc <- popStruc[rnames]
  gpData$geno <- matrix(unlist(gpData$geno),nrow=n)
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
    if((impute.type=="family" | impute.type=="beagleAfterFamily" | impute.type=="beagleAfterFamilyNoRand") & is.null(popStruc)) stop(paste("family information needed, but
    '",substitute(gpData),"$covar$family' is empty",sep=""))
    if((impute.type=="family" | impute.type=="beagleAfterFamily" | impute.type=="beagleAfterFamilyNoRand") & !is.null(popStruc)){
      if(any(is.na(popStruc))) warning("missing values in family information, imputation is likely to be incomplete")
      if(length(popStruc)!=n) stop("population structure must have equal length as obsersvations in genotypic data")
    }
  }

  # use same reference allele for all markers if not specified differently
  if(reference.allele[1]!="minor" & length(reference.allele)==1)     reference.allele <- rep(reference.allele,ncol(gpData$geno))
  knames <- cnames %in% keep.list

  #============================================================
  # step 1  - remove markers with more than nmiss fraction of missing values (optional, argument nmiss>0)
  #============================================================

  which.miss <- multiLapply(as.data.frame(is.na(gpData$geno)),mean, mc.cores=cores)
  if(!is.null(nmiss)){
    if(nmiss<0 | nmiss>1) stop("'nmiss' have to be in [0,1]")
    which.miss <- which.miss<=nmiss | knames
    gpData$geno <- gpData$geno[,which.miss]
    if(!(reference.allele[1]=="minor" | reference.allele[1]=="keep"))  reference.allele <- reference.allele[which.miss]
    if (verbose) cat("   step 1  :", sum(!which.miss),"marker(s) removed with >",nmiss*100,"% missing values \n")
    cnames <- cnames[which.miss]; knames <- knames[which.miss]
    # update map
    if(!is.null(gpData$map)) gpData$map <- gpData$map[rownames(gpData$map) %in% cnames,]
    rm(which.miss)
  } else if (any(which.miss==1)){
    which.miss <- which.miss!=1
    gpData$geno <- gpData$geno[,which.miss]
    if(!(reference.allele[1]=="minor" | reference.allele[1]=="keep"))  reference.allele <- reference.allele[which.miss]
    if (verbose) cat("   step 1  :", sum(!which.miss),"marker(s) removed with only missing values \n")
    cnames <- cnames[which.miss]; knames <- knames[which.miss]
    # update map
    if(!is.null(gpData$map)) gpData$map <- gpData$map[rownames(gpData$map) %in% cnames,]
    rm(which.miss)
  } else {
    if (verbose) cat("   step 1  : No markers removed due to fraction of missing values \n")
  }

  #============================================================
  # step 2  - coding alleles
  #============================================================

  if (verbose) cat("   step 2  : Recoding alleles \n")
  if(ploidy<3) midDose <- 1 else midDose <- .5
  if(gpData$info$codeGeno) {
    if(reference.allele[1]=="minor" | reference.allele[1]!="keep"){
      afCols <- cnames[colMeans(gpData$geno, na.rm=TRUE) > midDose]
      gpData$geno[, cnames%in%afCols] <-  rep(1, nrow(gpData$geno)) %*% t(rep(2, length(afCols))) - gpData$geno[, cnames%in%afCols]
      if(all(c("refer", "alter") %in% colnames(gpData$map))){
        gpData$map[cnames%in%afCols,c("refer", "alter")] <- gpData$map[cnames%in%afCols,c("alter", "refer")]
        if(!is.null(attr(gpData$geno, "identical")))
          attr(gpData$geno, "identical")[attr(gpData$geno, "identical")$kept %in% afCols, c("removed.refer", "removed.alter")] <-
             attr(gpData$geno, "identical")[attr(gpData$geno, "identical")$kept %in% afCols, c("removed.alter", "removed.refer")]
      } else {
        df.alleles <- matrix(rep(0:2, each=ncol(gpData$geno)), ncol=3)
        df.alleles[cnames%in%afCols,c(1,3)] <- df.alleles[cnames%in%afCols,c(3,1)]
        if(!is.null(attr(gpData$geno, "identical")))
          attr(gpData$geno, "identical")[attr(gpData$geno, "identical")$kept %in% afCols, c("removed.refer", "removed.alter")] <-
             attr(gpData$geno, "identical")[attr(gpData$geno, "identical")$kept %in% afCols, c("removed.alter", "removed.refer")]
        gpData$map <- cbind(gpData$map, df.alleles)
      }
    }
  } else {# codeGeno condition of gpData FALSE
    if(ploidy < 3){
      extract <- function(x,y){x[y]}
      colnames(gpData$geno) <- cnames
      if(is.numeric(gpData$geno)){
        alleles <- multiLapply(as.data.frame(gpData$geno),unique,mc.cores=cores)
      } else {
        alleles <- multiLapply(as.data.frame(gpData$geno),levels,mc.cores=cores)
      }
      names(alleles) <- cnames
      gpData$geno <- multiLapply(as.data.frame(gpData$geno), as.numeric, mc.cores=cores)
      df.allele <- data.frame(id=rownames(gpData$map), refer=NA, heter=NA, alter=NA, stringsAsFactors=FALSE)
      if(is.null(label.heter)) {
        df.allele$refer <- unlist(multiLapply(alleles, extract, 1, mc.cores=cores))
        df.allele$alter <- unlist(multiLapply(alleles, extract, 2, mc.cores=cores))
        df.allele$heter <- unlist(multiLapply(alleles, extract, 2, mc.cores=cores))
      } else {
        whereHetPos <- function(x, y=NULL){if(is.function(y)) z <- c((1:3)[y(x)], 3)
                                           else if(y=="alleleCoding") z <- c((1:3)[substr(x, 1, 1) != substr(x, nchar(x), nchar(x))],3)
                                           else z <- c((1:3)[x==y], 3)
                                           return(z[1])}
        hetPos <- numeric()
        if(length(label.heter) == 1){
          label.heter <- as.list(label.heter)
          hetPos <- c(unlist(multiLapply(alleles, whereHetPos, label.heter, mc.cores=cores)))
        } else {
          if(length(label.heter) != length(gpData$geno))
            label.heter <- rep(as.list(label.heter), ceiling(length(gpData$geno)/length(label.heter)))[1:length(gpData$geno)]
          for(i in unique(as.character(label.heter))){
            j <- label.heter[[match(i, as.character(label.heter))]]
            namWk <- names(alleles)[unlist(multiLapply(label.heter, identical, j, mc.cores=cores))]
            hetPos <- c(unlist(multiLapply(alleles[namWk], whereHetPos, j, mc.cores=cores)), hetPos)
          }
        }
        hetPos <- hetPos[names(alleles)]
        df.allele$heter <- unlist(multiLapply(alleles, extract, 2, mc.cores=cores))
        df.allele$refer <- unlist(multiLapply(alleles, extract, 1, mc.cores=cores))
        df.allele$alter <- unlist(multiLapply(alleles, extract, 3, mc.cores=cores))
        gpData$geno <- matrix(unlist(gpData$geno), ncol=length(gpData$geno), dimnames=list(1:n, names(gpData$geno)))
        gpData$geno[, hetPos==1] <- unlist(multiLapply(gpData$geno[, hetPos==1], function(x){2*(x%%2)+(x%/%2)}, mc.cores=cores))
        gpData$geno[, hetPos==3] <- unlist(multiLapply(gpData$geno[, hetPos==3], function(x){2*((1+x)%%2)+((1+x)%/%2)}, mc.cores=cores))
        df.allele[hetPos==1,2:3] <- df.allele[hetPos==1,3:2]
        df.allele[hetPos==3,4:3] <- df.allele[hetPos==3,3:4]
        alterMiss <- !is.na(df.allele$heter) & is.na(df.allele$alter)
        if(!all(alterMiss)){
          df.allele$alter[alterMiss] <- unlist(multiLapply(as.data.frame(t(df.allele[alterMiss,]), stringsAsFactors=FALSE),
                          function(x){alt <- unlist(strsplit(x[3],""))[!unlist(strsplit(x[3],"")) %in% unlist(strsplit(x[2],""))]
                                      mid <- substr(x[3], 2, nchar(x[3])-1)
                                      return(paste(alt,mid,alt,sep=""))}))
        }
      }
      gpData$geno <- gpData$geno-1
      if(reference.allele[1]=="minor"){
        afCols <- cnames[colMeans(gpData$geno, na.rm=TRUE)>midDose]
        gpData$geno[, cnames%in%afCols] <-  rep(1, nrow(gpData$geno)) %*% t(rep(2, length(afCols))) - gpData$geno[, cnames%in%afCols]
        df.allele[cnames%in%afCols,c(2,4)] <- df.allele[cnames%in%afCols,c(4,2)]
      }
      if(all.equal(colnames(gpData$geno), rownames(gpData$map)))
        gpData$map <- cbind(gpData$map, df.allele[, c("refer", "heter", "alter")])
      else gpData$alleles <- df.allele
    } else { # ploidy
      colnames(gpData$geno) <- cnames
      if(is.numeric(gpData$geno))
        alleles <- multiLapply(as.data.frame(gpData$geno),unique,mc.cores=cores)
      else {
        alleles <- multiLapply(as.data.frame(gpData$geno),levels,mc.cores=cores)
        if(!all(unlist(multiLapply(alleles, strsplit,"", cores=cores)%in% c("A","B")))) stop("Wrong coding for multiploid species. Only A and B is allowed!")
      }
      names(alleles) <- cnames
      gpData$geno <- multiLapply(as.data.frame(gpData$geno), as.numeric, mc.cores=cores)
      gpData$geno <- matrix((unlist(gpData$geno)-1)/ploidy, ncol=length(gpData$geno), dimnames=list(1:n, names(gpData$geno)))
      df.allele <- data.frame(id=rownames(gpData$map), refer=NA,alter=NA, stringsAsFactors=FALSE)
      df.allele$refer <- unlist(multiLapply(alleles, function(x,y){substr(x[y],1,1)}, 1, mc.cores=cores))
      df.allele$alter <- unlist(multiLapply(alleles, function(x,y){substr(x[y],nchar(x),nchar(x))}, ploidy+1, mc.cores=cores))
      if(reference.allele[1]=="minor"){
        afCols <- cnames[colMeans(gpData$geno, na.rm=TRUE)>midDose]
        gpData$geno[, cnames%in%afCols] <-  rep(1, nrow(gpData$geno)) %*% t(rep(1, length(afCols))) - gpData$geno[, cnames%in%afCols]
        df.allele[cnames%in%afCols,c(1,2)] <- df.allele[cnames%in%afCols,c(2,1)]
        if(!is.null(attr(gpData$geno, "identical")))
          attr(gpData$geno, "identical")[attr(gpData$geno, "identical")$kept %in% afCols, c("removed.refer", "removed.alter")] <-
             attr(gpData$geno, "identical")[attr(gpData$geno, "identical")$kept %in% afCols, c("removed.alter", "removed.refer")]
      }
      if(all.equal(colnames(gpData$geno), rownames(gpData$map)))
        gpData$map <- cbind(gpData$map, df.allele[, c("refer", "alter")])
      else gpData$alleles <- df.allele
    }
  }

  #============================================================
  # step 3  - Discarding markers for which the tester is not homozygous or values missing (optional, argument tester = "xxx")
  #============================================================

  if(!is.null(tester)){
    which.miss <- gpData$geno[rnames==tester,]!=label.heter&!is.na(gpData$geno[rnames==tester,])| knames
    gpData$geno <- gpData$geno[,which.miss]
    cnames <- cnames[which.miss]; knames <- knames[which.miss]
    if(sum(!which.miss) > 0){
      if (verbose) cat("   step 3 :",sum(!which.miss),"marker(s) discarded because heterozygousity at tester locus or \n          missing values of the tester\n")
    } else {
      if (verbose) cat("   step 3 : No marker(s) discarded because heterozygousity at tester locus or \n          missing values of the tester\n")
    }
    # update map
    if(!is.null(gpData$map)) gpData$map <- gpData$map[rownames(gpData$map) %in% cnames,]
  }

  #============================================================
  # step 4 - remove markers with minor allele frequency < maf  (optional, argument maf>0)
  #============================================================

  if(!is.null(maf)){
    if(maf<0 | maf>1) stop("'maf' must be in [0,1]")
    if(is.null(tester)){
        which.maf <- unlist(multiLapply(as.data.frame(gpData$geno),mean, na.rm=TRUE, mc.cores=cores))>=2*maf | knames
    } else {
        which.maf <- unlist(multiLapply(as.data.frame(gpData$geno),mean, na.rm=TRUE, mc.cores=cores))
        which.maf <- which.maf>=maf & which.maf<=1-maf | knames
    }
    if (verbose) cat("   step 4  :",sum(!which.maf),"marker(s) removed with maf <",maf,"\n")
    gpData$geno <- gpData$geno[,which.maf]
    cnames <- cnames[which.maf]; knames <- knames[which.maf]
    # update map
    if(!is.null(gpData$map)) gpData$map <- gpData$map[rownames(gpData$map) %in% cnames,]
     # update report list

  } else {
    if (verbose) cat("   step 4  : No markers discarded due to minor allele frequency \n")
  }
  #============================================================
  # step 5  - Discarding markers for which the tester has the minor allele
  #============================================================

  if(!is.null(tester)){
    which.miss <- gpData$geno[rnames==tester,] != 2 | knames
    gpData$geno <- gpData$geno[,which.miss]
    cnames <- cnames[which.miss]; knames <- knames[which.miss]
    if(sum(!which.miss) > 0){
      if (verbose) cat("   step 5  :",sum(!which.miss),"marker(s) discarded for which the tester has the minor allele\n")
    } else{
      if (verbose) cat("   step 5  : No marker(s) discarded for which the tester has the minor allele\n")
    }
    # update map
    if(!is.null(gpData$map)) gpData$map <- gpData$map[rownames(gpData$map) %in% cnames,]
  }

  #============================================================
  # step 6  - Discarding homozygout values of the minor allele and markers with more than nmiss values
  #============================================================

  if(!is.null(tester)){
    gpData$geno[gpData$geno == 2] <- NA
    gpData$geno <- matrix(as.numeric(gpData$geno), nrow = n)
    if(!is.null(nmiss)){
      which.miss <- multiLapply(as.data.frame(is.na(gpData$geno)),mean,na.rm=TRUE, mc.cores=cores) <= nmiss |knames
      gpData$geno <- gpData$geno[,which.miss]
      cnames <- cnames[which.miss]; knames <- knames[which.miss]
      if (verbose) cat("   step 6  :",sum(!which.miss),"marker(s) discarded with >",nmiss*100,"% false genotyping values \n")
      # update map
      if(!is.null(gpData$map)) gpData$map <- gpData$map[rownames(gpData$map) %in% cnames,]
    } else{
      if (verbose) cat("   step 6  : No markers discarded due to fraction of missing values \n")
    }
  }

  #============================================================
  # step 7  - imputing missing genotypes  (optional, argument impute=TRUE)
  #============================================================

  # initialize counter
    cnt1 <- rep(0,ncol(gpData$geno))   # for nr. of imputations with family structure
    cnt2 <- rep(0,ncol(gpData$geno))    # for nr. of beagle imputations
    cnt3 <- rep(0,ncol(gpData$geno))    # for nr. of random imputations
    names(cnt1) <- names(cnt2) <- names(cnt3) <- cnames

  # start of imputing
  if(impute){
    set.seed(SEED[1])
    # number of markers
    M <- ncol(gpData$geno)
    if(M==0) stop(" no markers remained after step 1 (to many missing values)")
    if (verbose) cat("   step 7  : Imputing of missing values \n")

    # number of missing values
    nmv <- sum(is.na(gpData$geno))

    ###########################################################################
    # if impute.type="fix", replace missing values according to specified value
    ###########################################################################
    if(impute.type=="fix"){
      gpData$geno[is.na(gpData$geno)] <- replace.value
      if (verbose) cat("   step 7a : Replace missing values by",replace.value," \n")
    }
    #########################################################
    # impute missing values according to population structure
    #########################################################
    if(impute.type %in% c("family", "beagleAfterFamily", "beagleAfterFamilyNoRand")){
      if (verbose) cat("   step 7b : Imputing of missing values by family information \n")
      # initialize counter (- number of heterozygous values)
      # loop over all markers
      probList <- list(c(1), c(.5,.5), c(.25,.5,.25))
      vec.cols <- (1:M)[is.na(colSums(gpData$geno, na.rm = FALSE))]
      nFam <- table(popStruc)
      vec.big <- popStruc %in% names(nFam)[nFam > minFam]
      ptm <- Sys.time()
      for (j in vec.cols){
        if(sum(!is.na(gpData$geno[,j]))>0){
          poptab <- table(popStruc[vec.big],gpData$geno[vec.big,j])
          rS <- rowSums(poptab)
          # compute otherstatistics
          major.allele <- unlist(attr(poptab,"dimnames")[[2]][apply(poptab,1,which.max)])
          # look if SNP is segregating  for this population
          polymorph <- apply(poptab,1,length) > 1 & (apply(poptab,1,min) != 0)
          polymorph2 <- rS > minFam
          polymorph[!polymorph2] <- TRUE
          # count missing values
          nmissfam <- tapply(is.na(gpData$geno[vec.big,j]),popStruc[vec.big],sum)
          # must be a named list
          names(major.allele) <- names(polymorph)
          # loop over all families
          for (i in rownames(poptab)[nmissfam > 0]){
            # impute values for impute.type="family" : all missing genotypes
            allTab <- table(gpData$geno[popStruc[vec.big] %in% i, j])
            if(length(allTab) == 1){
              gpData$geno[is.na(gpData$geno[,j]) & popStruc %in% i ,j] <- as.numeric(names(allTab))
              cnt1[j] <- cnt1[j] + nmissfam[as.character(i)]
            } else if(impute.type %in% c("family", "beagleAfterFamily")){
              if(length(allTab) == 0 & noHet) {
                allTab <- table(c(0,2))
              } else if(all(names(allTab) == c(0, 2)) & !noHet){
                allTab <- table(c(0,1,1,2))
              }
              if (impute.type=="family"|is.na(gpData$map$pos[j])){
                gpData$geno[is.na(gpData$geno[,j]) & popStruc %in% i ,j] <- ifelse(length(allTab)>1,
                                                               sample(as.numeric(names(allTab)),size=nmissfam[as.character(i)],prob=probList[[length(allTab)]],replace=TRUE),
                                                               as.numeric(names(allTab)))
                cnt3[j] <- cnt3[j] + nmissfam[as.character(i)]
              }
            }
          }
          if(j==ceiling(length(vec.cols)/50))
            if(verbose)  cat("         approximative run time for imputation by family information ",
                             paste(round(as.numeric(difftime(Sys.time(), ptm)*50)), digits=1, " ",
                             units(difftime(Sys.time(), ptm))," ... \n",sep=""))
        }   # end of if(sum(!is.na(gpData$geno[,j]))>0)
      } # end of marker loop
    }

    ###########################
    # run beagle for imputation
    ###########################
    if(impute.type %in% c("beagle","beagleAfterFamily","beagleNoRand","beagleAfterFamilyNoRand")){
      if (verbose) cat("   step 7c : Imputing of missing values by Beagle \n")
        #if(any(grep(" ",path.package()[grep("synbreed", path.package())])))
          #warning("The package is installed in folder ",path.package()[grep("synbreed", path.package())]," which contains a space. Torun beagle properly, please install the package to a differnt folder without spaces.")
        # use Beagle and impute NA for polymorphic families
        chr <-  unique(gpData$map$chr)
        chr <- chr[!is.na(chr)]
        if(!is.null(tester))
          gpData$geno <- gpData$geno*2
        rownames(gpData$geno) <- rnames
        colnames(gpData$geno) <- cnames
        cnt2 <- unlist(multiLapply(as.data.frame(is.na(gpData$geno)),sum, mc.cores=cores))
        pre <- paste(as.numeric(as.Date(Sys.time())), round(as.numeric(Sys.time())%%(24*3600)), sep="")
        if(!"beagle" %in% list.files()){
          dir.create("beagle")
        }
        markerTEMPbeagle <- discard.markers(gpData,whichNot=rownames(gpData$map[!is.na(gpData$map$pos),]))

        # write input files for beagle
        # create new directory "beagle" for beagle input and output files
        if(gpData$info$map.unit %in% c("Mb", "kb", "bp")){
          mapfile <- NULL
          while(any(duplicated(markerTEMPbeagle$map))) {markerTEMPbeagle$map$pos[duplicated(markerTEMPbeagle$map)] <- markerTEMPbeagle$map$pos[duplicated(markerTEMPbeagle$map)]+1}
        } else {
          mapfile <- data.frame(markerTEMPbeagle$map$chr, rownames(markerTEMPbeagle$map), markerTEMPbeagle$map$pos, markerTEMPbeagle$map$pos)
          if(!is.integer(mapfile[,1])) mapfile[,1] <- as.integer(as.factor(mapfile[,1]))
          if(gpData$info$map.unit == "M") { mapfile[, 4] <- 1000000 * (mapfile[, 3] <- mapfile[, 3] * 100 )}
          if(gpData$info$map.unit == "cM") mapfile[, 4] <- 1000000 * mapfile[, 3]
          while(any(duplicated(markerTEMPbeagle$map))) {mapfile[duplicated(markerTEMPbeagle$map), 4] <- mapfile[duplicated(markerTEMPbeagle$map), 4]+1}
          markerTEMPbeagle$map$pos <- mapfile[, 4]
          mapfile[,4] <- formatC(mapfile[,4], format="f", digits=0)
          mapfile[,1] <- paste("chr",mapfile[,1],sep="")
          write.table(mapfile, file=paste("beagle/run", pre, ".map", sep=""), col.names=FALSE, row.names=FALSE, quote=FALSE, na=".", sep="\t")
          mapfile <- paste(" map=beagle/run", pre, ".map ", sep="")
          markerTEMPbeagle$info$map.unit <- "bp"
        }
        markerTEMPbeagle$map$chr <- as.numeric(as.factor(markerTEMPbeagle$map$chr))
        write.vcf(markerTEMPbeagle,paste(file.path(getwd(),"beagle"),"/run",pre,"input.vcf", sep=""))
        if(noHet){
        output <- system(paste("java -Xmx3000m -jar ",
                         shQuote(paste(sort(path.package()[grep("synbreed", path.package())])[1], "/java/beagle.21Jan17.6cc.jar", sep="")),
                         # caution with more than one pacakge with names synbreed*, assume synbreed to be the first one
                         " gtgl=beagle/run", pre, "input.vcf out=beagle/run", pre, "out gprobs=true nthreads=", cores, mapfile, sep=""),
                         intern=!showBeagleOutput)
        } else {
        output <- system(paste("java -Xmx3000m -jar ",
                         shQuote(paste(sort(path.package()[grep("synbreed", path.package())])[1], "/java/beagle.21Jan17.6cc.jar", sep="")),
                         # caution with more than one pacakge with names synbreed*, assume synbreed to be the first one
                         " gtgl=beagle/run", pre, "input.vcf out=beagle/run", pre, "out gprobs=true nthreads=", cores, mapfile, sep=""),
                         intern=!showBeagleOutput)
        }
        # read data from beagle
        gz <- gzfile(paste("beagle/run",pre,"out.vcf.gz",sep=""))
        resTEMP <- read.vcf2matrix(file=gz, FORMAT="DS", IDinRow=TRUE, cores=cores)
        mode(resTEMP) <- "numeric"

        # convert dose to genotypes
        if(noHet){
          resTEMP[resTEMP<1] <- 0
          resTEMP[resTEMP>=1] <- 2
        } else {
          resTEMP <- round(resTEMP,0) # 0, 1, and 2
        }
        gpData$geno[,colnames(resTEMP)] <- resTEMP
      }

      #########################################################################
      # impute missing values with no population structure or missing positions
      #########################################################################
      if(impute.type %in% c("random", "beagle", "beagleAfterFamily", "family")){
        if (verbose) cat("   step 7d : Random imputing of missing values \n")
        # initialize counter (- number of heterozygous values)
        mImp <- unlist(multiLapply(as.data.frame(is.na(gpData$geno)),sum, mc.cores=cores))
        p <- unlist(multiLapply(as.data.frame(gpData$geno),mean,na.rm=TRUE, mc.cores=cores))/2
        ptm <- proc.time()[3]
        for (j in (1:M)[mImp>0]){
          cnt3[j] <-  cnt3[j] + mImp[j]
          # estimation of running time after the first iteration
          if(noHet){        # assuming only 2 homozygous genotypes
            gpData$geno[is.na(gpData$geno[,j]),j] <- sample(c(0,2),size=mImp[j],prob=c(1-p[j],p[j]),replace=TRUE)
          } else {                            # assuming 3 genotypes
            gpData$geno[is.na(gpData$geno[,j]),j] <- sample(c(0,1,2),size=mImp[j],prob=c((1-p[j])^2,2*p[j]*(1-p[j]),p[j]^2),replace=TRUE)
          }
          if(j==ceiling(M/100)) if(verbose) cat("         approximate run time for random imputation ",(proc.time()[3] - ptm)*99," seconds \n",sep=" ")
        }
        # update counter for Beagle, remove those counts which where imputed ranomly
        if(impute.type == "beagle") cnt2 <- cnt2-cnt3
      }
      if(!is.null(tester) & impute.type %in% c("random","beagle", "beagleAfterFamily")) gpData$geno <- gpData$geno/2

    #============================================================
    # step 8 - recoding
    #============================================================

    # recode again if allele frequeny changed to to imputing
    p <- unlist(multiLapply(as.data.frame(gpData$geno),mean,na.rm=TRUE, mc.cores=cores))
    if(any(p>1)&reference.allele[1]=="minor"){
      if (verbose) cat("   step 8  : Recode alleles due to imputation \n")
      gpData$geno[,which(p>1)] <- 2 - gpData$geno[,which(p>1)]
      p[which(p>1)] <- 2-p[which(p>1)]
    } else{
      if (verbose) cat("   step 8  : No recoding of alleles necessary after imputation \n")
    }
  }

  #============================================================
  # step 9 - remove markers with minor allele frequency < maf  (optional, argument maf>0)
  #============================================================

  if(!is.null(maf) & impute){
    if(maf<0 | maf>1) stop("'maf' must be in [0,1]")
    if(is.null(tester)) which.maf <- p>=2*maf | knames else
      which.maf <- p>=maf & p<=1-maf | knames
    if (verbose) cat("   step 9  :",sum(!which.maf),"marker(s) removed with maf <",maf,"\n")
    gpData$geno <- gpData$geno[,which.maf]
    cnames <- cnames[which.maf]; knames <- knames[which.maf]; p <- p[which.maf]
    # update map
    if(!is.null(gpData$map)) gpData$map <- gpData$map[rownames(gpData$map) %in% cnames,]
     # update report list


  } else {
    if (verbose & impute) cat("   step 9  : No markers discarded due to minor allele frequency \n")
  }

  #============================================================
  # step 10 - discard duplicated markers   (optional, argument keep.identical=FALSE)
  #============================================================

  if(!keep.identical){
    set.seed(SEED[2])
    colnames(gpData$geno) <- cnames
    cnms <- sample(1:ncol(gpData$geno))
    gpData$geno <- gpData$geno[,cnms]; cnames <- cnames[cnms]; knames[cnms]
    which.duplicated <- duplicated(gpData$geno,MARGIN=2)
    rev.which.duplicated <- duplicated(gpData$geno,MARGIN=2, fromLast=TRUE)
    rev.which.duplicated[which.duplicated] <- FALSE
    if(impute){
       if(sum(which.duplicated) >0){
        mat.ld <- multiCor(gpData$geno[, which.duplicated], gpData$geno[, rev.which.duplicated], use="pairwise.complete.obs", cores=cores)
        df.ld <- data.frame(kept=rep(colnames(mat.ld),each=nrow(mat.ld)),
                            removed=rep(rownames(mat.ld),ncol(mat.ld)),
                            ld=as.numeric(mat.ld),
                            stringsAsFactors=FALSE)
        df.ld <- df.ld[abs(df.ld$ld)>1-1e-14,]
        df.ld$ld <- NULL
        rm(mat.ld)
      } else df.ld <- data.frame(kept=as.character(), removed=as.character())
    } else {# end of imputed step
      if(!all(!is.na(gpData$geno))){
        if(sum(which.duplicated) >0){
          gpData$geno[is.na(gpData$geno)] <- 3
          mat.ld <- multiCor(gpData$geno[, which.duplicated], gpData$geno[, rev.which.duplicated], use="pairwise.complete.obs", cores=cores)
          print(mat.ld)
          df.ld <- data.frame(kept=rep(cnames[rev.which.duplicated], each=nrow(mat.ld)),
                              removed=rep(cnames[which.duplicated], ncol(mat.ld)),
                              ld=as.numeric(mat.ld),
                              stringsAsFactors=FALSE)
          df.ld <- df.ld[abs(df.ld$ld)>1-1e-14,]
          df.ld$ld <- NULL
          gpData$geno[gpData$geno==3] <- NA
          rm(mat.ld)
        } else df.ld <- data.frame(kept=as.character(), removed=as.character())
        which.miss <- unlist(multiLapply(as.data.frame(gpData$geno),function(x){sum(is.na(x))}, mc.cores=cores))>0
	    which.miss <- (1:length(which.miss))[which.miss] 	
        if(length(which.miss[which.miss]) == ncol(gpData$geno))
          which.miss <- which.miss[1:(length(which.miss)-1)] 	
        if(is.null(keep.list)){
          for(i in which.miss){ 	
            if(which.duplicated[i]) next 	
            for(j in ((i+1):ncol(gpData$geno))[!which.duplicated[(i+1):ncol(gpData$geno)]]){ 	
              if(all(gpData$geno[, i] == gpData$geno[, j], na.rm = TRUE)){ 	
                if(sum(is.na(gpData$geno[, i])) >= sum(is.na(gpData$geno[, j]))){ 	
                  which.duplicated[i] <- TRUE 	
                  df.ld <- rbind(df.ld, data.frame(kept=cnames[j], removed=cnames[i]))
                  break 	
                } else { 	
                  which.duplicated[j] <- TRUE 	
                  df.ld <- rbind(df.ld, data.frame(kept=cnames[i], removed=cnames[j]))
                }
              }
            }
          }
        } else {
          for(i in which.miss){ 	
            if(which.duplicated[i]) next 	
            for(j in ((i+1):ncol(gpData$geno))[!which.duplicated[(i+1):ncol(gpData$geno)]]){ 	
              if(all(gpData$geno[, i] == gpData$geno[, j], na.rm = TRUE)){
                if(knames[i]){# knames is logical vector for keep.list. Faster than testing if cnames[i] in keep.list!
                  if(knames[j]) next else which.duplicated[j] <- TRUE
                } else {
                  if(knames[j]){ which.duplicated[i] <- TRUE
                  } else {
                    if(sum(is.na(gpData$geno[, i])) >= sum(is.na(gpData$geno[, j]))){ 	
                      which.duplicated[i] <- TRUE 	
                      df.ld <- rbind(df.ld, data.frame(kept=cnames[j], removed=cnames[i]))
                      break 	
                    } else { 	
                      which.duplicated[j] <- TRUE 	
                      df.ld <- rbind(df.ld, data.frame(kept=cnames[i], removed=cnames[j]))
                    } # end choice of not keep.list elements
                  } # end of neither i or j in keep.list
                } # end of i not in the keep list
              } # end of to equal i and j proove
            } # end of the loop from i+1 to j
          } # end of loop through which.miss
        } # end of else step of is.null(keep.list) proove
        if(is.na(df.ld[1,1])) df.ld <- df.ld[-1,]
      } else {# end of missing value step
        if(sum(which.duplicated) >0){
          mat.ld <- multiCor(gpData$geno[, which.duplicated], gpData$geno[, rev.which.duplicated], use="pairwise.complete.obs", cores=cores)
          rownames(mat.ld) <- cnames[which.duplicated]
          colnames(mat.ld) <- cnames[rev.which.duplicated]
          df.ld <- data.frame(kept=rep(colnames(mat.ld), nrow(mat.ld)),
                              removed=rep(rownames(mat.ld), each=ncol(mat.ld)),
                              ld=as.numeric(mat.ld),
                              stringsAsFactors=FALSE)
          df.ld <- df.ld[abs(df.ld$ld)>1-1e-14,]
          df.ld$ld <- NULL
          rm(mat.ld)
        } else df.ld <- data.frame(kept=as.character(), removed=as.character())
      }
    } # end of not imputed step
    gpData$geno <- gpData$geno[, !which.duplicated]
    df.ld$removed.refer <- gpData$map[df.ld$removed, "refer"]
    df.ld$removed.alter <- gpData$map[df.ld$removed, "alter"]
    df.ld <- rbind(df.ldOld[, colnames(df.ld)], df.ld)
    df.ld$sort <- match(df.ld$kept, rownames(gpData$map))
    df.ld <- orderBy(~sort+removed, df.ld)
    df.ld$sort <- NULL
    attr(gpData$geno, "identical") <- df.ld
    cnames <- cnames[!which.duplicated]
    if (verbose) cat("   step 10 :",sum(which.duplicated),"duplicated marker(s) removed \n")
    # update map
    if(!is.null(gpData$map)) gpData$map <- gpData$map[rownames(gpData$map) %in% cnames,]

  } else{
    if (verbose) cat("   step 10 : No duplicated markers removed \n")
  }

  #============================================================
  # step 10a - discard markers for which only the tester is different
  #============================================================

  if(!is.null(tester)){
    which.fixed <- multiLapply(as.data.frame(gpData$geno), sum, cores=cores) == nrow(gpData$geno)-1 | knames
    gpData$geno <- gpData$geno[,!which.fixed]
    cnames <- cnames[!which.fixed]; knames <- knames[!which.fixed]
    if(!is.null(gpData$map)) gpData$map <- gpData$map[rownames(gpData$map) %in% cnames,]
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

  rownames(gpData$geno) <- rnames
  colnames(gpData$geno) <- cnames
  if(orgFormat == "matrix"){
    gpData$geno <- matrix(gpData$geno,nrow=n)
  }
  if(orgFormat == "data.frame"){
    gpData$geno  <- as.data.frame(gpData$geno)
  }

  if (verbose) cat("   End     :",ncol(gpData$geno),"marker(s) remain after the check\n")

  #============================================================
  # print summary of imputation
  #============================================================

  if(impute){
    cat("\n")
    cat("     Summary of imputation \n")
    cat(paste("    total number of missing values                :",nmv,"\n"))
    if(impute.type %in% c("family","beagleAfterFamily","familyNoRand","beagleAfterFamilyNoRand")) cat(paste("    number of imputations by family structure     :",sum(cnt1),"\n"))
    if(impute.type %in% c("beagle","beagleAfterFamily","beagleNoRand","beagleAfterFamilyNoRand")) cat(paste("    number of Beagle imputations                  :",sum(cnt2),"\n"))
    if(impute.type %in% c("beagle","random","family","beagleAfterFamily"))                        cat(paste("    number of random imputations                  :",sum(cnt3),"\n"))
  }

  if(!is.null(gpData$map)){
    gpData$map$sor <- substr(gpData$map$chr, nchar(as.character(gpData$map$chr)), nchar(as.character(gpData$map$chr)))
    if(any(unique(gpData$map$sor)[!is.na(unique(gpData$map$sor))] %in% 0:9)) gpData$map$sor <- 1
    # first order by rownames in alphabetical order (important for SNPs with the same position)
    gpData$map <- gpData$map[order(as.character(rownames(gpData$map))),]
    gpData$map <- orderBy(~sor+chr+pos,data=as.data.frame(gpData$map))
    gpData$map$sor <- NULL
    # sortcolumns in geno, too
    if(!is.null(attr(gpData$geno, "identical"))) attrG <- attr(gpData$geno, "identical") else attrG <- NULL
    gpData$geno <- gpData$geno[,rownames(gpData$map)]
    if(!is.null(attrG)) attr(gpData$geno, "identical") <- attrG
  }
  # overwrite original genotypic data
  if(orgFormat == "gpData") {
    gpData$info$codeGeno <- TRUE
    gpData$info$version <- paste("gpData object was coded by synbreed version", sessionInfo()$otherPkgs$synbreed$Version)
    gpData$info$Call <- infoCall
  }
  if(print.report){
    if (verbose) cat("  Writing report to file 'SNPreport.txt' \n")
    report.list <- data.frame(SNPname=cnames,reference=gpData$map$refer, alternative=gpData$map$alter,
                              MAF=round(colMeans(gpData$geno,na.rm=TRUE)/2,3),
                              impute.fam=cnt1[cnames], impute.beagle=cnt2[cnames],
                              impute.ran=cnt3[cnames])
    write.table(report.list,file="SNPreport.txt",quote=FALSE,row.names=FALSE)
   }

  # return a gpData object (or a matrix)
  if(reference.allele[1]=="keep"){
    gpData$geno[, cnames%in%afCols] <-  rep(1, nrow(gpData$geno)) %*% t(rep(2, sum(cnames%in%afCols))) - gpData$geno[, cnames%in%afCols]
  }
  return(gpData)
}
