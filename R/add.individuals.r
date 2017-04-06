# adding new individuals to gpData object

add.individuals <- function(gpData,pheno=NULL,geno=NULL,pedigree=NULL,covar=NULL,repl=NULL){
      # check if individuals are allready in data
      if (any(rownames(pheno) %in% gpData$covar$id) | any(rownames(geno) %in% gpData$covar$id) |
          any(pheno$ID %in% gpData$covar$id)){
        stop("some of the individuals of are already in ", substitute(gpData))
      }
      if(is.null(repl)) repl <- "repl"
      else if(!all(unique(pheno[, repl]) %in% dimnames(gpData$pheno)[[3]])) stop("Your values for replication is not in the dimnames of ", substitute(gpData$pheno))
      colnames(pheno)[colnames(pheno) == repl] <- "repl"
      # merge phenotypic data
      if(dim(gpData$pheno)[3] == 1) {
        if(!"ID" %in% colnames(pheno)) pheno$ID <- rownames(pheno)
        repl <- NULL
      } else {
        if(!"ID" %in% colnames(pheno))  stop("In", substitute(pheno), "the columns 'ID' and/or 'repl' are/is missing!")
      }
      if(!is.null(pheno))
        if(!all(!colnames(pheno)[colnames(pheno)!="ID"]%in%dimnames(gpData$pheno)[[2]] |
                !colnames(pheno)[colnames(pheno)!="ID"]%in%dimnames(gpData$phenoCovars)[[2]]))
          stop("different phenotypes (colnames) in '", substitute(gpData$pheno), "' or '",  substitute(gpData$phenoCovar),"' and '", substitute(pheno), "'")
      df.pheno <- gpData2data.frame(gpData, onlyPheno=TRUE, trait = dimnames(gpData$pheno)[[2]])
      if(!all(colnames(df.pheno) %in% colnames(pheno))) warning("Not all traits and covariates are available in the new data!")
      pheno[, colnames(df.pheno)[!colnames(df.pheno) %in% colnames(pheno)]] <- NA
      pheno <- pheno[, colnames(df.pheno)]
      pheno <- rbind(df.pheno, pheno)
      rm(df.pheno)
      if(dim(gpData$pheno)[3] == 1) {
        rownames(pheno) <- pheno$ID
        pheno$ID <- NULL
      }
      # merge genotypic data
      if(!is.null(geno)) if(any(colnames(geno)!=colnames(gpData$geno))) stop("different markers (colnames) in 'gpData$geno' and 'geno'")
      geno <- rbind(gpData$geno,geno)

      # merge pedigree
      pedigree <- rbind(gpData$pedigree,pedigree)
      # reorder if necessary
      if(!is.null(pedigree)) pedigree <- orderBy(~gener,pedigree)

      # merge covar  (not with first columns)

      if(any(colnames(covar) %in% c("genotyped","phenotyped","id"))) stop("not specify columns 'genotyped','phenotyped' and 'id' in 'covar' ")
      if(!is.null(covar)){
        cc <- data.frame(gpData$covar[,!colnames(gpData$covar) %in% c("genotyped","phenotyped","id")])
        rownames(cc) <- rownames(gpData$covar)
        colnames(cc) <- colnames(gpData$covar)[!colnames(gpData$covar) %in% c("genotyped","phenotyped","id")]
        covarUpdate <- rbind(cc,covar)
      }
      else covarUpdate <- gpData$covar

      # need id in rownames for create.gpData
      rownames(covarUpdate) <- c(as.character(gpData$covar$id),rownames(covar))

      # create new gpData object
      ret <- create.gpData(pheno=pheno,geno=geno,map=gpData$map,pedigree=pedigree,covar=covarUpdate,map.unit=gpData$info$map.unit,modCovar=dimnames(gpData$phenoCovars)[[2]],repeated=repl)
      return(ret)
}
