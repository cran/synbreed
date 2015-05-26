# discard markers from class 'gpData'

discard.markers <- function(gpData,which){
  if(class(gpData) != "gpData")
    stop(substitute(gpData), " is not a gpData-object!")
  # update geno
  gpData$geno <- gpData$geno[,!colnames(gpData$geno) %in% which]
  # update map
  gpData$map <- gpData$map[!rownames(gpData$map) %in% which,]
  return(gpData)
}

# discard individuals from class 'gpData'

discard.individuals <- function(gpData,which,keepPedigree=FALSE){
  if(class(gpData) != "gpData")
    stop(substitute(gpData), " is not a gpData-object!")
  if(!all(which %in% gpData$covar$id)) stop("Some individuals are not in the ", substitute(gpData), "-object!")
  if(!is.null(gpData$pheno)){
  # updata pheno
    phenoNames <- dimnames(gpData$pheno)
    phenoNames[[1]] <- phenoNames[[1]][!phenoNames[[1]] %in% which]
    phenoDim <- dim(gpData$pheno)
    phenoDim[1] <- length(phenoNames[[1]])
    gpData$pheno <- array(gpData$pheno[!rownames(gpData$pheno) %in% which, , ], dim = phenoDim)
    dimnames(gpData$pheno) <- phenoNames
  }
  # update geno
  gpData$geno <- subset(gpData$geno,!rownames(gpData$geno) %in% which)
  if(!is.null(gpData$phenoCovars)){
    phenoCovarsNames <- dimnames(gpData$phenoCovars)
    phenoCovarsNames[[1]] <- phenoCovarsNames[[1]][!phenoCovarsNames[[1]] %in% which]
    phenoCovarsDim <- dim(gpData$phenoCovars)
    phenoCovarsDim[1] <- length(phenoCovarsNames[[1]])
    gpData$phenoCovars <- array(gpData$phenoCovars[!rownames(gpData$phenoCovars) %in% which, , ], dim = phenoCovarsDim)
    dimnames(gpData$phenoCovars) <- phenoCovarsNames
  }

  # update pedigree
  if(!is.null(gpData$pedigree))
    if(!keepPedigree){
      gpData$pedigree <- subset(gpData$pedigree,!gpData$pedigree$ID %in% which)
      gpData$pedigree$Par1[gpData$pedigree$Par1 %in% which] <- 0
      gpData$pedigree$Par2[gpData$pedigree$Par2 %in% which] <- 0
      gpData$covar <- subset(gpData$covar,!gpData$covar$id %in% which)
    } else {
      gpData$covar[gpData$covar$id %in% which, c("phenotyped", "genotyped")] <- NA
    }
  # update covar
  gpData$covar <- subset(gpData$covar,!gpData$covar$id %in% which)
#  gpData$info$codeGeno <- FALSE

  return(gpData)
}
