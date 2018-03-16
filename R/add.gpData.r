#library(synbreedData)
#data(maize)
#gpData1 <- maize
#gpData2 <- maize
#gpData2$covar$id <- as.character(gpData2$covar$id)
#gpData2$covar$id[gpData2$covar$id %in% rownames(gpData2$geno)] <- paste(gpData2$covar$id[gpData2$covar$id %in% rownames(gpData2$geno)], 0, sep="")
#gpData2$pedigree$ID[gpData2$pedigree$ID %in% rownames(gpData2$geno)] <- paste(gpData2$pedigree$ID[gpData2$pedigree$ID %in% rownames(gpData2$geno)], 0, sep="")
#rownames(gpData2$geno) <- paste(rownames(gpData2$geno), 0, sep="")
#dimnames(gpData2$pheno)[[1]] <- paste(dimnames(gpData2$pheno)[[1]], 0, sep="")
#gpData2$pheno <- abind(gpData2$pheno, gpData2$pheno, along =2)
#dimnames(gpData2$pheno)[[2]] <- c("IDTrait", "Trait2")

add.gpData <- function(gpData1, gpData2){
  if(!is.null(gpData1$info$version)) stop(paste("Recode ", substitute(gpData1), "! You have used an old version to create/code ", substitute(gpData1), sep=""))
  if(substr(gpData1$info$version, 47, 50)<0.12) stop(paste("Recode ", substitute(gpData1), "! You have used an old version to create/code ", substitute(gpData1), sep=""))
  if(!is.null(gpData2$info$version)) stop(paste("Recode ", substitute(gpData1), "! You have used an old version to create/code ", substitute(gpData1), sep=""))
  if(substr(gpData2$info$version, 47, 50)<0.12) stop(paste("Recode ", substitute(gpData2), "! You have used an old version to create/code ", substitute(gpData2), sep=""))
  if(!is.null(gpData1$pheno))
    pheno1 <- gpData2data.frame(gpData1, onlyPheno=TRUE, trait=1:dim(gpData1$pheno)[2], stringsAsFactors=TRUE) else pheno1 <- NULL
  if(!is.null(gpData1$pheno))
    pheno2 <- gpData2data.frame(gpData2, onlyPheno=TRUE, trait=1:dim(gpData1$pheno)[2], stringsAsFactors=TRUE) else pheno2 <- NULL
  if(is.null(pheno1)){
    if(is.null(pheno2)) pheno <- NULL else{
      pheno <- pheno2
    }
  } else {
    if(is.null(pheno2)) pheno <- pheno1 else {
      pheno <- rbind(pheno1, pheno2)
    }
  }
  if(is.null(gpData1$geno)){
    if(is.null(gpData2$geno)) geno <- NULL else{
      geno <- gpData2$geno
    }
  } else {
    if(is.null(gpData2$geno)) geno <- gpData1$geno else {
      if(ncol(gpData1$geno)==ncol(gpData2$geno) & colnames(gpData1$geno)==colnames(gpData2$geno))
      geno <- rbind(gpData1$geno, gpData2$geno)
    }
  }
  return(0)
}
