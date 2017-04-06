# adding markers to gpData object

add.markers <- function(gpData,geno,map=NULL){
  class(gpData$map) <- "data.frame"
  # check if markers are allready in data
  if (any(colnames(geno) %in% colnames(gpData$geno)) | any(rownames(map) %in% colnames(gpData$geno))){
    stop("some of the markers of ", substitute(geno)," are allready in ", substitute(gpData))
  }
  if(is.null(gpData$map) & !is.null(map)) stop("There is no map available for ", substitute(gpData), "!")
  if(!all(rownames(geno) %in% gpData$covar$id)) stop("You like to put new individuals into the data set!\nUse add.individuals() for that!")
  if(!all(rownames(map) %in% colnames(geno))) stop("There are markers in the map, which don't have information in ", substitute(geno), "!")
  # take names form map if available
  if(is.null(colnames(geno)) & !is.null(map))
    if(ncol(geno) == nrow(map)) colnames(geno) <- rownames(map) else
      stop("Check the colnames of", substitute(geno), " and the rownames of ", substitute(map), "!")

  # merge genotypic data
  nmiss <- sum(!rownames(gpData$geno) %in% rownames(geno))
  if(nmiss > 0){
    geno <- rbind(matrix(NA, nrow = nmiss, ncol = ncol(geno)), geno)
    rownames(geno)[1:nmiss] <- rownames(gpData$geno)[!rownames(gpData$geno) %in% rownames(geno)]
  }
  geno <- geno[match(rownames(gpData$geno), rownames(geno)), ]
  gpData$geno <- cbind(gpData$geno,geno)
  # first column as rownames and delete first column
  # merge map
  if(is.null(map)){
    #map <- gpData$map[1:ncol(geno),]
    map <- data.frame(chr=rep(NA,ncol(geno)),pos=rep(NA,ncol(geno)))
    rownames(map) <- colnames(geno)
  } else if(nrow(map) != ncol(geno)){
    map[colnames(geno)[!colnames(geno) %in% rownames(map)],] <- NA
  }
  map <- rbind(gpData$map,map)
  # same approach as in create.gpData
  map$sor <- substr(map$chr, nchar(as.character(map$chr)), nchar(as.character(map$chr)))
  if(!all(!unique(map$sor)[!is.na(unique(map$sor))] %in% 0:9)) map$sor <- 1
  map <- map[order(as.character(rownames(map))),]
  map <- orderBy(~sor+chr+pos,data=map)
  map$sor <- NULL
  gpData$map <- map
  gpData$geno <- gpData$geno[, match(rownames(gpData$map), colnames(gpData$geno))]
  gpData$info$codeGeno <- FALSE
  # sortcolumns in geno, too
  gpData$geno <- gpData$geno[,rownames(gpData$map)]
  # create new gpData object
  class(gpData$map) <- "GenMap"
  return(gpData)
}
