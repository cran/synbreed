# read genomic prediction data
create.gpData <- function(pheno=NULL,geno=NULL,map=NULL,pedigree=NULL,family=NULL,covar=NULL,
                          reorderMap=TRUE,map.unit="cM",repeated=NULL,modCovar=NULL){
  
  # start with some checks on data
  # geno as matrix but not data.frame (storage) 
  if(!is.null(geno)){
    if(any(duplicated(rownames(geno)))) 
      stop(paste("In", substitute(geno), " are duplicated names of genotypes!"))
     if(class(geno) == "data.frame"){
        geno <- matrix(unlist(geno),nrow=nrow(geno),ncol=ncol(geno),dimnames=dimnames(geno))
        #if(any(duplicated(geno,MARGIN=1))) warning("individuals with dublicated genotypes")
     }
  }
  # test if positions in map are numeric
  if(!is.null(map$pos)) 
    if(!is.numeric(map$pos)) stop("Position informations have to be numeric values!")

 # test if chr in map is numeric or character
    if(!is.null(map$chr))    if(!(is.numeric(map$chr)|is.character(map$chr))) warning("Chromosome information should be numeric or character")


  # match geno and map
  if(!is.null(geno) & !is.null(map)){
    if(ncol(geno) != nrow(map)){  # different subsets of markers in geno and map
       if(is.null(colnames(geno)) & is.null(rownames(map))) stop("no marker names found in 'map' or 'geno'")
       else{ 
        # fill in gaps in map
        warning("not all markers in 'geno' mapped in 'map', gaps filled with 'NA' \n")
        map <- data.frame(chr=map$chr[match(colnames(geno),rownames(map))],pos=map$pos[match(colnames(geno),rownames(map))],row.names=colnames(geno))
       }
    }
    else{  # same subsets of markers in geno and map
      # missing colnames in geno
      if(is.null(colnames(geno))){  # no marker names in geno
        if (is.null(map)) stop("missing rownames in 'geno'")
        if(!is.null(rownames(map))){
          colnames(geno) <- rownames(map)
          warning("missing colnames in 'geno': assuming to be identical as rownames in 'map' \n")
        }
        else{ # default marker names
          warning("no marker names provide in 'geno' or 'map', using default names M1, M2, ... \n")
          colnames(geno) <- rownames(map) <- paste("M",1:ncol(geno),sep="")
        }
      }
      else{  # marker names in geno
      # missing rownames in map
         if(is.null(rownames(map))){
           rownames(map) <- colnames(geno)
           warning("missing rownames in 'map': assuming to be identical as colnames in 'geno' \n")  
         }
        # if(is.null(rownames(map)) & is.null(colnames(geno))){
        #  warning("missing marker names, setting default names M1, M2, ... ")
        #  colnames(geno) <- rownames(map) <- paste("M",1:ncol(geno),sep="")
        #} 
      }
    }
  }
  phenoCovars <- NULL
  attrModCovars <- NULL
  if(!is.null(pheno)){
    classList <- unlist(lapply(pheno, class))
    if(!all((classList[!names(classList) %in% repeated & !names(classList) %in% modCovar])[-1] %in% c("numeric", "integer"))) stop("Trait values have to be numeric!")
    # repeated measures? Use rownames of pheno as identifier for genotypes
    if(is.null(repeated)){
      add <- 10^ceiling(log10(nrow(pheno)))
      if(all(rownames(pheno) %in% 1:nrow(pheno))) rownames(pheno) <- add + as.numeric(rownames(pheno)) else add <- NULL
      if(dim(pheno)[2] ==1){# only a vector of traits
        phenoNames <- dimnames(pheno)
        arrPheno <- array(pheno[order(phenoNames[[1]]), ], dim = c(length(phenoNames[[1]]), 1, 1))
        dimnames(arrPheno) <- list(phenoNames[[1]][order(phenoNames[[1]])], phenoNames[[2]], "1")
      } else {# more than one trait, still unreplicated
        pheno <- pheno[order(rownames(pheno)), ]
        arrPheno <- array(as.matrix(pheno[, !(colnames(pheno) %in% modCovar)]), dim=c(dim(pheno[, !(colnames(pheno) %in% modCovar)]), 1))
        dimnames(arrPheno) <- list(rownames(pheno), colnames(pheno)[!(colnames(pheno) %in% modCovar)], "1")
      }
      if(!is.null(add)) dimnames(arrPheno)[[1]] <- as.numeric(dimnames(arrPheno)[[1]]) - add
      if(!is.null(modCovar)){
        arrModCovars <- array(1,  dim=c(dim(pheno[, colnames(pheno) %in% modCovar]), 1))
        dimnames(arrModCovars) <- list(dimnames(arrPheno)[[1]], colnames(pheno)[colnames(pheno) %in% modCovar], "1")
        for(i in colnames(pheno)[colnames(pheno) %in% modCovar])
          arrModCovars[, i, 1] <- pheno[, i]
      }
    } else {# a vector with replication idetifier is applied. The fist column is the identifier for genotypes
      dim3 <- data.frame(unique(pheno[, repeated]))
      colnames(dim3) <- repeated
      dim3 <- orderBy(as.formula(paste("~", paste(repeated, collapse = " + "))), data = dim3)
      for(i in 1:ncol(dim3)) dim3[, i] <- as.character(dim3[, i])
      rownam <- sort(unique(pheno[, 1]))
      if(!is.null(modCovar)) repeated <- unique(c(repeated, modCovar))
      arrPheno <- array(NA, dim = c(length(rownam), ncol(pheno)-(1+length(repeated)), nrow(dim3)))
      dimnames(arrPheno) <- list(rownam, (colnames(pheno)[!colnames(pheno) %in% repeated])[-1], as.character(apply(dim3, 1, paste, collapse = "_")))
      for(i in 1:nrow(dim3)){
        vec.bool <- apply(as.matrix(pheno[, colnames(dim3)]) == as.matrix(dim3[rep(i, nrow(pheno)), ]), 1, all)
        arrPheno[as.character(pheno[vec.bool, 1]), , i] <- as.matrix(pheno[vec.bool, (colnames(pheno)[!colnames(pheno) %in% repeated])[-1]])
      }
      if(!is.null(modCovar)){ # strip out covariates
        arrModCovars <- arrPheno[,rep(1, length(modCovar)), ]
        dimnames(arrModCovars)[[2]] <- colnames(pheno)[colnames(pheno) %in% modCovar]
        for(i in 1:nrow(dim3)){
          vec.bool <- apply(matrix(pheno[, colnames(dim3)] == dim3[rep(i, nrow(pheno)), ], ncol=ncol(dim3)), 1, all)
          arrModCovars[as.character(pheno[vec.bool, 1]), , i] <- as.matrix(pheno[vec.bool, colnames(pheno)[colnames(pheno) %in% modCovar]])
        }
      }
    }
    if(!is.null(modCovar)){ # take the correct class of the covariates
      phenoCovars <- arrModCovars
      attrModCovars <- classList[dimnames(arrModCovars)[[2]]]
      for(i in names(attrModCovars)){
        if(attrModCovars[i] == "numeric")
          attrModCovars[ i] <- "numeric"
        else
          attrModCovars[ i] <- "factor"
      }
    }
    pheno <- arrPheno
    rm(arrPheno)
  } 
  # match geno and pheno gpWheat1$pheno[1:5,]
  if(!is.null(geno) & !is.null(pheno)){
    if(is.null(dimnames(pheno)[[1]]) | is.null(rownames(geno))){
      if(dim(pheno)[1] == nrow(geno)){
        warning("assuming identical order of genotypes in 'pheno' and 'geno' \nControll the Output! There is no warranty of correctness!\n")  
        if(is.null(dimnames(pheno)[[1]])) dimnames(pheno)[[1]] <- rownames(geno)
        else rownames(geno) <- dimnames(pheno)[[1]]
        if(is.null(dimnames(pheno)[[1]]) & is.null(rownames(geno))) dimnames(pheno)[[1]] <- rownames(geno) <- paste("ID",10^ceiling(log10(nrow(geno)))+1:nrow(geno),sep="")
    }
    # now geno and pheno have rownames
    else stop("missing rownames for 'pheno' and 'geno'")
    }
   } 
    # sort geno and pheno by rownames (alphabetical order)
    if(!is.null(geno)){
      if(all(row.names(geno) %in% 1:nrow(geno))) 
        geno <- geno[order(as.numeric(row.names(geno))), ]
      else
        geno <- geno[order(row.names(geno)),]
    }

  # sort markers by chromosome and position within chromosome
  if(!is.null(map)){
    if(any(colnames(map) != c("chr","pos"))) stop("colnames of", substitute(map),"must be 'chr' and 'pos'")
    if (reorderMap){
      # first order in alphabetical oder (important for SNPs with the same position)
      map$sor <- substr(map$chr, nchar(as.character(map$chr)), nchar(as.character(map$chr)))
      if(!all(!unique(map$sor)[!is.na(unique(map$sor))] %in% 0:9)) map$sor <- 1
      map <- map[order(as.character(rownames(map))),]
      map <- orderBy(~sor+chr+pos,data=map)
      map$sor <- NULL
      # sortcolumns in geno, too
      geno <- geno[,rownames(map)]
    }
    class(map) <- c("GenMap", "data.frame")
  }

  # return object
  # geno as matrix
  if(!is.null(geno)) geno <- data.matrix(geno,TRUE)
  obj <- list(covar=NULL,pheno=pheno,geno=geno,map=map,pedigree=pedigree,phenoCovars=phenoCovars)
  
  # add information to element covar
  # sort all available individuals
  ids <- sort(unique(c(dimnames(obj$pheno)[[1]],rownames(obj$geno),as.character(obj$pedigree$ID)))) 
  if(all(ids %in% 1:length(ids))) ids <- sort(as.numeric(ids))

  if(is.null(covar)) obj$covar <- data.frame(id=ids,stringsAsFactors=FALSE)
  else obj$covar$id <- ids 

  obj$covar$phenotyped <- obj$covar$id %in% dimnames(obj$pheno)[[1]]
  obj$covar$genotyped <- obj$covar$id %in% rownames(obj$geno)
  
  # family information for genotyped indviduals  
  if(!is.null(family)){
    colnames(family)[1] <- "family"
    obj$covar <- merge(obj$covar,family,by.x="id",by.y=0,all=TRUE)
  } else obj$covar$family <- NA
  
  # add covar from arguments, if available 
  if(!is.null(covar)){
    if(is.null(rownames(covar))) stop("missing rownames in covar")    
    # do not use any existing columns named 'genotyped', 'phenotyped' or 'id'
    covar <- covar[!colnames(covar) %in% c("genotyped","phenotyped","id","family")]
    # merge with existing data
    if(!is.null(covar)) obj$covar <- merge(obj$covar,covar,by.x=1,by.y=0,all=TRUE)
    else  obj$covar <- obj$covar
  }
  
  # further information
  obj$info$map.unit <- map.unit
  obj$info$codeGeno <- FALSE
  obj$info$attrPhenoCovars <- attrModCovars

  # set class of sub-object pedigree
  if(!is.null(obj$pedigree)){
    if(!any(class(obj$pedigree)=="pedigree")) warning("object pedigree should be of class 'pedigree'")
  }
  # return object of class 'gpData'
  class(obj) <- "gpData"
  return(obj)
}
