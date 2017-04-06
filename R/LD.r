pairwiseLD <- function(gpData,chr=NULL,type=c("data.frame","matrix"),use.plink=FALSE,
                       ld.threshold=0,ld.window=99999,rm.unmapped=TRUE, cores=1){
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
  multiCor <- function(x, use="everything", method = c("pearson", "kendall", "spearman"), cores=1){
    if(cores==1) cor(x, use=use, method=method) else {
      method <- match.arg(method)
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
    }
  }

  # catch errors
  if(is.null(gpData$geno)) stop("no genotypic data available")
  if(!gpData$info$codeGeno) stop("use function 'codeGeno' before")
  if(is.null(gpData$map))  stop("no map information available")
  type <- match.arg(type)
  if(ld.threshold>0 & type=="matrix") warning("'ld.threshold not used for type='matrix'")
  if(any(is.na(gpData$geno))) stop("no missing values allowed, try to impute using 'codeGeno'")

  # extract information from gpData  (only use mapped markers)
  if(rm.unmapped){
    mapped <- !(is.na(gpData$map$chr) | is.na(gpData$map$pos))
    gpData$geno <- gpData$geno[,mapped]
    gpData$map <- gpData$map[mapped,]
  } else {
    gpData$map$chr <- as.character(gpData$map$chr)
    gpData$map$chr[is.na(gpData$map$pos)] <- "NA"
    gpData$map$chr <- as.factor(gpData$map$chr)
  }
  linkageGroup <- as.character(gpData$map$chr)
  pos <- gpData$map$pos
  names(pos) <- rownames(gpData$map)

  # select chromosomes if 'chr' is specified
  lg <- unique(linkageGroup)
  if(!is.null(chr)){
    lg <- chr
    if(any(chr=="all")) stop("option chr='all' not yet possible") #linkageGroup <- rep("all",length(linkageGroup))
    # NOTE: positions must be ordered consequtively within chromosomes
  }

  # initialize return LD value data.frame list
  retList <- list()
  # initialize return LD value matrix list
  retMat <- list()

  # loop over all chromosomes (linkage groups)
  for (i in 1:length(lg)){
    if(use.plink){ # i.e. if there are 3 genotypes
      # call PLINK to compute the LD as r2
      sel <- rownames(gpData$map)[gpData$map$chr!= lg[i]]
      gpTEMP <- discard.markers(gpData,which=sel)
      pre <- paste("chr",lg[i],sep="")
      write.plink(gpTEMP,type=type,ld.threshold=ld.threshold,ld.window=ld.window,prefix=pre)
      system(paste("plink --script ",pre,"plinkScript.txt",sep=""))
      # distances between markers
      if(type=="matrix"){
        distance <- as.matrix(dist(pos[linkageGroup == lg[i]],diag=FALSE,upper=FALSE))
      }
      # read data from PLINK
      if(type=="matrix") {
        ld.r2 <- as.matrix(read.table(paste(pre,".ld",sep="")))
        colnames(distance) <- rownames(distance) <- colnames(ld.r2 ) <- rownames(ld.r2 ) <- names(pos)[linkageGroup == lg[i]]
      }
      if(type=="data.frame"){
        ld.r2.df.plink <- read.table(paste(pre,".ld",sep=""),header=TRUE,stringsAsFactors=FALSE)
        #distance <- abs(pos[ld.r2.df.plink$SNP_A]-pos[ld.r2.df.plink$SNP_B])
        ld.r2.df <- with(ld.r2.df.plink,data.frame(marker1=SNP_A,marker2=SNP_B,r2=R2,dist=abs(BP_A-BP_B),stringsAsFactors=FALSE))
      }
      ld.r <- sqrt(ld.r2)
    } # end if(use.plink)
      else{  # i.e. if there are 2 genotypes (e.g. DH lines)
      # read information from data
      markeri <- gpData$geno[,linkageGroup==lg[i]]
      p <- ncol(markeri)
      mn <- colnames(markeri)
      posi <- pos[linkageGroup==lg[i]]
      ld.r <- multiCor(markeri,method="spearman",use="pairwise.complete.obs", cores=cores)
      if(type=="data.frame"){
        ld.ri <- ld.r[lower.tri(ld.r)]
        # index vectors for LD data.frame
        rowi <- rep(1:p,times=(p:1)-1)
        coli <- p+1 - sequence(1:(p-1))
        coli <- coli[length(coli):1]
        # distance between markers
        disti <- abs(posi[rowi] - posi[coli])
        ld.r2.df <- data.frame(marker1=mn[rowi],
                               marker2=mn[coli],
                               r=ld.ri,
                               r2=ld.ri**2,
                               dist=disti,
                               stringsAsFactors=FALSE)
        if(lg[i]=="NA") {
          ld.rini <- multiCor(gpData$geno[,linkageGroup!=lg[i]],markeri,method="spearman",use="pairwise.complete.obs", cores=cores)
          ld.r2.dfini <- data.frame(marker1=rep(colnames(ld.rini), each=nrow(ld.rini)),
                                    marker2=rep(rownames(ld.rini), ncol(ld.rini)),
                                    r=as.numeric(ld.rini),
                                    r2=as.numeric(ld.rini**2),
                                    dist=rep(NA,ncol(ld.rini)*nrow(ld.rini)),
                                    stringsAsFactors=FALSE)
          ld.r2.df <- rbind(ld.r2.df, ld.r2.dfini)
        }
      }
      if(type=="matrix"){
        # matrix of distances
        distance <- as.matrix(dist(pos[linkageGroup == lg[i]],diag=FALSE,upper=FALSE))
        colnames(distance) <- rownames(distance) <- colnames(gpData$geno)[linkageGroup == lg[i]]
      }
    }
    # create dataset with information from above in a data.frame
    if(type=="data.frame") retList[[i]] <- ld.r2.df
      # and as a matrix
    if(type=="matrix") retMat$LD[[i]] <- ld.r**2           # omit lower/upper triangle?
    if(type=="matrix") retMat$distance[[i]] <- distance
  }
  if(type=="data.frame") names(retList) <- paste("chr", lg, sep="_")
  if(type=="matrix") names(retMat$LD) <- names(retMat$distance) <- paste("chr", lg, sep="_")

  # return values
  #if(type=="both") return(list(dataFrame=retList,matrix=retMat))
  if(type=="data.frame") {
    class(retList) <- "LDdf"
    return(retList)
  }
  if(type=="matrix"){
    class(retMat) <- "LDmat"
    return(retMat)
  }
}
