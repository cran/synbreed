write.vcf <- function(gp,file){

  if(is.null(nrow(gp$geno)) | is.null(ncol(gp$geno))){
    print(str(gp$geno))
    stop("Wrong genotypic information!")
  }
  if(!gp$info$map.unit %in% c("bp", "kb", "Mb")) stop("You need basepairs as map positions to write a vcf-file!")
  if(gp$info$map.unit =="kb") gp$map$pos <- gp$map$pos * 1000
  if(gp$info$map.unit =="Mb") gp$map$pos <- gp$map$pos * 1000000
  if(any((gp$map$pos-round(gp$map$pos, digits=0)) >1e-6)) stop("Your map positions and the unit of the position do not fit!")
  geno <- as.data.frame(t(gp$geno), stringsAsFactors=FALSE)
  geno[geno==0] <- "0|0"
  geno[geno==1] <- "0|1"
  geno[geno==2] <- "1|1"
  bgl <- cbind(data.frame(CHROM=paste("chr",gp$map$chr, sep=""), POS=gp$map$pos,ID=rownames(gp$map), REF="A", ALT="G", QUAL=".", FILTER="PASS", INFO=".", FORMAT="GT", stringsAsFactors=FALSE),
               geno)
  if (any(grep(" ",colnames(geno)))) stop("no blanks allowed in IDs when running beagle")
  cat(file=file, '##fileformat=VCFv4.1\n##filedate=20120310\n##source="write.vcf of R-synbreed"\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n#')
  cat(file=file, paste(colnames(bgl), collapse="\t"), "\n", append=TRUE)
  write.table(bgl, file=file,
              quote=FALSE, col.names=FALSE, row.names=FALSE, append=TRUE, sep="\t", na=".|.")

}

read.vcf2matrix <- function(file, FORMAT="GT", coding=c("allele","ref"), IDinRow=TRUE){

  coding <- match.arg(coding)
  cnt=0
  while(scan(file=file, what="character", skip=cnt, nlines=1, quiet=TRUE)[1] !="#CHROM") cnt <- cnt+1
  Mnames <- scan(file, what="character", skip=cnt, nlines=1, quiet=TRUE)
  geno <- read.table(file=file, sep="\t", header=FALSE, skip=cnt+1, stringsAsFactors=FALSE)
  colnames(geno) <- Mnames
  rownames(geno) <- geno$ID
  ref <- geno$REF; alternative <- geno$ALT
  form <- unlist(strsplit(geno$FORMAT, ":"))
  geno <- geno[, !colnames(geno) %in% c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")]
  geno[1:nrow(geno), 1:ncol(geno)] <- unlist(lapply(geno, strsplit, ":"))[rep(form, ncol(geno)) == FORMAT]
  if(FORMAT == "GT" & coding == "allele"){
    geno[geno=="0|0"] <- rep(paste(ref,ref, sep="|"), ncol(geno))[geno=="0|0"]
    geno[geno=="1|0"] <- rep(paste(alternative,ref, sep="|"), ncol(geno))[geno=="1|0"]
    geno[geno=="0|1"] <- rep(paste(ref,alternative, sep="|"), ncol(geno))[geno=="0|1"]
    geno[geno=="1|1"] <- rep(paste(alternative,alternative, sep="|"), ncol(geno))[geno=="1|1"]
  }
  if(IDinRow) geno <- t(geno)
  return(geno)

}

read.vcf2list <- function(file, FORMAT="GT", coding=c("allele","ref"), IDinRow=TRUE){

  coding <- match.arg(coding)
  cnt=0
  while(scan(file=file, what="character", skip=cnt, nlines=1, quiet=TRUE)[1] !="#CHROM") cnt <- cnt+1
  Mnames <- scan(file, what="character", skip=cnt, nlines=1, quiet=TRUE)
  geno <- read.table(file=file, sep="\t", header=FALSE, skip=cnt+1, stringsAsFactors=FALSE)
  colnames(geno) <- Mnames
  rownames(geno) <- geno$ID
  ref <- geno$REF; alternative <- geno$ALT
  form <- unlist(strsplit(geno$FORMAT, ":"))
  map <- geno[, c("#CHROM", "POS")]
  colnames(map) <- c("chr", "pos")
  class(map) <- c("GenMap", "data.frame")
  geno <- geno[, !colnames(geno) %in% c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")]
  geno[1:nrow(geno), 1:ncol(geno)] <- unlist(lapply(geno, strsplit, ":"))[rep(form, ncol(geno)) == FORMAT]
  if(FORMAT == "GT" & coding == "allele"){
    geno[geno=="0|0"] <- rep(paste(ref,ref, sep="|"), ncol(geno))[geno=="0|0"]
    geno[geno=="1|0"] <- rep(paste(alternative,ref, sep="|"), ncol(geno))[geno=="1|0"]
    geno[geno=="0|1"] <- rep(paste(ref,alternative, sep="|"), ncol(geno))[geno=="0|1"]
    geno[geno=="1|1"] <- rep(paste(alternative,alternative, sep="|"), ncol(geno))[geno=="1|1"]
  }
  if(IDinRow) geno <- t(geno)
  return(list(geno=geno, map=map))

}
