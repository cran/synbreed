write.beagle <- function(gp,wdir=getwd(),prefix){

       # information from arguments
       geno <- gp$geno
       n <- nrow(geno)
       M <- ncol(geno)

       # prepare input file input.bgl
       if (any(grep(" ",rownames(geno)))) stop("no blanks allowed in rownames(geno) when running beagle")
       firstLine <- c("I","id",rep(rownames(geno),each=2))
       # assume coding AA, AB, BB
       allele1 <- substr(geno,1,1)
       allele2 <- substr(geno,2,2)
       hap <- rep(NA,2*n)
       hap[seq(1,(2*n*M),by=2)] <- allele1
       hap[seq(2,(2*n*M),by=2)] <- allele2
       # markers = rows
       haplotypes <-  matrix(hap,nrow=M,ncol=2*n,byrow=TRUE)
       # write input file
       input.bgl <- rbind(firstLine,cbind("M",colnames(geno),haplotypes))
       write.table(input.bgl,file=file.path(wdir,paste(prefix,"input.bgl",sep="")),quote=FALSE,col.names=FALSE,row.names=FALSE)

       # write marker file
       #alleles <- apply(gp$geno,2,unique)
       # make all equal length
       #alleles <- lapply(alleles,equal.length,max(unlist(lapply(alleles,length))))
       #alleles <- matrix(unlist(alleles),byrow=TRUE,nrow=M)
       write.table(cbind(rownames(gp$map),gp$map$pos,"A B"),file=file.path(wdir,paste(prefix,"marker.txt",sep="")),quote=FALSE,col.names=FALSE,row.names=FALSE)

}
