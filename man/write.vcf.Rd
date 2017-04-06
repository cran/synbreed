\name{write.vcf}
\alias{write.vcf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prepare genotypic data in vcf-Format
}
\description{
Create vcf-file for miscellaneous applications. Within the package it is used to write files for beagle usage.
}
\usage{
write.vcf(gp, file, unphased=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gp}{
\code{gpData} object with elements \code{geno} and \code{map}
}
  \item{file}{
\code{character}. Filename for writing the file.
}
  \item{unphased}{
\code{logical}. The default is TRUE. Than the seperator between the alleles is \code{"/"}, and the
possible codings are \code{"0/0"} for \code{0} in the genotype matrix, \code{"0/1"} for \code{1} and
\code{"1/1"} for {2}. For getting a phased output, use \code{unphased=FALSE}. Than the seperator is
\code{"|"}. For hetercygous genotypes you have to change the 1 to -1 if you like to get the coding
\code{"1|0"}, So posible codings in this case are \code{"0|0"} for \code{0} in the genotype matrix,
\code{"0|1"} for \code{1}, \code{"1|0"} for \code{-1} and \code{"1|1"} for {2}.
}
}

\details{
The function writes a vcf-file. The format of the output is "GT". Other formats are not supported.
}
\value{
No value is returned. Function creates files \code{[prefix]ingput.bgl} with genotypic data in Beagle input format and \code{[prefix]marker.txt} with marker information used by Beagle.
}
\author{
Hans-Juergen Auinger
}

\seealso{
\code{\link{read.vcf2matrix}}, \code{\link{read.vcf2list}}, \code{\link{codeGeno}}
}
\examples{
map <- data.frame(chr=c(1,1,1,1,1,2,2,2,2),pos=1:9)
geno <-  matrix(sample(c(0,1,2,NA),size=10*9,replace=TRUE),nrow=10,ncol=9)
colnames(geno) <- rownames(map) <- paste("SNP",1:9,sep="")
rownames(geno) <- paste("ID",1:10+100,sep="")

gp <- create.gpData(geno=geno,map=map)
gp1 <- discard.markers(gp,rownames(map[map$chr!=1,]))
\dontrun{write.vcf(gp1, file="test.vcf")}
}

\keyword{manip}
