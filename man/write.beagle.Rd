\name{write.beagle}
\alias{write.beagle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prepare genotypic data for Beagle
}
\description{
Create input file for Beagle software (Browning and Browning 2009) from an object of class \code{gpData}. This function is created for usage within function \code{codeGeno} to impute missing values.
}
\usage{
write.beagle(gp, wdir = getwd(), prefix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gp}{
\code{gpData} object with elements \code{geno} and \code{map}
}
  \item{wdir}{
\code{character}. Directory for Beagle input files
}
  \item{prefix}{
\code{character}. Prefix for Beagle input files.
}
}
\details{
The Beagle software must be used chromosomewise. Consequently, \code{gp} should contain only data from one chromosome (use \code{discard.markers, see Examples}).
}
\value{
No value is returned. Function creates files \code{[prefix]ingput.bgl} with genotypic data in Beagle input format and \code{[prefix]marker.txt} with marker information used by Beagle.
}
\references{
B L Browning and S R Browning (2009) A unified approach to genotype imputation and haplotype phase inference for large data sets of trios and unrelated individuals. Am J Hum Genet 84:210-22
}
\author{
Valentin Wimmer
}

\seealso{
\code{\link{codeGeno}}
}
\examples{
map <- data.frame(chr=c(1,1,1,1,1,2,2,2,2),pos=1:9)
geno <-  matrix(sample(c(0,1,2,NA),size=10*9,replace=TRUE),nrow=10,ncol=9)
colnames(geno) <- rownames(map) <- paste("SNP",1:9,sep="")
rownames(geno) <- paste("ID",1:10+100,sep="")

gp <- create.gpData(geno=geno,map=map)
gp1 <- discard.markers(gp,rownames(map[map$chr!=1,]))
\dontrun{write.beagle(gp1,prefix="test")}
}

\keyword{manip}

