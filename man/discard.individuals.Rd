\name{discard.individuals}
\alias{discard.individuals}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Subsets for objects of class gpData
}
\description{
The function produce subsets from an object of class \code{gpData} with reduced individuals. Individual information will be discarded from elements \code{geno}, \code{pheno}, \code{covar} and \code{pedigree}.
}
\usage{
discard.individuals(gpData, which, keepPedigree = FALSE, whichNot=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gpData}{
object of class \code{gpData}
}
  \item{which}{
character vector identifying names of individuals get discarded from a \code{gpData}-object.
}
  \item{keepPedigree}{
\code{logical}. Should the individual only be removed from elements \code{geno} and \code{pheno} but kept in the pedigree?
}
  \item{whichNot}{
character vector identifying names of individuals get kept in a \code{gpData}-object. Overwrites argument \code{which}!
}
}
\value{
Object of class \code{gpData}
}
\author{
Valentin Wimmer and Hans-Juergen Auinger
}

\seealso{
\code{\link{create.gpData}}, \code{\link{add.individuals}}, \code{\link{add.markers}}, \code{\link{discard.markers}}
}
\examples{

# example data
set.seed(311)
pheno <- data.frame(Yield = rnorm(10,200,5),Height=rnorm(10,100,1))
rownames(pheno) <- letters[1:10]
geno <- matrix(sample(c("A","A/B","B",NA),size=120,replace=TRUE,
prob=c(0.6,0.2,0.1,0.1)),nrow=10)
rownames(geno) <-  letters[1:10]
colnames(geno) <- paste("M",1:12,sep="")
# one SNP is not mapped (M5)
map <- data.frame(chr=rep(1:3,each=4),pos=rep(1:12))
map  <- map[-5,]
rownames(map) <- paste("M",c(1:4,6:12),sep="")
gp <- create.gpData(pheno=pheno,geno=geno,map=map)
summary(gp)

# discard genotypes with missing values in the marker matrix
gp3 <- discard.individuals(gp,names(which(rowSums(is.na(gp$geno))>0)))
summary(gp3)

\dontrun{
# add one new DH line to maize data
library(synbreedData)
data(maize)

# delete individual
maize2 <- discard.individuals(maize,rownames(maize$geno)[1:10])
summary(maize2)

}
}
\keyword{manip}
