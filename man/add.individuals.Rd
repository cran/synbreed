\name{add.individuals}
\alias{add.individuals}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Add new individuals to objects of class gpData
}
\description{
This function extends an object of class \code{gpData} by adding new phenotypes, genotypes and pedigree.
}
\usage{
add.individuals(gpData, pheno = NULL, geno = NULL,
                pedigree = NULL, covar = NULL, repl=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gpData}{
object of class \code{gpData} to be updated
}
  \item{pheno}{
\code{data.frame} with new rows for phenotypes with \code{rownames} indicating individuals. For repeated values the ID should be stored in a column with name \code{"ID"}.
}
  \item{geno}{
\code{matrix} with new rows for genotypic data  with \code{rownames} indicating individuals
}
  \item{pedigree}{
\code{data.frame} with new rows for pedigree data
}
  \item{covar}{
\code{data.frame} with new rows for \code{covar} information with \code{rownames} indicating individuals
}
  \item{repl}{
The column of the pheno \code{data.frame} for the replicated measures. If the values are not repeated or  this column is named \code{"repl"} this argument is not needed.
}
}
\details{
\code{colnames} in \code{geno}, \code{pheno} and \code{pedigree} must match existing names in \code{gpData} object.
}
\value{
object of class \code{gpData} with new individuals
}
\author{
Valentin Wimmer
}

\seealso{
\code{\link{add.markers}}, \code{\link{discard.individuals}}
}
\examples{
# add one new DH line to maize data
library(synbreedData)
data(maize)
newDHpheno <- data.frame(Trait=1000,row.names="newDH")
# simulating genotypic data
newDHgeno <- matrix(sample(c(0,1),ncol(maize$geno),replace=TRUE),nrow=1)
rownames(newDHgeno) <- "newDH"
# new pedigree
newDHpedigree <- data.frame(ID="newDH",Par1=0,Par2=0,gener=0)
# new covar information
newDHcovar <- data.frame(family=NA,DH=1,tbv=1000,row.names="newDH")

# add individual
maize2 <- add.individuals(maize,newDHpheno,newDHgeno,newDHpedigree,newDHcovar)
summary(maize2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
