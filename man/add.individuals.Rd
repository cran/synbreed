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

\seealso{\code{\link{add.markers}}, \code{\link{discard.individuals}}
}
\examples{
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

#new phenotypic data
newPheno <- data.frame(Yield=200,Height=100,row.names="newLine")
# simulating genotypic data
newGeno <- matrix(sample(c("A","A/B","B"),ncol(gp$geno),replace=TRUE),nrow=1)
rownames(newGeno) <- "newLine"
# new pedigree
newPedigree <- create.pedigree(ID="newLine",Par1=0,Par2=0,gener=0)

gp2 <- add.individuals(gp,pheno=newPheno,geno=newGeno,pedigree=newPedigree)


\dontrun{
# add one new DH line to maize data
library(synbreedData)
data(maize)
newDHpheno <- data.frame(Trait=1000,row.names="newDH")
# simulating genotypic data
newDHgeno <- matrix(sample(c(0,1),ncol(maize$geno),replace=TRUE),nrow=1)
rownames(newDHgeno) <- "newDH"
# new pedigree
newDHpedigree <-  create.pedigree(ID="newDH",Par1=0,Par2=0,gener=0)
# new covar information
newDHcovar <- data.frame(family=NA,DH=1,tbv=1000,row.names="newDH")

# add individual
maize2 <- add.individuals(maize,newDHpheno,newDHgeno,newDHpedigree,newDHcovar)
summary(maize2)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
