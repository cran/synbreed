\name{discard.markers}
\alias{discard.markers}
\alias{discard.individuals}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Subsets for objects of class gpData
}
\description{
Both functions could be used to get subsets from an object of class \code{gpData}. Use function \code{discard.markers} to discard markers from elements \code{geno} and \code{map} of an object of class \code{gpData}.  Use function \code{discard.individuals} to discard individuals from elements  \code{covar}, \code{pheno}, \code{geno} and \code{pedigree} of an object of class \code{gpData}.
}                      
\usage{
discard.markers(gpData, which)
discard.individuals(gpData, which, keepPedigree = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gpData}{
object of class \code{gpData}
}
  \item{which}{
character vector either identifying the \code{colnames} of markers in \code{geno} to discard (function \code{discard.markers}) or the \code{rownames} of individuals to discard (function \code{discard.individuals}).
}
  \item{keepPedigree}{
\code{logical}. Should the individual only be removed from elements \code{geno} and \code{pheno} but kept in the pedigree?
}
}
\value{
Object of class \code{gpData} 
}
\author{
Valentin Wimmer and Hans-Juergen Auinger
}

\seealso{
\code{\link{create.gpData}}
}
\examples{
\dontrun{
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

# remove unmapped SNP M5 (which has no postion in the map)
gp2 <- discard.markers(gp,"M5")
summary(gp2)

# discard genotypes with missing values in the marker matrix
gp3 <- discard.individuals(gp,names(which(rowSums(is.na(gp$geno))>0)))
summary(gp3)
}
 
}
\keyword{manip}
