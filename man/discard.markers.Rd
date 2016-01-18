\name{discard.markers}
\alias{discard.markers}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Subsets for objects of class gpData
}
\description{
The function produces subsets from an object of class \code{gpData} with reduced markers. Marker informartion will be discarded from elements \code{geno} and \code{map}
}
\usage{
discard.markers(gpData, which, whichNot=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gpData}{
object of class \code{gpData}
}
  \item{which}{
character vector identifying names of markers which get discarded in \code{geno} from a \code{gpData}-object.
}
  \item{whichNot}{
character vector identifying names of markers which get kept in \code{geno} from a \code{gpData}-object. Overwrites argument \code{which}!
}
}
\value{
Object of class \code{gpData}
}
\author{
Valentin Wimmer and Hans-Juergen Auinger
}

\seealso{
\code{\link{create.gpData}}, \code{\link{add.markers}}, \code{\link{add.individuals}}, \code{\link{discard.individuals}}
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

# remove unmapped SNP M5 (which has no postion in the map)
gp2 <- discard.markers(gp,"M5")
summary(gp2)

\dontrun{
# add one new DH line to maize data
library(synbreedData)
data(maize)

# delete markers
maize2 <- discard.individuals(maize,colnames(maize$geno)[1:50])
summary(maize2)
}
}
\keyword{manip}
