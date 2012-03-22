\name{gpData2data.frame}
\alias{gpData2data.frame}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Merge of phenotypic and genotypic data
}
\description{
Create a \code{data.frame} out of phenotypic and genotypic data in object of class \code{gpData} by merging datasets using the common id. The shared data set could either include individuals with phenotypes and genotypes (default) or additional unphenotyped or ungenotyped individuals. In the latter cases, the missing observations are filled by \code{NA}'s.
}
\usage{
gpData2data.frame(gpData,trait=1,onlyPheno=FALSE,all.pheno=FALSE,
                  all.geno=FALSE,repl=NULL,phenoCovars=TRUE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gpData}{
object of class \code{gpData}
}
  \item{trait}{
\code{numeric} or \code{character}. A vector with the names or numbers of the trait that should be extracted from pheno. Default is \code{1}.
}
  \item{onlyPheno}{
scalar \code{logical}. Only return phenotypic data.
}
  \item{all.pheno}{
scalar \code{logical}. Include all individuals with phenotypes in the \code{data.frame} and fill the genotypic data with \code{NA}.
}
  \item{all.geno}{
scalar \code{logical}. Include all individuals with genotypes in the \code{data.frame} and fill the phenotypic data with \code{NA}.
}
  \item{repl}{
\code{character} or \code{numeric}. A vector which contains names or numbers of replication that should be drawn from the phenotypic values and covariates. Default is \code{NULL}, i.e. all values are used.
}
\item{phenoCovars}{
\code{logical}. If \code{TRUE}, columns with the phenotypic covariables are attached from element \code{phenoCovars} to the \code{data.frame}. Only required for repeated measurements.
}
  \item{...}{
further arguments to be used in function \code{reshape}. The argument \code{times} could be useful to rename the levels of the grouping variable (such as locations or environments).
}
}
\details{
Argument \code{all.geno} can be used to predict the genetic value of individuals without phenotypic records using the \code{BLR} package. Here, the genetic value of individuals with \code{NA} as phenotype is predicted by the marker profile.


For multiple measures, phenotypic data in object \code{gpData} is arranged with replicates in an \code{array}. With \code{gpData2data.frame} this could be reshaped to "long" format with multiple observations in one column. In this case, one column for the phenotype and 2 additional columns for the \code{id} and the levels of the grouping variable (such as replications, years of locations in multi-environment trials) are added.
}
\value{
A \code{data.frame} with the individuals names in the first column, the phenotypes in the next column(s) and the marker genotypes in subsequent columns.
}
\author{
Valentin Wimmer and Hans-Juergen Auinger
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{create.gpData}}, \code{\link{reshape}}
}
\examples{
# example data with unrepeated observations
set.seed(311)

# simulating genotypic and phenotypic data 
pheno <- data.frame(Yield = rnorm(12,100,5),Height=rnorm(12,100,1))
rownames(pheno) <- letters[4:15]
geno <- matrix(sample(c("A","A/B","B",NA),size=120,replace=TRUE,
prob=c(0.6,0.2,0.1,0.1)),nrow=10)
rownames(geno) <-  letters[1:10]
colnames(geno) <- paste("M",1:12,sep="")
# different subset of individuals in pheno and geno

# create 'gpData' object
gp <- create.gpData(pheno=pheno,geno=geno)
summary(gp)
gp$covar

# as data.frame with individuals with genotypes and phenotypes
gpData2data.frame(gp,trait=1:2)
# as data.frame with all individuals with phenotypes
gpData2data.frame(gp,1:2,all.pheno=TRUE)
# as data.frame with all individuals with genotypes
gpData2data.frame(gp,1:2,all.geno=TRUE)

# example with repeated observations
set.seed(311)

# simulating genotypic and phenotypic data 
pheno <- data.frame(ID = letters[1:10], Trait = c(rnorm(10,1,2),rnorm(10,2,0.2),
                    rbeta(10,2,4)), repl = rep(1:3, each=10))
geno <- matrix(rep(c(1,0,2),10),nrow=10)
colnames(geno) <- c("M1","M2","M3")
rownames(geno) <-  letters[1:10]

# create 'gpData' object
gp <- create.gpData(pheno=pheno,geno=geno, repeated="repl")

# reshape of phenotypic data and merge of genotypic data,
# levels of grouping variable loc are named "a", "b" and "c"
gpData2data.frame(gp,onlyPheno=FALSE,times=letters[1:3])
                            

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
