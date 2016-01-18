\name{create.pedigree}
\alias{create.pedigree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create pedigree object
}
\description{
This function can be used to create a \code{pedigree} object.
}
\usage{
create.pedigree(ID, Par1, Par2, gener=NULL,sex=NULL,add.ancestors=FALSE,unknown=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ID}{
vector of unique IDs identifying individuals
}
  \item{Par1}{
vector of IDs identifying parent 1 (with animals: sire)
}
  \item{Par2}{
vector of IDs identifying parent 2 (with animals: dam)
}
  \item{gener}{
vector identifying the generation. If \code{NULL} gener will be 0 for unknown parents and max(gener(Par1),gener(Par2))+1 for generations 1,... .
}
  \item{sex}{
vector identifying the sex (female=0 and male=1).
}
  \item{add.ancestors}{
\code{logical}. Add ancestors which do not occur in \code{ID} to the pedigree.
}
  \item{unknown}{
value for unknown or missing ancestors.
}
}
\details{
Missing values for parents in the pedigree should be coded \code{NA}. 0 is treaded as unknown, too.
}
\value{
An object of class \code{pedigree}. Column \code{gener} starts from 0 and pedigree is sorted by generation.
}
\author{
Valentin Wimmer
}

\seealso{
\code{\link{plot.pedigree}, \link{add.pedigree}}
}
\examples{
# example with 9 individuals
id <- paste("ID", 1:9, sep="0")
par1 <- paste("ID", c(0,0,0,0,1,1,1,4,7), sep="")
par2 <- paste("ID", c("","","","",2,3,2,5,8), sep="")
gener <- c(0,0,0,0,1,1,1,2,3)

# create pedigree object (using argument gener)
ped <- create.pedigree(id,par1,par2,gener,unknown=c("ID0", "ID"))
ped
plot(ped)

# create pedigree object (without using argument gener)
ped2 <- create.pedigree(id,par1,par2,unknown=c("ID0", "ID"))
ped2
}
\keyword{manip}
