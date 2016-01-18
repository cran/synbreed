\name{add.pedigree}
\alias{add.pedigree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Merge \code{pedigree} object}
\description{
This function can be used to add some pedigree information to a existing \code{pedigree} object.
}
\usage{
add.pedigree(ped, IDadd, add.ancestors = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ped}{
\code{pedigree} object
}
  \item{IDadd}{
\code{pedigree} object
}
  \item{add.ancestors}{
\code{logical}. Add ancestors which do not occur in \code{ID} to the pedigree.
}
}
\details{
Missing values for parents in the pedigree should be coded with 0 for numeric ID or \code{NA} for character ID.
}
\value{
An object of class \code{pedigree}. Column \code{gener} starts from 0 and pedigree is sorted by generation.
}
\author{
Hans-Juergen Auinger
}

\seealso{
\code{\link{plot.pedigree}, \link{create.pedigree}}
}
\examples{
# example with 9 individuals
id <- paste("ID", 1:9, sep="0")
par1 <- paste("ID", c("","","","",1,1,1,4,7), sep="0")
par2 <- paste("ID", c("","","","",2,3,2,5,8), sep="0")
ped1 <- create.pedigree(id,par1,par2,unknown="ID0")

# create 2nd pedigree object
Id <- paste("ID", 10:16, sep="")
Par1 <- paste("ID", c("","",1,1,6,7,7), sep="0")
Par2 <- paste("ID", c("","",10,"08","09",11,14), sep="")
ped2 <- create.pedigree(Id,Par1,Par2,unknown=c("ID0", "ID"))
ped2

ped <- add.pedigree(ped1, ped2)
ped

}
