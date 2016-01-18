\name{plot.pedigree}
\alias{plot.pedigree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Visualization of pedigree
}
\description{
A function to visualize pedigree structure by a graph using the \code{igraph} package. Each genotype is represented as vertex and direct offsprings are linked by an edge.
}
\usage{
\method{plot}{pedigree}(x, effect = NULL,\dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
object of class \code{pedigree} or object of class \code{gpData} with element \code{pedigree}
}
  \item{effect}{
vector of length \code{nrow(pedigree)} with effects to plot on the x axis
}
  \item{...}{
Other arguments for function \code{igraph.plotting}
}
}
\details{
 The pedigree is structured top to bottom. The first generation is printed in the first line. Links over more than one generation are possible as well as genotypes with only one (known) parent. Usually, no structure in one generation is plotted. If an \code{effect} is given, the genotypes are ordered by this effect in the horizontal direction and a labeled axis is plotted at the bottom.
}
\value{
A named graph visualizing the pedigree structure. Color is used to distinguish sex.
}

\author{
Valentin Wimmer and Hans-Juergen Auinger
}
\note{
This function uses the plotting method for graphs in the library \code{igraph}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{create.pedigree}}, \code{\link{simul.pedigree}}
}
\examples{
id <- paste("ID", 1:9, sep="0")
par1 <- paste("ID", c("","","","",1,1,1,4,7), sep="0")
par2 <- paste("ID", c("","","","",2,3,2,5,8), sep="0")
ped1 <- create.pedigree(id,par1,par2,unknown="ID0")
ped1
plot(ped1)

# create 2nd pedigree object
Id <- paste("ID", 10:16, sep="")
Par1 <- paste("ID", c("","",1,1,6,7,7), sep="0")
Par2 <- paste("ID", c("","",10,"08","09",11,14), sep="")
ped2 <- create.pedigree(Id,Par1,Par2,unknown=c("ID0", "ID"))
ped2

ped <- add.pedigree(ped1, ped2)
plot(ped)
}
\keyword{hplot}
