\name{plot.relationshipMatrix}
\alias{plot.relationshipMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Heatmap for relationship Matrix
}
\description{
Visualization for objects of class \code{relationshipMatrix} using a heatmap of pairwise relatedness coefficients.
}
\usage{
\method{plot}{relationshipMatrix}(x, y=NULL, levelbreaks=NULL, axes=TRUE, cols=NULL, \ldots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Object of class \code{relationshipMatrix}
}
  \item{y}{
Optional for comparisons of objects of class \code{relationshipMatrix}
}
  \item{levelbreaks}{
\code{list} with one element for \code{x} and{y}. Define breaks in the color scheme of the plot.
If you make to many breaks, the color scheme repeats! If \code{y=NULL} this can be an
\code{vector}. If you like to have the same breaks or both relationship matrices, you can also use
just one vector.
}
  \item{axes}{
a \code{logical} value indicating whether axes should be drawn on the plot. Default is \code{TRUE}.
}
  \item{cols}{
a \code{list} with one element for each relationship matrix. Colors and the number of levelbreaks
should fit. But also if not, a plot is drawn. In case option \code{y=NULL}, \code{cols} can bei a
\code{vector}.
}
  \item{\dots}{
further graphical arguments passed to function \code{levelplot} in package \code{lattice}. To create equal colorkeys for two heatmaps, use \code{at=seq(from,to,length=9)}.
}
}

\author{
Valentin Wimmer and Hans-Juergen Auinger
}

\examples{
# small pedigree
ped <- simul.pedigree(gener=4,7)
gp <- create.gpData(pedigree=ped)
A <- kin(gp,ret="add")
plot(A)

# big pedigree
\dontrun{
library(synbreedData)
data(maize)
K <- kin(maize,ret="kin")
U <- kin(codeGeno(maize),ret="realized")/2
# equal colorkeys
plot(K,levelbreaks=seq(0,2,length=9))
plot(U,levelbreaks=seq(0,2,length=9))
}
}

\keyword{hplot}
