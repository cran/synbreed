\name{plotNeighbourLD}
\alias{plotNeighbourLD}
\title{
Plot neighbour linkage disequilibrium
}
\description{
A function to visualize Linkage Disequilibrium estimates between adjacent markers.
}
\usage{
plotNeighbourLD(LD, gpData, dense=FALSE, nMarker = TRUE,
                centr=NULL, file=NULL, fileFormat="pdf", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{LD}{
object of class \code{LDmat}, i.e the output of function \code{pairwiseLD} using argument \code{type="matrix"}.
}
  \item{gpData}{
object of class \code{gpData} with object \code{map} or a \code{data.frame} with columns 'chr' (specifying the chromosome of the marker) and 'pos' (position of the marker within chromosome measured with genetic or physical distances)
}
  \item{dense}{
\code{logical}. Should density visualization for high-density genetic maps be used?
}
  \item{nMarker}{
\code{logical}. Print number of markers for each chromosome?
}
  \item{centr}{
\code{numeric} vector. Positions for the centromeres in the same order as chromosomes in \code{map}. If \code{"maize"}, centromere positions of maize in Mbp are used.
}
  \item{file}{
Optionally a path to a file where the plot is saved to
}
  \item{fileFormat}{
\code{character}. At the moment two file formats are supported: pdf and png. Default is \code{"pdf"}.
}
  \item{\dots}{
further graphical arguments for function \code{plot}
}

}
\details{
The plot is similar to \code{plotGenMap} with the option \code{dense=TRUE}, but here the LD between adjacent markers is plotted along the chromosomes.
}
\value{
Plot of neighbour LD along each chromosome. One chromosome is displayed from the first to the last marker.
}

\author{
Theresa Albrecht and Hans-Juergen Auinger
}

\seealso{
\code{\link{plotGenMap}}, \code{\link{pairwiseLD}}
}
\examples{
\dontrun{
library(synbreedData)
data(maize)
maize2 <-codeGeno(maize)
LD <- pairwiseLD(maize2,chr=1:10,type="matrix")
plotNeighbourLD(LD,maize2,nMarker=FALSE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}

