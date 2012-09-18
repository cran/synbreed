\name{plotGenMap}
\alias{plotGenMap}
\alias{plot.GenMap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot marker map
}
\description{
A function to visualize low and high-density marker maps.
}
\usage{
\method{plot}{GenMap}(x, dense = FALSE, nMarker = TRUE, bw=1, centr=NULL, file=NULL, fileFormat="pdf",...)

plotGenMap(map, dense = FALSE, nMarker = TRUE, bw=1, centr=NULL, file=NULL, fileFormat="pdf",...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
object of class \code{GenMap}, i. e. the \code{map} object in a \code{gpData}-object
}
  \item{map}{
object of class \code{gpData} with object \code{map} or a \code{data.frame} with columns 'chr' (specifying the chromosome of the marker) and 'pos' (position of the marker within chromosome measured with genetic or physical distances)
}                                                                                  
  \item{dense}{
\code{logical}. Should density visualization for high-density genetic maps be used?
}
  \item{nMarker}{
\code{logical}. Print number of markers for each chromosome?
}
  \item{bw}{
\code{numeric}. Bandwidth to use for  \code{dense=TRUE} to control the resolution (default = 1 [map unit]).
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
In the low density plot, the unique positions of markers are plotted as horizontal lines. In the  high-density plot, the distribution of the markers is visualized as a heatmap of density estimation together with a color key. In this case, the number of markers within an interval of equal bandwidth \code{bw} is counted. The high density plot is typically useful if the number of markers exceeds 200 per chromosome on average.
}
\value{
Plot of the marker positions within  each chromosome. One chromosome is displayed from the first to the last marker.
}

\author{
Valentin Wimmer and Hans-Juergen Auinger
}


\seealso{
\code{\link{create.gpData}}
}
\examples{
\dontrun{
# low density plot
data(maize)
plotGenMap(maize)

# high density plot  
data(mice)
plotGenMap(mice,dense=TRUE,nMarker=FALSE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
