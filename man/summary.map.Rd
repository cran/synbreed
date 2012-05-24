\name{summaryGenMap}
\alias{summaryGenMap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary of marker map information
}
\description{
This function can be used to summarize information from a marker map in an object of class \code{gpData}.  Return value is a \code{data.frame} with one row for each chromosome and one row summarizing all chromosomes.
}
\usage{
summaryGenMap(map)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{map}{
\code{data.frame} with columns \code{chr} and \code{pos} or a \code{gpData} object with element \code{map}
}
}
\details{
Summary statistics of differences are based on euclidian distances between markers with non-missing position in \code{map}, i.e. \code{pos!=NA}.
}
\value{
A \code{data.frame} with one row for each chromosome and the intersection of all chromosomes and columns
  \item{noM}{number of markers}
  \item{range}{range of positions, i.e. difference between first and last marker}
  \item{avDist}{avarage distance of markers}
  \item{maxDist}{maximum distance of markers}
  \item{minDist}{minimum distance of markers}
}
\author{
Valentin Wimmer
}
\seealso{
\code{\link{create.gpData}}
}
\examples{
\dontrun{
data(maize)
summaryGenMap(maize)
}
}
