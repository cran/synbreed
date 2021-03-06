\name{plot.LDdf}
\alias{plot.LDdf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot function for class \code{LDdf}
}
\description{
The function visualises wheter the LD between adjacent values or visualization of pairwise Linkage Disequilibrium (LD) estimates generated by function \code{pairwiseLD} versus marker distance. A single plot is generated for every chromosome.
}
\usage{
\method{plot}{LDdf}(x, gpData, plotType = "dist", dense = FALSE, nMarker = TRUE,
                    centr = NULL, chr = NULL, type = "p",breaks = NULL, n = NULL, 
                    file = NULL, fileFormat = "pdf", onefile = TRUE, colL = 2, 
                    colD = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Object of class \code{LDdf}, i.e the output of function \code{pairwiseLD} with argument \code{type="data.frame"}.
}
  \item{gpData}{
Object of class \code{gpData} with object \code{map}
}
  \item{plotType}{
You can decide, if you like to have a plot with the LD of the neighbouring markers (option \code{"neighbour"}), or you like to have a scatter plot of distance and LD (default option \code{"dist"}).
}
  \item{dense}{
For \code{plotType="neighbour"}, \code{logical}. Should density visualization for high-density genetic maps be used?
}
  \item{nMarker}{
For \code{plotType="neighbour"}, \code{logical}. Print number of markers for each chromosome?
}
  \item{centr}{
For \code{plotType="neighbour"}, \code{numeric} vector. Positions for the centromeres in the same order as chromosomes in \code{map}. If \code{"maize"}, centromere positions of maize in Mbp are used.
}
  \item{chr}{
For \code{plotType="dist"}, \code{numeric} scalar or vector. Return value is a plot for each chromosome in \code{chr}. Note: Remember to add in a batch-script one empty line for each chromosome, if you use more than one chromosome!
}  
  \item{type}{
For \code{plotType="dist"}, character string to specify the type of plot. Use \code{"p"} for a scatterplot, \code{"bars"} for stacked bars or \code{"nls"} for scatterplot together with nonlinear regression curve according to Hill and Weir (1988). 
}
  \item{breaks}{
For \code{plotType="dist"}, \code{list} containing breaks for stacked bars (optional, only for \code{type="bars"}). Components are \code{dist} with breaks for distance on x-axis and \code{r2} for breaks on for r2 on y-axis. By default, 5 equal spaced categories for dist and r2 are used.  
}
  \item{n}{
For \code{plotType="dist"}, \code{numeric}. Number of observations used to estimate LD. Only required for \code{type="nls"}.
}
  \item{file}{
Optionally a path to a file where the plot is saved to
}
  \item{fileFormat}{
\code{character}. At the moment two file formats are supported: pdf and png. Default is \code{"pdf"}.
}
  \item{onefile}{
\code{logical}. If \code{fileFormat = "pdf"} you can decide, if you like to have all graphics in one file or in multiple files.
}
  \item{colL}{
The color for the line if \code{type="nls"} is used. In other cases without a meaning.
}
  \item{colD}{
The color for the dots in the plot of \code{type="nls"} and \code{type="p"}
}
  \item{\dots}{
further graphical arguments for function \code{plot}
}
}
\details{
For more Details see at \code{\link{plotNeighbourLD}} or \code{\link{LDDist}} 
}
%%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%%}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
 Hans-Juergen Auinger
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{plotNeighbourLD}}, \code{\link{LDDist}}, \code{\link{plotGenMap}}, \code{\link{pairwiseLD}}
}
%%\examples{
%%
%%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
