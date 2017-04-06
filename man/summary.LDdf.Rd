\name{summary.LDdf}

\Rdversion{1.1}

\alias{summary.LDdf}
\alias{summary.LDmat}

\title{
Summary for LD objects
}
\description{
\code{Summary} method for class ''LDdf''  and  ''LDmat''
}
\usage{
\method{summary}{LDdf}(object,cores=1,...)
\method{summary}{LDmat}(object,cores=1,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
object of class \code{LDdf} or \code{LDmat} which is the output of function \code{pairwiseLD} and argument \code{type="data.frame"}  or   \code{type="matrix"}
}
  \item{cores}{
\code{numeric}. Specifies the number of cores for parallel computing.
}
  \item{...}{
not used
}
}
\details{
Returns for each chromosome: Number of markers; mean, minimum and maximum LD measured as \eqn{r^2}{r2};  fraction of markers with \eqn{r^2>0.2}{r2>0.2};  average and maximum distance of number of markers
}
\author{
Valentin Wimmer and Hans-Juergen Auinger
}

\seealso{
\code{\link{pairwiseLD}}, ~~~
}
\examples{
  \dontrun{
  library(synbreed)
  data(maize)
  maizeC <- codeGeno(maize)
  maizeLD <- pairwiseLD(maizeC,chr=1:10,type="data.frame")
  maizeLDm <- pairwiseLD(maizeC,chr=1:10,type="matrix")
  summary(maizeLD)
  summary(maizeLDm)
  }
}

\keyword{methods}
