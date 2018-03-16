\name{add.gpData}
\alias{add.gpData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Join two \code{gpData} objects
}
\description{
Function for joining two \code{gpData} objects
}
\usage{
add.gpData(gpData1, gpData2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gpData1}{
A \code{gpData} object with at least elements \code{geno} and \code{map}
}
  \item{gpData2}{
Second \code{gpData} object with at least elements \code{geno} and \code{map}
}
}

\details{
The function writes a vcf-file. The format of the output is "GT". Other formats are not supported.
}
\value{
It is returned a \code{gpData} object, which contains \code{gpData1} and \code{gpData2}
}
\author{
Hans-Juergen Auinger
}

\seealso{ \code{\link{create.gpData}} \code{\link{codeGeno}}
}
\examples{
\dontrun{add.gpData(maize, maize)}
}
