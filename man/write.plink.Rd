\name{write.plink}
\alias{write.plink}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prepare data for PLINK
}
\description{
Create input files and corresponding script for PLINK (Purcell et al. 2007) to estimate pairwise LD through function \code{pairwiseLD}.
}
\usage{
write.plink(gp, wdir = getwd(), prefix = paste(substitute(gp)),
            ld.threshold = 0, type = c("data.frame", "matrix"),
            ld.window=99999)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gp}{
\code{gpData} object with elements \code{geno} and \code{map}
}
  \item{wdir}{
\code{character}. Directory for PLINK input files
}
  \item{prefix}{
\code{character}. Prefix for PLINK input files.
}
  \item{ld.threshold}{
\code{numeric}. Threshold for the LD used in PLINK.
}
  \item{type}{
\code{character}. Specifies the type of return value for PLINK.
}
  \item{ld.window}{
\code{numeric}. Window size for pairwise differences which will be reported by PLINK (only for \code{use.plink=TRUE}) to thin the output dimensions. Only SNP pairs with a distance < \code{ld.window} are reported (default = 99999).
}
}
\value{
No value returned. Files \code{prefix.map}, \code{prefix.ped} and \code{prefixPlinkScript.txt} are created in the working directory
}
\references{
Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira MAR, Bender D, Maller J, Sklar P, de Bakker PIW, Daly MJ & Sham PC (2007)
PLINK: a toolset for whole-genome association and population-based linkage analysis. American Journal of Human Genetics, 81.
}
\author{
Valentin Wimmer
}

\seealso{
\code{\link{pairwiseLD}}
}
\examples{
\dontrun{write.plink(maize,type="data.frame")}
}

\keyword{manip}

