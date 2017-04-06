\name{read.vcf2matrix}
\alias{read.vcf2matrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read data of a vcf-file to a matrix
}
\description{
To easily read genomic data in vcf-Format to a matrix. Function \code{codeGeno} uses \code{read.vcf2matrix} with imputing by beagle.
}
\usage{
read.vcf2matrix(file, FORMAT = "GT", coding = c("allele", "ref"), IDinRow = TRUE, cores=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{
\code{character}. The name of the file which the data are to be read from.
}
  \item{FORMAT}{
\code{character}. The default is \code{"GT"}. If there are more formats in your vcf-file you can
decide which one you like to have in your output matrix.
}
  \item{coding}{
This option has only an effect with \code{FORMAT="GT"}. \code{allele} gives you back the alles as
defined as REF and ALT in your vcf-file. \code{ref} gives you back \code{"0"} for the reference
allele and \code{"1"} for the alternative allele.
}
  \item{IDinRow}{
\code{logical}. Default is \code{TRUE}, this means the genotypes are in the rows and the markers in
the column. For \code{FALSE} it is the other way round.
}
  \item{cores}{
\code{numeric}. Specifies the number of cores for parallel computing.
}
}
 \value{
A matrix (\code{\link[base]{matrix}}) containing a representation of the data in the file.
}
\author{
Hans-Juergen Auinger
}

\seealso{
\code{\link{write.vcf}}, \code{\link{read.vcf2list}}
}
\examples{
\dontrun{
library(synbreedData)
data(maize)
maize$info$map.unit <- "kb"
maize <- codeGeno(maize)
write.vcf(maize, "maize.vcf")
geno <- read.vcf2matrix("maize.vcf")
}
}
