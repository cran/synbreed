\name{[.relationshipMatrix}
\alias{[.relationshipMatrix}
\docType{data}
\title{
Extract or replace part of relationship matrix
}
\description{
Extract or replace part of an object of class \code{relationshipMatrix}.
}
\usage{
\method{[}{relationshipMatrix}(x,\dots)
}
\arguments{
  \item{x}{
object of class ''relationshipMatrix''
}
  \item{...}{
indices
}
}
\examples{
\dontrun{
data(maize)
U <- kin(codeGeno(maize),ret="realized")
U[1:3,1:3]
}
}
