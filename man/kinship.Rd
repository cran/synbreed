\name{kin}
\alias{kin}
\title{
Relatedness based on pedigree or marker data
}
\description{     
This function implements different measures of relatedness between individuals in an object of class \code{gpData}: (1) Expected relatedness based on pedigree and (2) realized relatedness based on marker data. See 'Details'.  The function uses as first argument an object of class \code{gpData}. An argument \code{ret} controls the type of relatedness coefficient.
}
\usage{
kin(gpData,ret=c("add","kin","dom","gam","realized","realizedAB","sm","sm-smin"),
    DH=NULL, maf=NULL)
}
\arguments{
  \item{gpData}{
object of class \code{gpData}
}
  \item{ret}{
\code{character}. The type of relationship matrix to be returned. See 'Details'.
}
  \item{DH}{
\code{logical} vector of length \eqn{n}. \code{TRUE} or 1 if individual is a doubled-haploid (DH) line and \code{FALSE} or 0 otherwise. Only used for pedigree based relatedness coefficients
}
  \item{maf}{
    \code{numeric} vector of length equal the number of markers. Supply values for the \eqn{p_i}{pi} of each marker, which were used to correct the allele counts in \code{ret="realized"} and \code{ret="realizedAB"}. If not specified, \eqn{p_i}{pi} equals the minor allele frequency of each locus.
  }
}
\details{

\bold{Pedigree based relatedness (return arguments \code{"add"}, \code{"kin"}, \code{"dom"}, and \code{"gam"})}

Function \code{kin} provides different types of measures for pedigree based relatedness. An element \code{pedigree} must be available in the object of class \code{gpData}. In all cases, the first step is to build the gametic relationship. The gametic relationship is of order 2\eqn{n} as each individual has two alleles (e.g. individual \eqn{A} has alleles \eqn{A1} and \eqn{A2}). The gametic relationship
is defined as the matrix of probabilities that two alleles are identical by descent (IBD).  Note that the diagonal elements of the gametic relationship matrix are 1. The off-diagonals of individuals with unknown or unrelated parents in the pedigree are 0. If \code{ret="gam"} is specified, the gametic relationship matrix constructed by pedigree is returned.

The gametic relationship matrix can be used to construct other types of relationship matrices. If \code{ret="add"}, the additive numerator relationship matrix is returned. The additive relationship of individuals A (alleles \eqn{A1,A2}) and B (alleles \eqn{B1,B2}) is given by  the entries of the gametic relationship matrix
\deqn{0.5\cdot \left[(A1,B1) + (A1,B2) + (A2,B1) + (A2,B2)\right],}{0.5*[(A1,B1) + (A1,B2) + (A2,B1) + (A2,B2)],}
where \eqn{(A1,B1)} denotes the element [A1,B1] in the gametic relationship matrix. If \code{ret="kin"}, the kinship matrix is returned which is half of the additive relationship matrix.

If \code{ret="dom"}, the dominance relationship matrix is returned. The dominance relationship matrix between individuals A (\eqn{A1,A2}) and B (\eqn{B1,B2}) in case of no inbreeding is given by
\deqn{\left[(A1,B1) \cdot (A2,B2) + (A1,B2) \cdot (A2,B1)\right],}{[(A1,B1) * (A2,B2) + (A1,B2) * (A2,B1)],} 
where \eqn{(A1,C1)} denotes the element [A1,C1] in the gametic relationship matrix. 

\bold{Marker based relatedness (return arguments \code{"realized"},\code{"realizedAB"}, \code{"sm"}, and \code{"sm-smin"})}

Function \code{kin} provides different types of measures for marker based relatedness. An element \code{geno} must be available in the object of class \code{gpData}. Furthermore, genotypes must be coded by the number of copies of the minor allele, i.e. function \code{codeGeno} must be applied in advance.

If \code{ret="realized"}, the realized relatedness between individuals is computed according to the formulas in Habier et al. (2007) or vanRaden (2008)
\deqn{U = \frac{ZZ'}{2\sum p_i(1-p_i)}}{ZZ'/(2\sum pi(1-pi))} where \eqn{Z=W-P}, \eqn{W} is the marker matrix, \eqn{P} contains the allele frequencies multiplied by 2, \eqn{p_i}{pi} is the allele frequency of marker \eqn{i}, and the sum is over all loci.

If \code{ret="realizedAB"}, the realized relatedness between individuals is computed according to the formula in Astle and Balding (2009)
\deqn{U = \frac{1}{M} \sum \frac{(w_i-2p_i)(w_i-2p_i)'}{2p_i(1-p_i)}}{1/M sum((wi-2pi)(wi-2pi)'/(2pi(1-pi)))} where \eqn{w_i}{wi} is the marker genotype, \eqn{p_i}{pi} is the allele frequency at marker locus \eqn{i}, and \eqn{M} is the number of marker loci, and the sum is over all loci.

If \code{ret="sm"}, the realized relatedness between individuals is computed according to the simple matching coefficient (Reif et al. 2005). The simple matching coefficient counts the number of shared alleles across loci. It can only be applied to homozygous inbred lines, i.e. only genotypes 0 and 2. To account for loci that are alike in state but not identical by descent (IBD),  Hayes and Goddard (2008) correct the simple matching coefficient by the minimum of observed simple matching coefficients
\deqn{\frac{s-s_{min}}{1-s_{min}}}{s-smin/(1-smin)}
where \eqn{s} is the matrix of simple matching coefficients. This formula is used with argument \code{ret="sm-smin"}.
}
\value{
An object of class "relationshipMatrix".
}
\author{
Valentin Wimmer and Theresa Albrecht, with contributions by Yvonne Badke
}
\references{

Habier D, Fernando R, Dekkers J (2007). The Impact of Genetic Relationship information
on Genome-Assisted Breeding Values. Genetics, 177, 2389 -- 2397.

vanRaden, P. (2008). Efficient methods to compute genomic predictions.
Journal of Dairy Science, 91:4414 -- 4423.

Astle, W., and D.J. Balding (2009). Population Structure and Cryptic Relatedness in Genetic Association Studies. Statistical Science, 24(4), 451 -- 471.

Reif, J.C.; Melchinger, A. E. and Frisch, M. Genetical and mathematical properties of similarity and dissimilarity coefficients applied in plant breeding
and seed bank management. Crop Science, January-February 2005, vol. 45, no. 1, p. 1-7.

Rogers, J., 1972 Measures of genetic similarity and genetic distance. In Studies in genetics
VII, volume 7213. Univ. of Texas, Austin

Hayes, B. J., and M. E. Goddard. 2008. Technical note: Prediction of breeding values using marker derived relationship matrices.
J. Anim. Sci. 86

}     
\seealso{
\code{\link{plot.relationshipMatrix}} 
}
\examples{

#=========================
# (1) pedigree based relatedness
#=========================
\dontrun{
data(maize)
K <- kin(maize,ret="kin")
plot(K)
}

#=========================
# (2) marker based relatedness
#=========================
\dontrun{
data(maize)
U <- kin(codeGeno(maize),ret="realized")
plot(U)
}


### Example for Legarra et al. (2009), J. Dairy Sci. 92: p. 4660
id <- 1:17
par1 <- c(0,0,0,0,0,0,0,0,1,3,5,7,9,11,4,13,13)
par2 <- c(0,0,0,0,0,0,0,0,2,4,6,8,10,12,11,15,14)
ped <- create.pedigree(id,par1,par2)
gp <- create.gpData(pedigree=ped)

# additive relationship
A <- kin(gp,ret="add")
# dominance relationship
D <- kin(gp,ret="dom")
}
