\name{pairwiseLD}
\alias{pairwiseLD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pairwise LD between markers
}
\description{
Estimate pairwise Linkage Disequilibrium (LD) between markers measured as \eqn{r^2}{r2} using an object of class \code{gpData}. For the general case, a gateway to the software PLINK (Purcell et al. 2007) is established to estimate the LD. A within-R solution is available for marker data with only 2 genotypes, i.e. homozgous inbred lines. Return value is an object of class \code{LDdf} which is a \code{data.frame} with one row per marker pair or an object of class \code{LDMat} which is a \code{matrix} with all marker pairs. Additionally, the euclidian distance between position of markers is computed and returned.
}
\usage{
pairwiseLD(gpData, chr = NULL, type = c("data.frame", "matrix"),
           use.plink=FALSE, ld.threshold=0,
           ld.window=99999, rm.unmapped = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gpData}{
object of class \code{gpData} with elements \code{geno} and \code{map}
}                                 
  \item{chr}{
\code{numeric} scalar or vector. Return value is a list with pairwise LD of all markers for each chromosome in \code{chr}.
}                                      
  \item{type}{
\code{character}. Specifies the type of return value (see 'Value').
}
  \item{use.plink}{
\code{logical}. Should the software PLINK be used for the computation?
}
  \item{ld.threshold}{
\code{numeric}. Threshold for the LD to thin the output. Only pairwise LD>\code{ld.threshold} is reported when PLINK is used. This argument can only be used for \code{type="data.frame"}.
}
  \item{ld.window}{
\code{numeric}. Window size for pairwise differences which will be reported by PLINK (only for \code{use.plink=TRUE}) to thin the output dimensions. Only SNP pairs with a distance < \code{ld.window} are reported (default = 99999).
}
  \item{rm.unmapped}{
\code{logical}. Remove markers with unknown postion in \code{map} before using PLINK?
}
}
\details{
The function \code{write.plink} is called to prepare the input files and the script for PLINK. The executive PLINK file \code{plink.exe} must be available (e.g. in the working directory or through path variables). The function \code{pairwiseLD} calls PLINK and reads the results. The evaluation is performed separately for every chromosome. The measure for LD is \eqn{r^2}{r2}. This is defined as
\deqn{D= p_{AB} - p_Ap_B }{D = p(AB) - p(A)*p(B)} and
\deqn{r^2=\frac{D^2}{p_Ap_Bp_ap_b}}{r2=D2/p(A)p(B)p(a)p(b)}
where \eqn{p_{AB}}{p(AB)} is defined as the observed frequency of haplotype \eqn{AB}, \eqn{p_A=1-p_a} and \eqn{p_B=1-p_b} the observed frequencies of alleles \eqn{A} and \eqn{B}. 
If the number of markers is high, a threshold for the LD can be used to thin the output. In this case, only pairwise LD above the threshold is reported (argument \code{--ld-window-r2 in PLINK}).

Default PLINK options used
--no-parents
--no-sex 
--no-pheno
--allow-no-sex
--ld-window-kb 99999

}
\value{
For \code{type="data.frame"} an object of class \code{LDdf} with one element for each chromosome is returned. Each element is a \code{data.frame} with columns \code{marker1}, \code{marker2}, \code{r2} and \code{distance} for all \eqn{p(p-1)/2} marker pairs (or thinned, see 'Details').

For \code{type="matrix"} an object of class \code{LDmat} with one element for each chromosome is returned. Each element is a list of 2: a \eqn{p \times p}{p x p} \code{matrix} with pairwise LD and the corresponding \eqn{p \times p}{p x p} \code{matrix} with pairwise distances.

}
\references{
Hill WG, Robertson A (1968). Linkage Disequilibrium in Finite Populations. Theoretical
and Applied Genetics, 6(38), 226 - 231.

Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira MAR, Bender D, Maller J, Sklar P, de Bakker PIW, Daly MJ & Sham PC (2007)
PLINK: a toolset for whole-genome association and population-based linkage analysis. American Journal of Human Genetics, 81.
}
\author{
Valentin Wimmer
}

\seealso{
\code{\link{LDDist}}, \code{\link{LDMap}}
}
\examples{
\dontrun{
data(maize)
maizeC <- codeGeno(maize)
maizeLD <- pairwiseLD(maizeC,chr=1,type="data.frame")
}
}
