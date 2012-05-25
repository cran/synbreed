\name{codeGeno}
\alias{codeGeno}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Recode genotypic data, imputation of missing values and preselection of markers
}
\description{
This function combines all algorithms for processing of marker data within \code{synbreed} package.
Raw marker data is a matrix with elements of arbitrary format (e.g. alleles coded as pair of observed alleles "A/T","G/C", ... , or by genotypes "AA", "BB", "AB"). The function is limited to biallelic markers with a maximum of 3 genotypes per locus. Raw data is recoded into the number of copies of the minor allele, i.e. 0, 1 and 2.
Imputation of missing values can be done by random sampling from allele distribution, the \code{Beagle} software or family information (see details). Additional preselection of markers can be carried out according to the minor allele frequency and/or fraction of missing values.
}
\usage{
codeGeno(gpData,impute=FALSE,impute.type=c("random","family",
         "beagle", "beagleAfterFamily","fix"), replace.value=NULL, maf=NULL,
         nmiss=NULL, label.heter="AB", keep.identical=TRUE, verbose=FALSE, 
         minFam=5, showBeagleOutput=FALSE, tester=NULL, print.report=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gpData}{
object of class \code{gpData} with arbitrary coding in element \code{geno}. Missing values have to be coded as \code{NA}.
}
  \item{impute}{
\code{logical}. Should missing value be replaced by imputing?
}
  \item{impute.type}{
\code{character} with one out of \code{"fix"}, \code{"random"} , \code{"family"}, \code{"beagle"}, \code{"beagleAfterFamily"} (default = \code{"random"}).  See details.
}
  \item{replace.value}{
\code{numeric} scalar to replace missing values in case \code{impute.type="fix"}.
}
  \item{maf}{
\code{numeric} scalar. Threshold to discard markers due to the minor allele frequency (MAF). Markers with a MAF < \code{maf} are discarded, thus  \code{maf} in [0,0.5]. If \code{map} in \code{gpData} is available, markers are also removed from \code{map}.
}
  \item{nmiss}{
\code{numeric} scalar.  Markers with more than \code{nmiss} fraction of missing values are discarded, thus  \code{nmiss} in [0,1]. If \code{map} in \code{gpData} is available, markers are also removed from \code{map}.
}
  \item{label.heter}{
This is either a scalar or vector of characters to identify heterozygous genotypes or a function returning \code{TRUE} if an element of the marker matrix is the heterozygous genotype. Defining a function is useful, if number of unique heterozygous genotypes is large, i.e. if genotypes are coded by alleles. If the heterozygous genotype is coded like "A/T","G/C", ... then \code{label.heter="alleleCoding"} can be used. Note that heterozygous values must be identified unambiguously by \code{label.heter}. Use \code{label.heter=NULL} if there are only homozygous genotypes, i.e. in DH lines, to speed up computation and restrict imputation to values 0 and 2.
}
  \item{keep.identical}{
\code{logical}. Should duplicated markers be kept? NOTE: From a set of identical markers (with respect to the non-missing alleles) the one with the smallest number of missing values is kept. For those with an identical number of missing values, the first one is kept and all others are removed.
}
 \item{verbose}{
\code{logical}. If \code{TRUE} verbose output is generated during the steps of the algorithm. This is useful to obtain numbers of discarded markers due to different criteria.
 }
 \item{minFam}{
For \code{impute.type} \code{family} and \code{beagleAfterFamily}, each family should have at least \code{minFam} members with available information for a marker to impute missing values according to the family. The default is 5.
 }
 \item{showBeagleOutput}{
\code{logical}. Would you like to see the output of the Beagle software package? The default is \code{FALSE}.
 }
 \item{tester}{
 This option is in testing mode at the moment.
 }
  \item{print.report}{
 \code{logical}. Should a file \code{SNPreport.txt} be generated containing further information on SNPs. This includes SNP name, original coding of major and minor allele, MAF and number of imputed values.
 }
}
\details{
Coding of genotypic data is done in the following order (depending on choice of arguments; not all steps are performed):

1. Discarding markers with fraction > \code{nmiss} of missing values

2. Recoding alleles from character/factor/numeric into the number of copies of the minor alleles, i.e. 0, 1 and 2. In \code{codeGeno}, in the first step heterozygous genotypes are coded as 1. From the other genotypes, the less frequent genotype is coded as 2 and the remaining genotype as 0. Note that  function \code{codeGeno} will terminate with an error whenever more than three genotypes are found.

2.1 Discarding duplicated markers if \code{keep.identical=FALSE} before starting of the imputing step. From identical marker based on pairwise complete oberservations one is discarded randomly. For getting identical results use the function \code{set.seed()} before \code{code.geno()}.

3. Replace missing values by \code{replace.value}  or impute missing values according to one of the following methods:

Imputing is done according to \code{impute.type}
\describe{
\item{"family"}{
This option is only suitable for homozygous individuals (such as doubled-haploid lines) structured in families.
Suppose an observation \eqn{i} is missing (NA) for a marker \eqn{j} in family \eqn{k}. If marker \eqn{j} is fixed in family \eqn{k}, the imputed value will be the fixed allele. If marker \eqn{j} is segregating for the population \eqn{k},
the value is 0 with probability of 0.5 and 2 with probability of 0.5. To use this algorithm, family information has to be stored as variable \code{family} in list element \code{covar} of an object of class \code{gpData}. This column should contain a \code{character} or \code{numeric} to identify family of all genotyped individuals.}
\item{"beagle"}{ Use Beagle Genetic Analysis Software Package (Browning and Browning 2007; 2009) to
    infer missing genotypes. If you use this option, please cite the original papers in publications. Beagle uses a HMM to reconstruct missing genotypes by the flanking
    markers. The beagle executive file \code{beagle.jar} (version 3.3.1) is in the directory \code{exec} of the package.
    Function \code{codeGeno} will create a directory \code{beagle} for Beagle
    input and output files (if it does not exist) and run Beagle with default settings. The information on marker
    position is taken from element \code{map}. Indeed, the postion in \code{map$pos} must be
    available for all markers. By default, three genotypes 0, 1, 2 are imputed. To restrict the
    imputation only to homozygous genotypes, use \code{label.heter=NULL}.}
\item{"beagleAfterFamily"}{
In the first step, missing genotypes are imputed according to the algorithm with \code{impute.type="family"}, but only for markers that are fixed within the family. Moreover, markers with a missing position (\code{map$pos=NA}) are imputed using the algorithm of \code{impute.type="family"}. In the second step, the remaining genotypes are imputed by Beagle.}
\item{"random"}{The missing values for a marker \eqn{j} are sampled from the marginal allele distribution of marker \eqn{j}. With 2 possible genotypes (to force this option, use \code{label.heter=NULL}), i.e. 0 and 2, values are sampled from distribution with probabilities \eqn{P(x=0)=1-p} and \eqn{P(x=2)=p}, where \eqn{p} is the minor allele frequency of marker \eqn{j}.  In the standardd case of 3 genotypes, i.e. with heterozygous genotypes, values are sampled from distribution \eqn{P(x=0)=(1-p)^2}, \eqn{P(x=1)=p(1-p)} and \eqn{P(x=2)=p^2} assuming Hardy-Weinberg equilibrium for all loci.}
\item{"fix"}{All missing values are imputed by \code{replace.value}. Note that only 0, 1 or 2 should be chosen.}
}

4. Recoding of alleles after imputation, if necessary due to changes in allele frequencies caused by the imputed alleles

5. Discarding markers with a minor allele frequency of <= \code{maf}

6. Discarding duplicated markers if \code{keep.identical=FALSE}. From identical marker based on pairwise complete oberservations one is discarded randomly. For getting identical results use the function \code{set.seed()} before \code{code.geno()}.

7. Restoring original data format (\code{gpData}, \code{matrix} or \code{data.frame})

Information about imputing is reported after a call of \code{codeGeno}.

Note: Beagle is included in the synbreed package. Once required, Beagle is called using \code{path.package()}. 
}
\value{
An object of class \code{gpData} containing the recoded marker matrix. If \code{maf} or \code{nmiss} were specified or \code{keep.identical=FALSE}, dimension of \code{geno} and \code{map} may be reduced due to selection of markers.  The genotype which is homozygous for the minor allele is coded as 2, the other homozygous genotype is coded as 0 and heterozygous genotype is coded as 1.
}

\references{
S R Browning and B L Browning (2007) Rapid and accurate haplotype phasing and missing data inference for whole genome association studies using localized haplotype clustering. Am J Hum Genet 81:1084-1097

B L Browning and S R Browning (2009) A unified approach to genotype imputation and haplotype phase inference for large data sets of trios and unrelated individuals. Am J Hum Genet 84:210-223
}

\author{
Valentin Wimmer and Hans-Juergen Auinger 
}


\examples{
# create marker data for 9 SNPs and 10 homozygous individuals
snp9 <- matrix(c(
  "AA",   "AA",   "AA",   "BB",   "AA",   "AA",   "AA",   "AA",  NA,
  "AA",   "AA",   "BB",   "BB",   "AA",   "AA",   "BB",   "AA",  NA,
  "AA",   "AA",   "AB",   "BB",   "AB",   "AA",   "AA",   "BB",  NA,
  "AA",   "AA",   "BB",   "BB",   "AA",   "AA",   "AA",   "AA",  NA,
  "AA",   "AA",   "BB",   "AB",   "AA",   "BB",   "BB",   "BB",  "AB",
  "AA",   "AA",   "BB",   "BB",   "AA",   NA,     "BB",   "AA",  NA,
  "AB",   "AA",   "BB",   "BB",   "BB",   "AA",   "BB",   "BB",  NA,
  "AA",   "AA",    NA,    "BB",    NA,    "AA",   "AA",   "AA",  "AA",
  "AA",    NA,     NA,    "BB",   "BB",   "BB",   "BB",   "BB",  "AA",
  "AA",    NA,    "AA",   "BB",   "BB",   "BB",   "AA",   "AA",  NA),
  ncol=9,byrow=TRUE)

# set names for markers and individuals
colnames(snp9) <- paste("SNP",1:9,sep="")
rownames(snp9) <- paste("ID",1:10+100,sep="")

# create object of class 'gpData'
gp <- create.gpData(geno=snp9)

# code genotypic data
gp.coded <- codeGeno(gp,impute=TRUE,impute.type="random")

# comparison
gp.coded$geno
gp$geno

# example with heterogeneous stock mice
\dontrun{
data(mice)
summary(mice)
# heterozygous values must be labeled  (may run some seconds)
mice.coded <- codeGeno(mice,label.heter=function(x) substr(x,1,1)!=substr(x,3,3))


# example with maize data and imputing by family
data(maize)
# first only recode alleles
maize.coded <- codeGeno(maize,label.heter=NULL)

# set 200 random chosen values to NA
set.seed(123)
ind1 <- sample(1:nrow(maize.coded $geno),200)
ind2 <- sample(1:ncol(maize.coded $geno),200)
original <- maize.coded$geno[cbind(ind1,ind2)]

maize.coded$geno[cbind(ind1,ind2)] <- NA
# imputing of missing values by family structure
maize.imputed <- codeGeno( maize.coded,impute=TRUE,impute.type="family",label.heter=NULL)


# compare in a cross table
imputed <- maize.imputed$geno[cbind(ind1,ind2)]
(t1 <- table(original,imputed) )
# sum of correct replacements
sum(diag(t1))/sum(t1)

# compare with random imputation
maize.random <- codeGeno(maize.coded,impute=TRUE,impute.type="random",label.heter=NULL)
imputed2 <- maize.random$geno[cbind(ind1,ind2)]
(t2 <- table(original,imputed2) )
# sum of correct replacements
sum(diag(t2))/sum(t2)
}
}
\keyword{manip}
