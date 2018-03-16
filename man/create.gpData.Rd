\name{create.gpData}
\alias{create.gpData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create genomic prediction data object
}
\description{
This function combines all raw data sources in a single, unified data object of class \code{gpData}. This is a \code{list} with elements for phenotypic, genotypic, marker map, pedigree and further covariate data. All elements are optional.
}
\usage{
create.gpData(pheno = NULL, geno = NULL, map = NULL, pedigree = NULL,
              family = NULL, covar = NULL, reorderMap = TRUE, map.unit = "cM",
              repeated  = NULL, modCovar = NULL, na.string="NA", cores=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pheno}{
\code{data.frame} with individuals organized in rows and traits organized in columns. For
unrepeated measures unique \code{rownames} should identify individuals. For repeated measures, the first
column identifies individuals and a second column indicates repetitions (see also argument \code{repeated}).}
  \item{geno}{
\code{matrix} with individuals organized in rows and markers organized in columns. Genotypes could be coded arbitrarily. Missing values should be coded as \code{NA}. Colums or rows with only missing values not allowed. Unique \code{rownames} identify individuals and unique \code{colnames} markers. If no \code{rownames} are available, they are taken from element \code{pheno} (if available and if dimension matches).  If no \code{colnames} are used, the \code{rownames} of \code{map} are used if dimension matches.
}
  \item{map}{
\code{data.frame} with one row for each marker and two columns (named \code{chr} and \code{pos}). First columns gives the chromosome (\code{numeric} or \code{character} but not \code{factor}) and second column the position  on the chromosome in centimorgan or the physical distance relative to the reference sequence in basepairs. Unique \code{rownames} indicate the marker names which should match with marker names in \code{geno}. Note that order and number of markers must not be identical with the order in \code{geno}. If this is the case, gaps in the map are filled with \code{NA} to ensure the same number and order as in element \code{geno} of the resulting \code{gpData} object.
}
  \item{pedigree}{
Object of class \code{pedigree}.
}
  \item{family}{
\code{data.frame} assigning individuals to families with names of individuals in \code{rownames} This information could be used for replacing of missing values with function \code{codeGeno}.
}
  \item{covar}{
\code{data.frame} with further covariates for all individuals that either appear in \code{pheno}, \code{geno} or \code{pedigree$ID},  e.g. sex or age. \code{rownames} must be specified to identify individuals. Typically this element is not specified by the user.
}
  \item{reorderMap}{
\code{logical}. Should markers in \code{geno} and \code{map} be reordered by chromosome number and position within chromosome according to \code{map} (default = \code{TRUE})?
}
  \item{map.unit}{
\code{character}. Unit of position in \code{map}, i.e. 'cM' for genetic distance or 'bp' for physical distance (default = 'cM').
}
  \item{repeated}{
  This column is used to identify the replications of the phenotypic values. The unique values become the names of the third dimension of the pheno object in the \code{gpData}. This argument is only required for repeated measurements.
}
  \item{modCovar}{
  \code{vector} with \code{colnames} which identify columns with covariables in \code{pheno}. This argument is only required for repeated measurements.
}
  \item{na.string}{
\code{character} or vector of \code{characters}. You can specify values with which \code{NA} is coded in your geno object. In case you read missing values from a file not as missing, but as character strings. It can be specified more than one value for missings in a vector. Default is \code{"NA"}.
}
  \item{cores}{
 \code{numeric}. Here you can specify the number of cores you like to use.
}
}
\details{
The class \code{gpData} is designed to provide a unified framework for data related to genomic prediction analysis. Every data source can be omitted. In this case, the corresponding argument must be \code{NULL}. By default (argument \code{reorderMap}), markers in \code{geno} are ordered by their position in \code{map}. Individuals are ordered in alphabetical order.

An object of class \code{gpData} can contain different subsets of individuals or markers  in  the elements \code{pheno}, \code{geno} and \code{pedigree}. In this case the \code{id} in \code{covar} comprises all individuals that either appear in \code{pheno}, \code{geno} and \code{pedigree}. Two additional columns in \code{covar} named \code{phenotyped} and \code{genotyped} are automatically generated to identify individuals that appear in the corresponding \code{gpData} object.
}
\value{
Object of class \code{gpData} which is a \code{list} with the following elements
\item{covar}{\code{data.frame} with information on individuals}
\item{pheno}{\code{array} (individuals x traits x replications) with phenotypic data}
\item{geno}{\code{matrix} marker matrix containing genotypic data. Columns (marker) are in the same order as in \code{map} (if \code{reorderMap=TRUE}.) }
\item{pedigree}{object of class \code{pedigree}}
\item{map}{\code{data.frame} with columns 'chr' and 'pos' and markers sorted by 'pos' within 'chr'}
\item{phenoCovars}{\code{array} with phenotypic covariates}
\item{info}{\code{list} with additional information on data (coding of data, unit in \code{map})
 From synbreed version 0.11-11 on the function \code{\link{codeGeno}} adds here the package version which was used to do the coding. There are differences in codings between version 0.10-11 and 0.11-0!
}
}
\author{
Valentin Wimmer and Hans-Juergen Auinger with contributions be Peter VandeHaar
}
\note{
In case of missing row names or column names in one item, information is substituted from other elements (assuming the same order of individuals/markers) and a warning specifying the assumptions is returned. Please check them carefully.
}

\seealso{
\code{\link{codeGeno}}, \code{\link{summary.gpData}}, \code{\link{gpData2data.frame}}
}
\examples{
set.seed(123)
# 9 plants with 2 traits
n <- 9  # only for n > 6
pheno <- data.frame(Yield = rnorm(n,200,5), Height=rnorm(n,100,1))
rownames(pheno) <- letters[1:n]

# marker matrix
geno <- matrix(sample(c("AA","AB","BB",NA),size=n*12,replace=TRUE,
prob=c(0.6,0.2,0.1,0.1)),nrow=n)
rownames(geno) <-  letters[n:1]
colnames(geno) <- paste("M",1:12,sep="")

# genetic map
# one SNP is not mapped (M5) and will therefore be removed
map <- data.frame(chr=rep(1:3,each=4),pos=rep(1:12))
map  <- map[-5,]
rownames(map) <- paste("M",c(1:4,6:12),sep="")

# simulate pedigree
ped <- simul.pedigree(3,c(3,3,n-6))

# combine in one object
gp <- create.gpData(pheno,geno,map,ped)
summary(gp)


# 9 plants with 2 traits , 3 replcations
n <- 9  #
pheno <- data.frame(ID = rep(letters[1:n],3), rep = rep(1:3,each=n),
                    Yield = rnorm(3*n,200,5), Height=rnorm(3*n,100,1))

# combine in one object
gp2 <- create.gpData(pheno,geno,map,repeated="rep")
summary(gp2)
}
\keyword{manip}
