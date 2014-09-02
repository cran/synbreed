\name{gpMod}
\alias{gpMod}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Genomic predictions models for objects of class gpData
}
\description{
This function fits genomic prediction models based on phenotypic and genotypic data in an object of class \code{gpData}. The possible models are Best Linear Unbiased Prediction (BLUP) using a pedigree-based or a marker-based genetic relationship matrix and Bayesian Lasso (BL) or Bayesian Ridge regression (BRR).
BLUP models are fitted using the REML implementation of the \code{regress} package (Clifford and McCullagh, 2012). The Bayesian regression models are fitted using the Gibbs-Sampler of the \code{BLR} package (de los Campos and Perez, 2010). The covariance structure in the BLUP model is defined by an object of class \code{relationshipMatrix}. The training set for the model fit consists of all individuals with phenotypes and genotypes. All data is restricted to individuals from the training set used to fit the model.


} \usage{
gpMod(gpData, model=c("BLUP","BL","BRR"), kin=NULL, predict=FALSE, trait=1,
      repl=NULL, markerEffects=FALSE, fixed=NULL, random=NULL, ...) }
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gpData}{
object of class \code{gpData}
}
  \item{model}{
\code{character}. Type of genomic prediction model. \code{"BLUP"} indicates best linear unbiased prediction (BLUP) using REML for both pedigree-based (P-BLUP) and marker-based (G-BLUP) model. \code{"BL"} and \code{"BRR"} indicate Bayesian Lasso and Bayesian Ridge Regression, respectively.
}
  \item{kin}{
object of class \code{relationshipMatrix} (only required for \code{model = "BLUP"}). Use a pedigree-based kinship to evaluate P-BLUP or a marker-based kinship to evaluate G-BLUP. For \code{"BL"} and \code{"BRR"}, also a kinship structure may be used as additional polygenic effect \eqn{u} in the Bayesian regression models (see \code{BLR} package).
}
  \item{predict}{
\code{logical}. If \code{TRUE}, genetic values will be predicted for genotyped but not phenotyped individuals. Default is \code{FALSE}. Note that this option is only meaningful for marker-based models. For pedigree-based model, please use function \code{predict.gpMod}.
}
  \item{trait}{
\code{numeric} or \code{character}. A vector with names or numbers of the traits to fit the model
}
  \item{repl}{
\code{numeric} or \code{character}. A vector with names or numbers of the repeated values of \code{gpData$pheno} to fit the model
}
  \item{markerEffects}{
  \code{logical}.
 Should marker effects be estimated for a G-BLUP model, i.e. RR-BLUP? In this case, argument \code{kin} is ignored (see Details). Plose note, that in this case also the variance components pertaining to model G-BLUP are reported instead of those from the G-BLUP model (see vignette). If the variance components are committed to \code{crossVal}, it must be guaranteed that there also the RR-BLUP model is used, e.g. no \code{cov.matrix} object should be specified.
}
  \item{fixed}{
  A formula for fixed effects. The details of model specification are the same as for \code{lm} (only right hand side required). Only for \code{model="BLUP"}.
}
  \item{random}{
  A formula for random effects of the model. Specifies the matrices to include in the covariance structure. Each term is either a symmetric matrix, or a factor. Independent Gaussian random effects are included by passing the corresponding block factor. For mor details see \code{\link[regress]{regress}}. Only for \code{model="BLUP"}
}
  \item{\dots}{
further arguments to be used by the genomic prediction models, i.e. prior values and MCMC options for the \code{BLR} function (see \code{\link[BLR]{BLR}}) or parameters for the REML algorithm in regress.
}
}
\details{
By default, an overall mean is added to the model. If no \code{kin} is specified and  \code{model = "BLUP"}, a G-BLUP model will be fitted. For BLUP, further fixed and random effects can be added through the arguments \code{fixed} and \code{random}.

The marker effects \eqn{\hat{m}}{hatm} in the RR-BLUP model (available with \code{markerEffects}) are calculated as \deqn{\hat{m}= X'G^{-1}\hat{g}}{m = X'*Ginv*ghat} with \eqn{X} being the marker matrix, \eqn{G=XX'} and \eqn{hat{g}}{ghat} the vector of predicted genetic values.

Only a subset of the individuals - the training set - is used to fit the model. This contains all individuals with phenotypes and genotypes. If \code{kin} does not match the dimension of the training set (if, e.g. ancestors are included), the respective rows and columns from the trainings set are choosen.
}
\value{
Object of class \code{gpMod} which is a list of
\item{fit}{The model fit returned by the genomic prediction method}
\item{model}{The model type, see 'Arguments'}
\item{y}{The phenotypic records for the individuals in the training set}
\item{g}{The predicted genetic values for the individuals in the training set}
\item{m}{Predicted SNP effects (if available)}
\item{kin}{Matrix \code{kin}}
}
\note{
The verbose output of the \code{BLR} function is written to a file \code{BLRout.txt} in the working directory to prevent the screen output from overload.
}
\references{
Clifford D, McCullagh P (2012). regress: Gaussian Linear Models with Linear Covariance
Structure. R package version 1.3-8, URL http://www.csiro.au.

Gustavo de los Campos and Paulino Perez Rodriguez, (2010). BLR:
Bayesian Linear Regression. R package version 1.2.
http://CRAN.R-project.org/package=BLR

}
\author{
Valentin Wimmer, Hans-Juergen Auinger and Theresa Albrecht
}
\seealso{
\code{\link{kin}}, \code{\link{crossVal}}
}
\examples{
\dontrun{
library(synbreedData)
data(maize)
maizeC <- codeGeno(maize)

# pedigree-based (expected) kinship matrix
K <- kin(maizeC,ret="kin",DH=maize$covar$DH)

# marker-based (realized) relationship matrix
# divide by an additional factor 2
# because for testcross prediction the kinship of DH lines is used
U <- kin(maizeC,ret="realized")/2
# BLUP models
# P-BLUP
mod1 <- gpMod(maizeC,model="BLUP",kin=K)
# G-BLUP
mod2 <- gpMod(maizeC,model="BLUP",kin=U)

# Bayesian Lasso
prior <- list(varE=list(df=3,S=35),lambda = list(shape=0.52,rate=1e-4,value=20,type='random'))
mod3 <- gpMod(maizeC,model="BL",prior=prior,nIter=6000,burnIn=1000,thin=5)

summary(mod1)
summary(mod2)
summary(mod3)}
}
