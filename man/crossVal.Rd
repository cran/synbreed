\name{crossVal}
\alias{crossVal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Cross validation of different prediction models
}
\description{Function for the application of the cross validation procedure on prediction
models with fixed and random effects. Covariance matrices must be committed to the function
and variance components can be committed or reestimated with ASReml or the BLR function.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
crossVal(gpData, trait=1, cov.matrix = NULL,  k = 2, Rep = 1, Seed = NULL,
         sampling = c("random", "within popStruc", "across popStruc","commit"),
         TS=NULL,ES=NULL, varComp = NULL, popStruc = NULL, VC.est = c("commit",
         "ASReml","BRR","BL"),verbose=FALSE,...)
}
\arguments{
  \item{gpData}{Object of class \code{gpData}
}
  \item{trait}{\code{numeric} or \code{character}. The name or number of the trait in the \code{gpData} object to be used as trait.
}
  \item{cov.matrix}{\code{list} including covariance matrices for the random effects. Size and order of rows and columns should be equal to rownames of \code{y}.
If no covariance is given, an identity matrix and marker genotypes are used for a marker regression. In general, a covariance matrix should be non-singular and positive definite to be invertible, if this is not the case, a constant of \code{1e-5} is added to the diagonal elements of the covariance matrix.
}
  \item{k}{\code{numeric}. Number of folds for k-fold cross validation, thus \code{k} should be in [2,\code{nrow(y)}] (default=2).
}
  \item{Rep}{\code{numeric}. Number of replications (default = 1).
}
  \item{Seed}{\code{numeric}. Number for \code{set.seed()} to make results reproducable.
}
  \item{sampling}{Different sampling strategies can be \code{"random"}, \code{"within popStruc"} or \code{"across popStruc"}. If sampling is \code{"commit"} test sets have to specified in TS (see Details).
}
  \item{TS}{A (optional) list of vectors with IDs for the test set in each fold within a list of replications, same layout as output for \code{id.TS} .
}
  \item{ES}{A (optional) list of IDs for the estimation set in each fold within each replication.
}
  \item{varComp}{A \code{vector} of variance components for the random effects, which has to be specified if \code{VC.est="commit"}. The first variance components should be
the same order as the given covariance matrices, the last given variance component is for the residuals.
}
  \item{popStruc}{Vector of length \code{nrow(y)} assigning individuals to a population structure. If no \code{popStruc} is defined, family information of \code{gpData} is used. Only required for options \code{sampling="within popStruc"} or  \code{sampling="across popStruc"}
}
  \item{VC.est}{Should variance components be reestimated with "\code{ASReml}" or with Bayesian Ridge Regression "\code{BRR}" or Bayesian Lasso "\code{BL}" of the \code{BLR} package within the estimation set of each fold in the cross validation? If
\code{VC.est="commit"}, the variance components have to be defined in \code{varComp}. For \code{ASReml}, ASReml software has to be installed on the system.
}
  \item{verbose}{Logical. Whether output shows replications and folds.
}
  \item{\dots}{
further arguments to be used by the genomic prediction models, i.e. prior values and MCMC options for the \code{BLR} function (see \code{\link[BGLR]{BLR}}).
}
}
\details{
In cross validation the data set is splitted into an estimation (ES) and a test set (TS). The effects are estimated with the ES and used to predict observations in the TS.
For sampling into ES and TS, k-fold cross validation is applied, where the data set is splitted into k subsets and k-1 comprising the ES and 1 is the TS, repeated for each subset.

To account for the family structure (Albrecht et al. 2011), \code{sampling} can be defined as:
\describe{
\item{random}{ Does not account for family structure, random sampling within the complete data set}
\item{within popStruc}{ Accounts for within population structure information, e.g. each family is splitted into k subsets}
\item{across popStruc}{ Accounts for across population structure information, e.g. ES and TS contains a set of complete families}
}
The following mixed model equation is used for \code{VC.est="commit"}:

\deqn{\bf y=\bf{Xb}+\bf{Zu}+\bf e}{y=Xb+Zu+e} with \deqn{\bf u \sim N(0,G\sigma^2_u)}{u=N(0,Gsigma2u)}
     gives the mixed model equations
\deqn{\left(\begin{array}{cc} \bf X'\bf X & \bf X'\bf Z \\ \bf Z'\bf X & \bf Z'\bf Z + \bf G^{-1}\frac{\sigma^2_e}{\sigma^2_u} \end{array} \right) \left( \begin{array}{c} \bf b \\ \bf  u \end{array}\right) = \left(\begin{array}{c}\bf  X'\bf y \\ \bf Z'\bf y \end{array} \right)}{(X'X,X'Z,Z'X,ZZ+GIsigma2e/sigma2u)(b,u)=(X'y,Z'y)}
}
\value{
An object of class \code{list} with following items:
\item{bu }{Estimated fixed and random effects of each fold within each replication.}
\item{n.DS}{Size of the data set (ES+TS) in each fold.}
\item{y.TS}{Predicted values of all test sets within each replication.}
\item{n.TS}{Size of the test set in each fold.}
\item{id.TS}{List of IDs of each test sets within a list of each replication.}
\item{PredAbi}{Predictive ability of each fold within each replication calculated as correlation coefficient \eqn{r(y_{TS},\hat y_{TS})}.}
\item{rankCor}{Spearman's rank correlation of each fold within each replication calculated between \eqn{y_{TS}} and \eqn{\hat y_{TS}}.}
\item{mse}{Mean squared error of each fold within each replication calculated between \eqn{y_{TS}} and \eqn{\hat y_{TS}}.}
\item{bias}{Regression coefficients of a regression of the observed values on the predicted values in the TS.
A regression coefficient \eqn{< 1} implies inflation of predicted values, and a coefficient of \eqn{> 1} deflation of predicted values.}
\item{m10}{Mean of observed values for the 10\% best predicted of each replication. The k test sets are pooled within each replication.}
\item{k}{Number of folds}
\item{Rep}{Replications}
\item{sampling}{Sampling method}
\item{Seed}{Seed for \code{set.seed()}}
\item{rep.seed}{Calculated seeds for each replication}
\item{nr.ranEff}{Number of random effects}
\item{VC.est.method}{Method for the variance components (\code{committed} or \code{reestimated with ASReml/BRR/BL})}
}
\references{
Albrecht T, Wimmer V, Auinger HJ, Erbe M, Knaak C, Ouzunova M, Simianer H, Schoen CC (2011) Genome-based prediction of testcross values in maize. Theor Appl Genet 123:339-350

Mosier CI (1951) I. Problems and design of cross-validation 1. Educ Psychol Measurement 11:5-11

Crossa J, de los Campos G, Perez P, Gianola D, Burgueno J, et al. (2010) Prediction of genetic values of quantitative traits in plant breeding using pedigree and molecular markers, Genetics 186:713-724

Gustavo de los Campos and Paulino Perez Rodriguez, (2010). BLR:
Bayesian Linear Regression. R package version 1.2.
http://CRAN.R-project.org/package=BLR

}
\author{Theresa Albrecht
}
\seealso{
\code{\link{summary.cvData}}}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
# loading the maize data set
\dontrun{
library(synbreedData)
data(maize)
maize2 <- codeGeno(maize)
U <- kin(maize2,ret="realized")
# cross validation
cv.maize  <- crossVal(maize2,cov.matrix=list(U),k=5,Rep=1,
            Seed=123,sampling="random",varComp=c(26.5282,48.5785),VC.est="commit")
cv.maize2 <- crossVal(maize2,k=5,Rep=1,
             Seed=123,sampling="random",varComp=c(0.0704447,48.5785),VC.est="commit")
# comparing results, both are equal!
cv.maize$PredAbi
cv.maize2$PredAbi
summary(cv.maize)
summary(cv.maize2)
}
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
