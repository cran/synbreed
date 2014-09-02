\name{MME}
\alias{MME}
\title{
Mixed Model Equations
}
\description{
Set up Mixed Model Equations for given design matrices, i.e. variance components for random effects must be known.
}
\usage{
MME(X, Z, GI, RI, y)
}
\arguments{
  \item{X}{
Design matrix for fixed effects
}
  \item{Z}{
Design matrix for random effects
}
  \item{GI}{
Inverse of (estimated) variance-covariance matrix of random (genetic) effects multplied by the ratio of residual to genetic variance
}
  \item{RI}{
Inverse of (estimated) variance-covariance matrix of residuals (without multiplying with a constant, i.e. \eqn{\sigma^2_e}{sigma2e})
}
  \item{y}{
Vector of phenotypic records
}
}
\details{
The linear mixed model is given by
\deqn{\bf y = \bf X \bf b + \bf Z \bf u + \bf e }{y = Xb + Zu +e}
with \eqn{\bf u \sim N(0,\bf G)}{u = N(0,G)} and \eqn{\bf e \sim N(0,\bf R)}{e = N(0,R)}.
Solutions for fixed effects \eqn{b} and random effects \eqn{u} are obtained by solving the corresponding mixed model equations (Henderson, 1984)
\deqn{\left(\begin{array}{cc} \bf X'\bf R^{-1}\bf X & \bf X'\bf R^{-1}\bf Z \\ \bf Z'\bf R^{-1}\bf X & \bf Z'\bf R^{-1}\bf Z + \bf G^{-1} \end{array}\right) \left(\begin{array}{c} \bf  \hat b \\ \bf  \hat u \end{array}\right) = \left(\begin{array}{c}\bf  X'\bf  R^{-1} \bf y \\ \bf Z'\bf R^{-1}\bf y \end{array}\right)}{(X'RIX,X'RIZ,Z'RIX,ZRIZ+GI)(bhat,uhat)=(X'RIy,Z'RIy)}
Matrix on left hand side of mixed model equation is denoted by LHS and matrix on the right hand side of MME is denoted as RHS. Generalized Inverse of LHS equals prediction error variance matrix. Square root of diagonal values multiplied with \eqn{\sigma^2_e}{sigma2e} equals standard error of prediction. Note that variance components for fixed and random effects are not estimated by this function but have to be specified by the user, i.e. \eqn{G^{-1}}{GI} must be multiplied with shrinkage factor \eqn{\frac{\sigma^2_e}{\sigma^2_g}}{sigma2e/sigma2g}.
}
\value{
A list with the following arguments
\item{b}{Estimations for fixed effects vector}
\item{u}{Predictions for random effects vector}
\item{LHS}{left hand side of MME}
\item{RHS}{right hand side of MME}
\item{C}{Generalized inverse of LHS. This is the prediction error variance matrix}
\item{SEP}{Standard error of prediction for fixed and random effects}
\item{SST}{Sum of Squares Total}
\item{SSR}{Sum of Squares due to Regression}
\item{residuals}{Vector of residuals}
}
\references{
Henderson, C. R. 1984. Applications of Linear Models in Animal
Breeding. Univ. of Guelph, Guelph, ON, Canada.
}
\author{
Valentin Wimmer
}
\seealso{
\code{\link[regress]{regress}}, \code{\link{crossVal}}
}
\examples{
\dontrun{
library(synbreedData)
data(maize)

# realized kinship matrix
maizeC <- codeGeno(maize)
U <- kin(maizeC,ret="realized")/2

# solution with gpMod
m <- gpMod(maizeC,kin=U,model="BLUP")

# solution with MME
diag(U) <- diag(U) + 0.000001  # to avoid singularities
# determine shrinkage parameter
lambda <- m$fit$sigma[2]/ m$fit$sigma[1]
# multiply G with shrinkage parameter
GI <- solve(U)*lambda
y <- maizeC$pheno[,1,]
n <- length(y)
X <- matrix(1,ncol=1,nrow=n)
mme <- MME(y=y,Z=diag(n),GI=GI,X=X,RI=diag(n))

# comparison
head(m$fit$predicted[,1]-m$fit$beta)
head(mme$u)
}
}

