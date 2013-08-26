\name{predict.gpMod}
\alias{predict.gpMod}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prediction for genomic prediction models.
}
\description{
S3 \code{predict} method for objects of class \code{gpMod}. A genomic prediction model is used 
to predict the genetic performance for e.g. unphenotyped individuals using an object of class 
\code{gpMod} estimated by a training set.
}
\usage{
\method{predict}{gpMod}(object,newdata,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
object of class \code{gpMod} which is the model used for the prediction. If the model includes a \code{relationshipMatrix}, this must include both the individuals in the training data used for fitting \code{gpMod} and those which sould be predicted in \code{newdata} (see example below).
}
  \item{newdata}{
for \code{model="BL"} and \code{"BRR"} an object of class \code{gpData} with the marker data of 
the unphenotyped individuals. For \code{model="BLUP"} a \code{character} vector with the names 
of the individuals to predict. If \code{newdata=NULL},  the genetic performances of the 
individuals for the training set are returned.
}
\item{...}{not used
}
}
\details{
For models, \code{model="RR"} and \code{"BL"}, the prediction for the unphenotyped individuals
is given by \deqn{\hat{g}'=\hat{\mu} + W'\hat{m}}{ghat=muhat + W mhat} with the estimates taken 
from the \code{gpMod} object. For the prediction using \code{model="BLUP"}, the full 
relationship matrix including individuals of the training set and the prediction set must be 
specified in the \code{gpMod}. This model is used to predict the unphenotyped individuals of 
the prediction set by solving the corresponding mixed model equations using the variance 
components of the fit in \code{gpMod}.
}
\value{
a named vector with the predicted genetic values for all individuals in \code{newdata}.
}
\references{
Henderson C (1977) Best linear unbiased prediction of breeding values not in the model for 
records. Journal of Dairy Science 60:783-787

Henderson CR (1984). Applications of linear models in animal breeding. University of Guelph.
}
\author{
Valentin Wimmer
}

\seealso{
\code{\link{gpMod}}
}
\examples{
# Example from Henderson (1977)
dat <- data.frame(y=c(132,147,156,172),time=c(1,2,1,2),animal=c(1,2,3,4))
ped <- create.pedigree(ID=c(6,5,1,2,3,4),Par1=c(0,0,5,5,1,6),Par2=c(0,0,0,0,6,2))
gp <- create.gpData(pheno=dat,pedigree=ped)
A <- kin(gp,ret="add")

# assuming h2=sigma2u/(sigma2u+sigma2)=0.5
# no REML fit possible due to the limited number of observations
y <- c(132,147,156,172)
names(y) <- paste(1:4)
mod1 <- list(fit=list(sigma=c(1,1)),kin=A,model="BLUP",y=y,m=NULL)
# matrix A included all individuals (including those which should be predicted)
class(mod1) <- "gpMod"
predict(mod1,c("5","6"))

# prediction 'by hand'
X <- matrix(1,ncol=1,nrow=4)
Z <- diag(6)[-c(1,2),]
AI <- solve(A)
RI <- diag(4)

res <- MME(X,Z,AI,RI,y)
res$b + res$u[1:2]

%# prediction of genetic performance of the last 50 individuals in the maize data set
%data(maize)
%maizeC <- codeGeno(maize)
%U <- kin(maizeC,ret="realized")
%maizeC2 <- discard.individuals(maizeC,rownames(maizeC$pheno)[1201:1250])
%modU <- gpMod(maizeC2,model="BLUP",kin=U)
%predict(modU,rownames(maizeC$pheno)[1201:1250])
}
