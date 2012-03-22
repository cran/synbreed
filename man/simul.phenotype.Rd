\name{simul.phenotype}
\Rdversion{1.1}
\alias{simul.phenotype}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulation of a field trial with single trait
}
\description{
Simulates observations from a field trial using an animal model.
The field trial consists of multiple locations and randomized complete block design within locations. 
A single quantitative trait is simulated according to the model \code{Trait ~ id(A) + block + loc + e}.
}
\usage{
simul.phenotype(pedigree = NULL, A = NULL, mu = 100, vc = NULL,
                Nloc = 1, Nrepl = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pedigree}{
object of class "pedigree" 
}
  \item{A}{
object of class "relationshipMatrix"
}
  \item{mu}{
\code{numeric}; Overall mean of the trait.
}
  \item{vc}{
list containing the variance components. \code{vc} consists of elements \code{sigma2e}, \code{sigma2a}, \code{sigma2l}, \code{sigma2b} with the variance components of the residual, the additive genetic effect,
the location effect and the block effect.
}
  \item{Nloc}{
\code{numeric}. Number of locations in the field trial.
}
  \item{Nrepl}{
\code{Numeric}. Number of complete blocks within location.
}
}
\details{
Either \code{pedigree} or \code{A} must be specified. If \code{pedigree} is given, pedigree information is used to set up numerator relationship matrix with function \code{kinship}. If unrelated individuals should be used for simulation,
use identity matrix for \code{A}. True breeding values for \eqn{N} individuals is  simulated according to following distribution
\deqn{tbv \sim N(0,\bf A \sigma_a^2)}{tbv = chol(A)*sigma2a*rnorm(N,0,1)}       
Observations are simulated according to
\deqn{y \sim N(mu + tbv + block +loc,\sigma^2_e)}{}
If no location or block effects should appear, use \code{sigma2l=0} and/or \code{sigma2b=0}.
}
\value{
A \code{data.frame} with containing the simulated values for trait and the following variables
\item{ID}{Factor identifying the individuals. Names are extracted from \code{pedigree} or rownames of \code{A}}
\item{Loc}{Factor for Location}
\item{Block}{Factor for Block within Location}
\item{Trait}{Trait observations}
\item{TBV}{Simulated values for true breeding values of individuals}
Results are sorted for locations within individuals.
}
\author{
Valentin Wimmer
}
\seealso{
\code{\link{simul.pedigree}}
}
\examples{
\dontrun{
ped <- simul.pedigree(gener=5)
varcom <- list(sigma2e=25,sigma2a=36,sigma2l=9,sigma2b=4)
# field trial with 3 locations and 2 blocks within locations
data.simul <- simul.phenotype(ped,mu=10,vc=varcom,Nloc=3,Nrepl=2)
head(data.simul)
# analysis of variance
anova(lm(Trait~ID+Loc+Loc:Block,data=data.simul))
}
}

