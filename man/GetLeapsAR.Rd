\name{GetLeapsAR}
\alias{GetLeapsAR}

\title{ Select lags for Best Subset ARp Model }

\description{
The subset ARp model is the usual subset model, for example see Tong (1977).
This function is used by SelectModel for model identification for
ARp models.
}

\usage{
GetLeapsAR(z, lag.max = 15, Criterion = "UBIC", Best = 3, Candidates=5, t="default", ExactQ=FALSE)
}

\arguments{
  \item{z}{ts object or vector containing time series }
  \item{lag.max}{maximum order of the AR }
  \item{Criterion}{default UBIC, other choices are "AIC", "BIC", "EBIC", "BICq", "GIC" }
  \item{Best}{the number of based selected. Ignore with "GIC". }
  \item{Candidates}{number of models initially selected using the approximate criterion }
  \item{t}{tuning parameter, EBIC, BICq, GIC  }
  \item{ExactQ}{exhaustive numeration using exact likelihood. Still under
    under development. NOT AVAILABLE IN THIS VERSION }
}

\details{
The \R function leaps in the \R package leaps is used to compute the subset
regression model with the smallest residual sum of squares containing
1, \dots, lag.max parameters.  The mean is always included, so the only
parameters considered are the phi coefficients. After the best models
containing 1, \dots, lag.max parameters are selected the models are individually
refit to determine the exact likelihood function
for each selected model. Based on this likelihood the UBIC/BIC/AIC is
computed and then the best models are selected. The UBIC criterion was
developed by Chen and Chen (2007).
  The EBIC using a tuning parameter, G, where 0 <= G <= 1.
  The BICq takes a tuning parameter, Q, where 0 < Q < 1.
  The GIC takes a tuning parameter, p, where 0<p<0.25.
}

\value{
When 'Criterion' is one of UBIC, AIC, BIC, EBIC, BICq, a list with components:
  \item{p}{lags present in model}
  \item{UBIC}{approximate UBIC (Chen & Chen, 2007), if Criterion=="UBIC"}
  \item{AIC}{approximate AIC (McLeod and Zhang, 2006a, eqn. 15), if Criterion=="AIC"}
  \item{BIC}{approximate BIC (McLeod and Zhang, 2006a, eqn. 15), if Criterion=="BIC"}
  \item{EBIC}{approximate EBIC (McLeod and Zhang, 2006a, eqn. 15), if Criterion=="EBIC"}
  \item{BICq}{approximate BICq, if Criterion=="BICq"}
  \item{GIC}{approximate GIC, if Criterion=="GIC"}

}

\references{ 
Tong, H. (1977)
Some comments on the Canadian lynx data.
Journal of the Royal Statistical Society A 140, 432-436.

Chen, J. and Chen, Z. (2008).
Extended Bayesian Information Criteria for Model Selection with Large Model Space.
Biometrika.

Changjiang Xu and A. I. McLeod (2010).
Bayesian information criterion with Bernoulli prior.
Submitted for publication.

Changjiang Xu and A. I. McLeod (2010).
Model selection using generalized information criterion.
Submitted for publication.

}

\author{ A.I. McLeod }

\note{Requires leaps package. 
Since the least-squares is used, the number of observations
depends on 'lag.max'. Hence different subsets may be chosen
depending on the 'lag.max. See example below.}

\section{Warning }{
AIC and BIC values produced are not comparable to AIC and BIC produced
by SelectModel for ARz models. 
However comparable AIC/BIC values are produced when the selected models
are fit by FitAR.
}

\seealso{
 \code{\link{SelectModel}}, 
 \code{\link{GetFitARpLS}},
 \code{\link[leaps]{leaps}} 
}

\examples{
#Example 1: Simple Example
#for the log(lynx) Tong (1977) selected an ARp(1,2,4,10,11)
#using the AIC and a subset selection algorithm. Our more exact
#approach shows that the ARp(1,2,3,4,10,11) has slightly lower
#AIC (using exact likelihood evaluation).  
z<-log(lynx)
GetLeapsAR(z, lag.max=11)
GetLeapsAR(z, lag.max=11, Criterion="BIC")

#Example 2: Subset autoregression depends on lag.max!
#Because least-squares is used, P=lag.max observations are
#  are deleted. This causes different results depending on lag.max.
#This phenomenon does not happen with "ARz" subset models
#ARp models depend on lag.max
GetLeapsAR(z, lag.max=15, Criterion="BIC")
GetLeapsAR(z, lag.max=20, Criterion="BIC")

#Example 3: Comparing GIC with BIC, AIC, UBIC and BICq
z <- log(lynx)
GetLeapsAR(z, lag.max=15, Criterion="BIC", Best=1)
GetLeapsAR(z, lag.max=15, Criterion="AIC", Best=1)
GetLeapsAR(z, lag.max=15, Criterion="UBIC", Best=1)
GetLeapsAR(z, lag.max=15, Criterion="BICq", Best=1, t=0.25)
GetLeapsAR(z, lag.max=15, Best=1, Criterion="GIC", t=0.01)
ans<-GetLeapsAR(z, lag.max=15, Best=3, Criterion="GIC", t=0.001)
plot(ans)
}

\keyword{ ts }
