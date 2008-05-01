\name{GetLeapsAR}
\alias{GetLeapsAR}

\title{ Select lags for Best Subset ARp Model }

\description{
The subset ARp model is the usual subset model, for example see Tong (1977).
This function is used by SelectModel for model identification for
ARp models.
}

\usage{
GetLeapsAR(z, lag.max = 15, Criterion = "UBIC", Best = 3, Candidates=5)
}

\arguments{
  \item{z}{ ts object or vector containing time series }
  \item{lag.max}{ maximum order of the AR }
  \item{Criterion}{ default UBIC, other choices are "AIC" or "BIC" }
  \item{Best}{ the number of based selected based on highest AIC/BIC }
  \item{Candidates}{ number of models initially selected using the approximate criterion }
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
}

\value{
a list with components
  \item{NumParameters} {Number of parameters in model}
  \item{UBIC} {approximate UBIC (Chen & Chen, 2007), if Criterion=="UBIC"}
  \item{AIC} {approximate AIC (McLeod and Zhang, 2006a, eqn. 15), if Criterion=="AIC"}
  \item{BIC} {approximate BIC (McLeod and Zhang, 2006a, eqn. 15), if Criterion=="BIC"}
  \item{p} {lags present}
}

\references{ 
Tong, H. (1977)
Some comments on the Canadian lynx data.
Journal of the Royal Statistical Society A 140, 432-436.


Chen, J. and Chen, Z. (2007).
Extended Bayesian Information Criteria for Model Selection with Large Model Space.
Preprint.

 
}

\author{ A.I. McLeod }

\note{Requires leaps package}

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
#for the log(lynx) Tong (1977) selected an ARp(1,2,4,10,11)
#using the AIC and a subset selection algorithm. Our more exact
#approach shows that the ARp(1,2,3,4,10,11) has slightly lower
#AIC (using exact likelihood evaluation).  
z<-log(lynx)
GetLeapsAR(z, lag.max=11)
GetLeapsAR(z, lag.max=11, Criterion="BIC")
}

\keyword{ ts }
