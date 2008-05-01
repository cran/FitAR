\name{VarianceRacfARp}
\alias{VarianceRacfARp}
\title{ Covariance Matrix Residual Autocorrelations for ARp }
\description{
The ARp subset model is defined by taking a subset of the parameters
in the regular AR(p) model. With this function one can obtain the
standard deviations of the residual autocorrelations which can
be used for diagnostic checking with \code{\link{RacfPlot}}.
}
\usage{
VarianceRacfARp(phi, lags, MaxLag, n)
}
\arguments{
  \item{phi}{vector of AR coefficients  }
  \item{lags}{ lags in subset AR }
  \item{MaxLag}{covariance matrix for residual autocorrelations at 
  lags 1,...,m, where m=MaxLag is computes  }
  \item{n}{length of time series  }
}
\details{
The covariance matrix for the residual autocorrelations
is derived in McLeod (1978, eqn. 15) for the general ARMA case.
McLeod (1978, eqn. 35) specializes this result to the subset
case.
}
\value{
The m-by-m covariance matrix of residual autocorrelations at lags 
1,\dots,m, where m = MaxLag.
}
\references{
McLeod, A.I. (1978), 
On the distribution and applications of residual autocorrelations 
in Box-Jenkins modelling, 
\emph{Journal of the Royal Statistical Society B}, \bold{40}, 296-302.
  }
\author{ A.I. McLeod }

\seealso{ 
\code{\link{VarianceRacfAR}}, 
\code{\link{VarianceRacfARz}}, 
\code{\link{RacfPlot}} 
}
\examples{
#the standard deviations of the first 5 residual autocorrelations
#to a subset AR(1,2,6) model fitted to Series A is
v<-VarianceRacfARp(c(0.36,0.23,0.23),c(1,2,6), 5, 197)
sqrt(diag(v))
}
\keyword{ ts }

