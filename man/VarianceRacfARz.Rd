\name{VarianceRacfARz}
\alias{VarianceRacfARz}
\title{ Covariance Matrix Residual Autocorrelations for ARz }

\description{
The ARz subset model is defined by taking a subset of the 
partial autocorrelations (zeta parameters) in the AR(p) model.
With this function one can obtain the
standard deviations of the residual autocorrelations which can
be used for diagnostic checking with \code{\link{RacfPlot}}.
}

\usage{
VarianceRacfARz(zeta, lags, MaxLag, n)
}
\arguments{
  \item{zeta}{ zeta parameters (partial autocorrelations) }
  \item{lags}{ lags in model }
  \item{MaxLag}{covariance matrix for residual autocorrelations at 
  lags 1,...,m, where m=MaxLag is computes  }
  \item{n}{length of time series  }
}
\details{
The covariance matrix of the residual autocorrelations in the
subset ARz case is derived in McLeod and Zhang (2006, eqn. 16)
}

\value{
The m-by-m covariance matrix of residual autocorrelations at lags 
1,\dots,m, where m = MaxLag.
}

\references{ 
McLeod, A.I. and Zhang, Y. (2006).  
Partial autocorrelation parameterization for subset autoregression.  
\emph{Journal of Time Series Analysis}, \bold{27}, 599-612.
}
\author{ A.I. McLeod and Y. Zhang}

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


