\name{VarianceRacfAR}
\alias{VarianceRacfAR}
\title{ Covariance Matrix Residual Autocorrelations for AR }
\description{
  Computes the variance-covariance matrix for the residual autocorrelations
in an AR(p).
}
\usage{
VarianceRacfAR(phi, MaxLag, n)
}
\arguments{
  \item{phi}{vector of AR coefficients  }
  \item{MaxLag}{covariance matrix for residual autocorrelations at 
  lags 1 ,..., m, where m = MaxLag is computes  }
  \item{n}{length of time series  }
}
\details{
The covariance matrix for the residual autocorrelations
is derived in McLeod (1978, eqn. 15) for the general ARMA case.
With this function one can obtain the
standard deviations of the residual autocorrelations which can
be used for diagnostic checking with \code{\link{RacfPlot}}.
}
\value{
The m-by-m covariance matrix of residual autocorrelations at lags 
1, ..., m, where m = MaxLag.
}
\references{
McLeod, A.I. (1978), 
On the distribution and applications of residual autocorrelations 
in Box-Jenkins modelling, 
\emph{Journal of the Royal Statistical Society B}, \bold{40}, 296--302
  }
\author{ A.I. McLeod }
\note{
The derivation assumes normality of the innovations,
mle estimation of the parameters and a known mean-zero time series.
It is easily seen that the same result still holds for
IID innovations with mean zero and finite variance,
any first-order efficient estimates of the parameters including
the AR coefficients and mean.
  }

\seealso{ 
\code{\link{VarianceRacfARp}}, 
\code{\link{VarianceRacfARz}}, 
\code{\link{RacfPlot}} 
}
\examples{
VarianceRacfAR(0.5,5,100)
}
\keyword{ ts }
