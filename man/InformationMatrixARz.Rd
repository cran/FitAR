\name{InformationMatrixARz}
\alias{InformationMatrixARz}
\title{ Fisher Information Matrix Subset Case, ARz }
\description{
Computes the large-sample Fisher information matrix per observation for the
AR coefficients in a subset AR when parameterized by the partial 
autocorrelations.
}
\usage{
InformationMatrixARz(zeta, lags)
}
\arguments{
  \item{zeta}{ vector of coefficients, ie. partial autocorrelations
  at lags specified in the argument \code{lags}  }
  \item{lags}{ lags in subset model, same length as zeta argument }
}
\details{
The details of the computation are given in 
McLeod and Zhang (2006, eqn 13).
\code{FitAR} uses \code{InformationMatrixARz} to obtain estimates
of the standard errors of the estimated parameters in the subset
AR model when partial autocorrelation parameterization is used.
}
\value{
a p-by-p Toeplitz matrix, p=length(zeta)
}
\references{ 
McLeod, A.I. and Zhang, Y. (2006).  
Partial autocorrelation parameterization for subset autoregression.  
Journal of Time Series Analysis, 27, 599-612.
}

\author{ A.I. McLeod and Y. Zhang}

\seealso{ 
\code{\link{FitAR}}, 
\code{\link{InformationMatrixAR}}, 
\code{\link{InformationMatrixARp}} 
}
\examples{
#Information matrix for ARz(1,4) with parameters 0.9 and 0.9.
InformationMatrixARz(c(0.9, 0.9), lags=c(1,4))
}
\keyword{ ts }
