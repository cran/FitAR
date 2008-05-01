\name{DetAR}
\alias{DetAR}
\title{ Covariance Determinant of AR(p) }
\description{
Computes the covariance determinant of p successive observations
from an AR(p) process with unit innovation variance.
}
\usage{
DetAR(phi)
}
\arguments{
  \item{phi}{vector of AR coefficients}
}
\details{
The AR coefficients are transformed to PACF and then the
determinant is computed as a product of PACF terms as given
in McLeod and Zhang (2006, eqn. 4).
}
\value{
Determinant
}
\references{ 
McLeod, A.I. and Zhang, Y. (2006).  
Partial autocorrelation parameterization for subset autoregression.  
Journal of Time Series Analysis, 27, 599-612.
}
\author{ A.I. McLeod and Y. zhang}

\seealso{ \code{\link{FastLoglikelihoodAR}} }
\examples{
DetAR(c(0.1,0.1,0.1))
}
\keyword{ ts }
