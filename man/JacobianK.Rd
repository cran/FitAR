\name{JacobianK}
\alias{JacobianK}
\title{ Internal Utility Function }
\description{
The matrix defined in eqn. (10) of McLeod and Zhang (2006).
Used by the function \code{\link{Jacobian}}.
}
\usage{
JacobianK(zeta, k)
}
\arguments{
  \item{zeta}{ partial autocorrelations }
  \item{k}{ k-th Jacobian }
}
\value{
Matrix
}
\references{ 
McLeod, A.I. and Zhang, Y. (2006).  
Partial autocorrelation parameterization for subset autoregression.  
Journal of Time Series Analysis, 27, 599-612.
}

\author{ A.I. McLeod }

\seealso{ \code{\link{Jacobian}} }
\examples{
JacobianK(rep(0.8,4),3)
}
\keyword{ ts }

