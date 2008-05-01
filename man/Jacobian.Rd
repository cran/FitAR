\name{Jacobian}
\alias{Jacobian}
\title{ Jacobian AR-coefficients to Partial Autocorrelations }
\description{
This is more or less and internal routine used by InformationMatrixZeta
but it is described in more details since it may be useful in other
computations.
}
\usage{
Jacobian(zeta)
}
\arguments{
  \item{zeta}{ partial autocorrelation parameters }
}
\details{
The computation is described in detail in 
McLeod and Zhang (2006, Section 2.2)
}
\value{
square matrix of order length(zeta)
}

\references{ 
McLeod, A.I. and Zhang, Y. (2006).  
Partial autocorrelation parameterization for subset autoregression.  
Journal of Time Series Analysis, 27, 599-612.
}

\author{ A.I. McLeod }

\seealso{ \code{\link{InformationMatrixARz}} }
\examples{
#In McLeod and Zhang (2006, p.603) a symbolic example is given
# for the AR(4).  
#
 Jacobian(rep(0.8,4))

}
\keyword{ ts }
