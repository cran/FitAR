\name{coef.FitAR}
\alias{coef.FitAR}
\title{ Display Estimated Parameters from Output of "FitAR" }
\description{
Method function to display fitted parameters, their standard
errors and Z-ratio for AR models fit with FitAR.
}
\usage{
\method{coef}{FitAR}(object, ...)
}
\arguments{
  \item{object}{ \code{obj} the output from \code{FitAR} }
  \item{...}{optional parameters}
}
\value{
A matrix is returned. The columns of the matrix are labeled
MLE, sd and Z-ratio.  The rows labels indicate the AR coefficients
which were estimated followed by mu, the estimate of mean.
}
\references{ 
McLeod, A.I. and Zhang, Y. (2006).  
Partial autocorrelation parameterization for subset autoregression.  
Journal of Time Series Analysis, 27, 599-612.
}
\author{ A.I. McLeod and Y. Zhang}

\examples{
# Fit subset AR to SeriesA
 data(SeriesA) 
 outA<-FitAR(SeriesA, c(1,2,7), ARModel="ARz")
 coef(outA)
#
 outALS<-FitAR(SeriesA, c(1,2,7), ARModel="ARp")
 coef(outALS)

  }
}
\keyword{ ts }
