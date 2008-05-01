\name{print.FitAR}
\alias{print.FitAR}
\title{ Print Method for "FitAR" Object }
\description{
A terse summary is given.
}
\usage{
\method{print}{FitAR}(x, ...)
}
\arguments{
  \item{x}{ object of class "FitAR" }
  \item{...}{ optional arguments }
}
\value{
A terse summary is displayed
}
\references{ 
McLeod, A.I. and Zhang, Y. (2006).  
Partial autocorrelation parameterization for subset autoregression.  
Journal of Time Series Analysis, 27, 599-612.
}
\author{ A.I. McLeod and Y. Zhang}

\seealso{ 
\code{\link{summary.FitAR}}
 }
\examples{
data(SeriesA)
FitAR(SeriesA, c(1,2,6,7))
}
\keyword{ ts }
