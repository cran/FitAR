\name{Boot.ts}
\alias{Boot.ts}
\title{ Parametric Time Series Bootstrap }
\description{
  An AR(p) model is fit to the time series using the AIC and then
  it is simulated.
}
\usage{
\method{Boot}{ts}(obj, R=1, ...)
}
\arguments{
  \item{obj}{a time series, class "ts" }
  \item{R}{number of bootstrap replicates }
  \item{\dots}{optional arguments }
}
\value{
A time series or vector.
}
\references{
Davison, A.C. and Hinkley, D.V. (1997), Bootstrap Methods and
Their Application. Cambridge University Press.
}
\author{ A.I. McLeod and Y. Zhang}
\note{
Parametric and nonparametric time series bootstraps are discussed
by Davison and Hinkley (1997, Ch.8.2).
}

\seealso{ \code{\link{Boot.FitAR}}.
Nonparametric bootstrap for time series is available
in the function \code{tsboot} in the library \code{boot}.
 }
\examples{
layout(matrix(c(1,2,1,2),ncol=2))
TimeSeriesPlot(SeriesA)
TimeSeriesPlot(Boot(SeriesA),main="Bootstrap of Series A")
}
\keyword{ ts }
