\name{plot.FitAR}
\alias{plot.FitAR}
\title{ Plot Method for "FitAR" Object }
\description{
Diagnostic plots: portmanteau p-values; residual autocorrelation plot;
normal probability plot and Jarque-Bera test; spectral density function
}
\usage{
\method{plot}{FitAR}(x, terse=TRUE, ...)
}
\arguments{
  \item{x}{ object of class "FitAR" }
  \item{terse}{ if TRUE, only one graph is produced, otherwise many diagnostic plots.}
  \item{...}{ optional arguments }
}
\value{
No value is returned. Plots are produced as side-effect.
}

\references{ 
McLeod, A.I. and Zhang, Y. (2006).  
Partial autocorrelation parameterization for subset autoregression.  
Journal of Time Series Analysis, 27, 599-612.
}

\author{ A.I. McLeod and Y. Zhang}

\note{When terse=FALSE, numerous graphs are produced. Turn on recording to be able to
page back and forth between them.}

\seealso{ 
\code{\link{summary.FitAR}}, 
\code{\link{FitAR}}, 
\code{\link{JarqueBeraTest}}
\code{\link{RacfPlot}}
\code{\link{LBQPlot}}
 }

\examples{
data(SeriesA)
obj<-FitAR(SeriesA, c(1,2,6,7))
plot(obj)
}
\keyword{ ts }



