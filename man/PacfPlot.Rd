\name{PacfPlot}
\alias{PacfPlot}
\title{ Plot Partial Autocorrelations and Limits }
\description{
The sample partial autocorrelations and their individual 95 percent confidence intervals are plotted
under the assumption the model is contained in an AR(P), where P is a specified maximum order.
}
\usage{
PacfPlot(z, lag.max = 15, ...)
}
\arguments{
  \item{z}{ time series }
  \item{lag.max}{ maximum lag, P }
  \item{\dots}{ optional parameters passed through to plot. }
}
\details{
The Burg algorithm is used to estimate the PACF.
 
}
\value{
No value is returned. Graphical output is produced as side-effect.
}
\references{ 
McLeod, A.I. and Zhang, Y. (2006).  
Partial autocorrelation parameterization for subset autoregression.  
Journal of Time Series Analysis, 27, 599-612.
}
\author{ A.I. McLeod and Y. Zhang}

\seealso{ 
\code{\link{ar.burg}}
\code{\link{pacf}}
 }
\examples{
#For the log(lynx) series and taking lag.max=15, the PacfPlot and
# the minimum BIC subset selection produce the same result.
z<-log(lynx)
PacfPlot(z)
SelectModel(z,lag.max=15,ARModel="ARz",Best=1,Criterion="BIC")	
}
\keyword{ ts }
