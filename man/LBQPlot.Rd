\name{LBQPlot}
\alias{LBQPlot}
\title{ Plot Ljung-Box Test P-value vs Lag }
\description{
The Ljung-Box portmanteau p-value is plotted vs lag.
}
\usage{
LBQPlot(res, lag.max = 30, StartLag = k + 1, k = 0, SquaredQ = FALSE)
}
\arguments{
  \item{res}{ residuals }
  \item{lag.max}{ maximum lag }
  \item{StartLag}{ starting lag }
  \item{k}{ number of AR parameters fit }
  \item{SquaredQ}{ default, SquaredQ = FALSE, regular autocorrelations. If SquaredQ = TRUE use
  autocorrelations of squared residuals. }
}
\value{
Plot is produced as a side-effect. No output
}
\references{ 
Ljung, G.M. and Box, G.E.P. (1978)
On a measure of lack of fit in time series models.
Biometrika 65, 297-303.

McLeod, A.I. and Zhang, Y. (2006).  
Partial autocorrelation parameterization for subset autoregression.  
Journal of Time Series Analysis, 27, 599-612.

}
\author{A.I. McLeod  and Y. Zhang}
\note{ This function is normally invoked when plot.FitAR
is used. }


\seealso{
 \code{\link{plot.FitAR}}, 
 \code{\link{FitAR}}
}
\examples{
#fit subset AR and plot diagnostic check
 out<-FitAR(SeriesA, c(1,2,7), ARModel="ARp")
 res<-resid(out)
 LBQPlot(res)
#note that plot produces LBQPlot and RacfPlot
 plot(out)
}
\keyword{ ts }
