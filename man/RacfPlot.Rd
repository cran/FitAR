\name{RacfPlot}
\alias{RacfPlot}
\title{ Residual Autocorrelation Plot }
\description{
Residual autocorrelation plot for "FitAR" objects.  This
plot is useful for diagnostic checking models fit with the
function \code{\link{FitAR}}.
}
\usage{
RacfPlot(obj, lag.max = 1000, SquaredQ=FALSE, ylab="")
}

\details{
The standard deviations of the residual autocorrelations are
obtained from McLeod (1978, eqn.16) or
McLeod and Zhang (2006, eqn.16).
Simultaneous confidence bounds are shown and constructed
using the Bonferonni approximation as suggested by
Hosking and Ravishanker (1993)
}

\arguments{
  \item{obj}{ output from FitAR }
  \item{lag.max}{ maximum lag. Set to 1000 since minimum of this value and the value
  in the \code{obj} is used. }
  \item{SquaredQ}{ default is FALSE. For squared residual autocorrelations, set to TRUE}
  \item{ylab}{y-axis label}
}
\value{
Plot is produced as a side-effect. No output
}
\references{
Hosking, J.R.M. and Ravishanker, N. (1993)
Approximate simultaneous significance intervals for residual autocorrelations of
autoregressive-moving average time series models.
Journal of Time Series Analysis 14, 19-26.

McLeod, A.I. (1978), 
On the distribution and applications of residual autocorrelations 
in Box-Jenkins modelling, 
Journal of the Royal Statistical Society B 40, 296-302.

McLeod, A.I. and Zhang, Y. (2006).  
Partial autocorrelation parameterization for subset autoregression.  
Journal of Time Series Analysis, 27, 599-612.
}
\author{A.I. McLeod and Y. Zhang }
\note{ This function is normally invoked when plot.FitAR
is used. }


\seealso{
 \code{\link{plot.FitAR}}, 
 \code{\link{FitAR}}, 
}
\examples{
#fit subset AR and plot diagnostic check
 data(SeriesA)
 out<-FitAR(SeriesA, c(1,2,7), ARModel="ARp")
 RacfPlot(out)
#note that plot produces LBQPlot and RacfPlot
 plot(out)
#check squared residuals
 RacfPlot(out, SquaredQ=TRUE)

}
\keyword{ ts }






