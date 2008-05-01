\name{AcfPlot}
\alias{AcfPlot}
\title{Basic ACF Plotting }
\description{
Produces theoretical correlation plot
}
\usage{
AcfPlot(g, LagZeroQ= TRUE, ylab=NULL, main=NULL, ...)
}
\arguments{
  \item{g}{ vector of autocorrelations at lags 1,..,length(g) }
  \item{LagZeroQ}{ start plot at lag zero with g[0]=1 }
  \item{ylab}{ vertical axis label }
  \item{main}{ plot title }
  \item{\dots}{optional graphical parameters }
}
\value{
No value. Plot is produced via plot function.
}
\author{ A.I. McLeod and Y. Zhang}
\seealso{ \code{\link{acf}} }
\examples{
#
#Simple example, plot acf for AR(1)
 phi<-0.8
 maxLag<-20
 g<-phi^(1:maxLag)
 AcfPlot(g)
 AcfPlot(g, LagZeroQ=FALSE)
#
# Plot the sample inverse partial autocorrelations. 
# On the basis of this plot, Cleveland (1972) suggested an ARp(1,2,7)
# for this data
"InverseAcf" <-
function(z, p=15){
g<-TacvfMA(GetFitARpLS(z-mean(z),1:p)$phiHat, lag.max=p)
g/g[1]
}
#
data(SeriesA)
AcfPlot(InverseAcf(SeriesA),LagZeroQ=FALSE)
}
\keyword{ ts }

