\name{PlotARSdf}
\alias{PlotARSdf}
\title{ Plot AR or ARMA Spectral Density }
\description{
Constructs a plot of the AR spectral density function.
}
\usage{
PlotARSdf(phi = NULL, theta = NULL, units = "radial", logSdf = FALSE, InnovationVariance = 1, main = NULL, sub = NULL, lwd=3, col="blue", ...)
}
\arguments{
  \item{phi}{ AR Coefficients}
  \item{theta}{ MA Coefficients}
  \item{units}{ default is "radial" }
  \item{logSdf}{ default is FALSE otherwise log sdf is plotted }
  \item{InnovationVariance}{ innovation variance, default is 1 }
  \item{main}{ optional plot title }
  \item{sub}{ optional subtitle }
  \item{lwd}{ optional lwd for plot, default lwd=3. }
  \item{col}{ optional col for plot. Default "blue". }
  \item{...}{ optional arguments }
}
\details{
The spectral density function is symmetric and defined in
(-pi, pi) but plotted over (0, pi).  
If units are not "radial", it is plotted over (0, 0.5).
}
\value{
Plot is produced using plot.
Matrix with 2 columns containing the frequencies and spectral density is returned invisibly.
}

\seealso{ \code{\link{ARSdf}} }
\examples{
#AR(1)
PlotARSdf(0.8)
#MA(1)
PlotARSdf(theta=0.8)
#ARMA(1,1)
PlotARSdf(0.9,0.5)
#white noise
PlotARSdf()
}
\keyword{ ts }
