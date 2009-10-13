\name{SimulateGaussianAR}
\alias{SimulateGaussianAR}
\title{ Autoregression Simulation }
\description{
  Simulate a mean-zero stationary Gaussian AR(p) time series.
}
\usage{
SimulateGaussianAR(phi, n = 100, InnovationVariance = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phi}{ vector containing AR coefficients }
  \item{n}{ length of time series }
  \item{InnovationVariance}{ innovation variance }
}
\details{
The p initial values are simulated using the appropriate multivariate distribution
as was suggested in McLeod (1975).  The R function rnorm() is used.
}
\value{
A vector of length n, the simulated series
}
\references{
McLeod, A.I. (1975), Derivation of the theoretical autocorrelation function of 
autoregressive moving-average time series, \emph{Applied Statistics} \bold{24}, 255--256.  
Percival, D.B. and Walden, A.T. (1993), \emph{Spectral Analysis for Physical Applications}.
}

\author{ A.I. McLeod }
\note{If the process is non-stationary, then random initial values are
used determined by the first p innovations.}

\seealso{ \code{\link{Boot.FitAR}} }
\examples{
#Percival and Walden (1993, p.46) illustrated a time series with a 
#very peaked spectrum with the AR(4) with coefficients
#c(2.7607,-3.8106,2.6535,-0.9238) with NID(0,1) innovations.
#
z<-SimulateGaussianAR(c(2.7607,-3.8106,2.6535,-0.9238),1000)
library(lattice)
TimeSeriesPlot(z)
}
\keyword{ ts }
