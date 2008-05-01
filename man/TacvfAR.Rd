\name{TacvfAR}
\alias{TacvfAR}
\title{ Theoretical Autocovariance Function of AR }
\description{
The theoretical autocovariance function of an AR(p) with unit
variance is computed. This algorithm has many applications.
In this package it is used for the computation of the information
matrix, in simulating p initial starting values for AR
simulations and in the computation of the exact mle for the mean.
}
\usage{
TacvfAR(phi, lag.max = 20)
}
\arguments{
  \item{phi}{ vector of AR coefficients }
  \item{lag.max}{ computes autocovariances lags 0,1,...,maxlag }
}
\details{
The algorithm given by McLeod (1975) is used.

The built-in R function ARMAacf could also be used but it is quite
complicated and apart from the source code, the precise algorithm
used is not described.  The only reference given for ARMAacf is the
Brockwell and Davis (1991) but this text does not give any detailed
exact algorithm for the general case.

Another advantage of TacvfAR over ARMAacf is that it will be easier
for to translate and implement this algorithm in other computing
environments such as MatLab etc. since the code is entirely written in R.
}
\value{
Vector of length = (lag.max+1) containing the autocovariances at lags 0,...,lag.max is returned.
}
\references{
McLeod, A.I. (1975), 
Derivation of the theoretical autocorrelation function of autoregressive 
moving-average time series. 
Applied Statistics, 24, 255-256.
} 
\author{ A.I. McLeod }


\seealso{ 
\code{\link{ARMAacf}}, 
\code{\link{InformationMatrixAR}}, 
\code{\link{GetARMeanMLE}}, 
\code{\link{SimulateGaussianAR}} 
 }
\examples{
#calculate and plot the autocorrelations from an AR(2) model
# with parameter vector c(1.8,-0.9).
 g<-TacvfAR(c(1.8,-0.9),20)
 AcfPlot(g/g[1], LagZeroQ=FALSE)

}
\keyword{ ts }
