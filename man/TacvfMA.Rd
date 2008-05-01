\name{TacvfMA}
\alias{TacvfMA}
\title{ Theoretical Autocovariances for Moving Average Process }
\description{
The theoretical autocovariance function of a MA(q) with unit
variance is computed.

}
\usage{
TacvfMA(theta, lag.max = 20)
}
\arguments{
  \item{theta}{ q parameters in MA(q)}
  \item{lag.max}{ number of lags required. }
}
\details{
The first q+1 values are determined using a matrix multiplication -
avoiding a loop. The remaining values set to zero.
}
\value{
Vector of length q+1 containing the autocovariances at lags 0,1,...,lag.max
}
\references{
McLeod, A.I. and Zhang, Y. (2006),  
Partial autocorrelation parameterization for subset autoregression.  
\emph{Journal of Time Series Analysis}, \bold{27}, 599-612.
}

\author{ A.I.McLeod and Y. Zhang}
\note{See Details in \code{\link{TacvfAR}} for why we prefer
to use this algorithm instead of \code{\link{ARMAacf}}  }

\seealso{ 
\code{\link{ARMAacf}}, 
\code{\link{TacvfAR}} 
}
\examples{
TacvfMA(c(1.8,-0.9), 10)
}
\keyword{ ts }
