\name{ARToMA}
\alias{ARToMA}
\title{Coefficients In Infinite Moving Average Expansion}
\description{
 A stationary-causal AR(p) can be written as a general linear process (GLP).
 This function obtains the moving-average expansion out to the L-th lag, 
 z[t] = a[t]+psi[1]*a[t-1]+...+psi[L]*a[t-L].
}
\usage{
ARToMA(phi, lag.max)
}
\arguments{
  \item{phi}{AR Coefficient vector}
  \item{lag.max}{maximum lag }
}
\details{
The coefficients are computed recursively as indicated in Box and Jenkins (1970).
}
\value{
Vector of length L+1 containing, (1,psi[1],...,psi[L])
}
\references{ Box and Jenkins (1970), Time Series Analysis, Forecasting & Control }
\author{ A.I. McLeod and Y. Zhang }


\seealso{ \code{\link{InvertibleQ}} }
\examples{
ARToMA(0.5,20)
ARToMA(c(0.2,0.5), 15)
}
\keyword{ ts }
