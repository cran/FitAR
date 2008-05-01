\name{BackcastResidualsAR}
\alias{BackcastResidualsAR}
\title{Innovation Residuals in AR }
\description{
Obtains the residuals (estimated innovations).  The residuals for t=1,...,p
are obtained using the backforecasting algorithm of Box and Jenkins (1970).
}
\usage{
BackcastResidualsAR(y, phi, Q = 100, demean=TRUE)
}
\arguments{
  \item{y}{a time series or vector}
  \item{phi}{AR coefficients, lags 1,...,p}
  \item{Q}{for backcasting, the AR is approximated by an MA(Q)}
  \item{demean}{ subtract sample mean }
}
\details{
The backforecasting algorithm is described in detail in the book of Box and Jenkins (1970).
The idea is to compute the expected value of the innovation assuming a high-order MA(q).
}
\value{
Vector of residuals
}
\references{Box and Jenkins (1970). Time Series Analysis: Forecasting and Control.} 
\author{A.I. McLeod and Y. Zhang}

\note{No check is done that the AR is causal-stationary. }
\seealso{ \code{\link{InvertibleQ}},
\code{\link{FitAR}} 
}
\examples{
#compare residuals obtained using backcasting with fitted parameters and
# the residuals extracted from output of FitAR.  They are identical.
p<-11
out<-FitAR(log(lynx), p)
phi<-out$phiHat #fitted parameters
resphi<-BackcastResidualsAR(log(lynx), phi)
sum(abs(resphi-resid(out)))
}
\keyword{ ts }
