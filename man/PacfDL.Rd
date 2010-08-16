\name{PacfDL}
\alias{PacfDL}
\title{ Partial Autocorrelations via Durbin-Levinson }
\description{
Given autocovariances, the partial autocorrelations 
and/or autoregressive coefficients in an AR
may be determined
using the Durbin-Levinson algorithm. If the autocovariances are sample
autocovariances, this is equivalent to using the Yule-Walker equations.
But as noted below our function is more general than the built-in R
functions.  
}
\usage{
PacfDL(c, LinearPredictor = FALSE)
}
\arguments{
  \item{c}{ autocovariances at lags 0,1,\dots,p = length(c)-1 }
  \item{LinearPredictor}{ if TRUE, AR coefficients are also determined
  using the Yule-Walker method }
}
\details{
The Durbin-Levinson algorithm is described in many books on time
series and numerical methods, for example Percival and Walden (1993,
eqn 403).    
}
\value{
If LinearPredictor = FALSE, vector of length p = length(c)-1 containing the partial autocorrelations
at lags 1,\dots,p. Otherwise a list with components:
  \item{Pacf }{vector of partial autocorrelations}
  \item{ARCoefficients }{vector of AR coefficients}
  \item{ResidualVariance }{residual variance for AR(p)}

}
\references{
Percival, D.B. and Walden, A.T. (1993).  
Spectral Analysis For Physical Applications, 
Cambridge University Press.  
}
\author{A.I. McLeod and Y. Zhang }

\section{Warning }{
Stationarity is not tested.
} 

\note{Sample partial autocorrelations can also be computed with the
\code{\link{acf}} function and Yule-Walker estimates can
be computed with the \code{\link{ar}} function.
Our function \code{\link{PacfDL}} provides more flexibility since
then input c may be any valid autocovariances not just the usual
sample autocovariances.  For example, we can determine the minimum
mean square error one-step ahead linear predictor of order p for
theoretical autocovariances from a fractional arma or other linear
process. 
}

\seealso{ 
\code{\link{acf}}, 
\code{\link{ar}}

}
\examples{
#first define a function to compute the Sample Autocovariances
 sacvf<-function(z, lag.max){
 c(acf(z, plot=FALSE, lag.max=lag.max)$acf)*(length(z)-1)/length(z)
 }
#now compute PACF and also fit AR(7) to SeriesA
 ck<-sacvf(SeriesA, 7)
 PacfDL(ck)
 PacfDL(ck, LinearPredictor = TRUE)
#compare with built-in functions
 pacf(SeriesA, lag.max=7, plot=FALSE)
 ar(SeriesA, lag.max=7, method="yw")
#fit an optimal linear predictor of order 10 to MA(1)
 g<-TacvfMA(0.8,5)
 PacfDL(g, LinearPredictor=TRUE)
#
#Compute the theoretical pacf for MA(1) and plot it
 ck<-c(1,-0.4,rep(0,18))
 AcfPlot(PacfDL(ck)$Pacf)
 title(main="Pacf of MA(1), r(1)=-0.4")

}
\keyword{ ts }
