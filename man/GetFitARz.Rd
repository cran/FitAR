\name{GetFitARz}
\alias{GetFitARz}

\title{ Exact MLE for AR(p) and Subset ARz -- Short Version }
\description{
Obtain the exact MLE for AR(p) or subset ARz model.
This function is used by \code{FitAR} and \code{FitARz}.
One might prefer to use GetFitARz for applications such
as bootstrapping since it is faster.
}
\usage{
GetFitARz(z, pvec, MeanValue=0, ...)
}
\arguments{
  \item{z}{ time series }
  \item{pvec}{lags included in AR model. If pvec = 0, white noise model assumed. }
  \item{MeanValue}{by default it is assumed the mean of z is 0 }
  \item{\dots}{ optional arguments passed through to optim }
}
\details{
The built-in function optim is used to obtain the MLE estimates for
an AR or subset AR. First "BFGS" is tried.  This usually works fine.
In the rare cases where convergence is not obtained, "Nelder-Mead" is
used. A warning message is given if this happens.
}
\value{
  \item{loglikelihood}{value of maximized loglikelihood}
  \item{zetaHat}{estimated zeta parameters}
  \item{phiHat}{estimated phi parameters}
  \item{convergence}{result from optim}
  \item{pvec }{lags of estimated AR coefficient}
  \item{algorithm}{"BFGS" or "Nelder-Mead" }
}


\author{ A.I. McLeod and Y. Zhang}

\references{ 
McLeod, A.I. and Zhang, Y. (2006).  
Partial autocorrelation parameterization for subset autoregression.  
Journal of Time Series Analysis, 27, 599-612.

McLeod, A.I. and Zhang, Y. (2008a). Faster ARMA Maximum Likelihood Estimation, 
Computational Statistics and Data Analysis
52-4, 2166-2176.  
DOI link: \url{http://dx.doi.org/10.1016/j.csda.2007.07.020}. 

McLeod, A.I. and Zhang, Y. (2008b, Submitted).
Improved Subset Autoregression: With R Package.
Journal of Statistical Software.
}

\seealso{ 
\code{\link{FitAR}}, 
\code{\link{FitARz}}, 
\code{\link{FitARp}}, 
\code{\link{GetFitARpMLE}}, 
\code{\link{RacfPlot}} 
}

\examples{
#compare results from GetFitARz and FitAR
z<-log(lynx)
z<-z - mean(z)
GetFitARz(z, c(1,2,8))
out<-FitAR(log(lynx), c(1,2,8), ARModel="ARz")
out
coef(out)
}
\keyword{ ts }

