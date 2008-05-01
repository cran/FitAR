\name{GetFitARpMLE}
\alias{GetFitARpMLE}

\title{ Exact MLE for subset ARp Models}
\description{
Uses built-in function \code{arima} to fit subset ARp model, that is,
the subset model is formed by constraining some coefficients to zero.
}
\usage{
GetFitARpMLE(z, pvec)
}
\arguments{
  \item{z}{ time series }
  \item{pvec}{lags included in AR model. If pvec = 0, white noise model assumed. }
}
\details{
Due to the optimization algorithms used by \code{arima}, this method
is not very reliable.  The optimization may simply fail. 
Example 1 shows it working but in Example 2 below it fails.
}

\value{
a list with components:
  \item{loglikeliihood }{the exact loglikelihood} 
  \item{phiHat }{estimated AR parameters}
  \item{constantTerm }{constant term in the linear regression}
  \item{pvec }{lags of estimated AR coefficient}
  \item{res }{the least squares regression residuals}
  \item{InvertibleQ}{True, if the estimated parameters are in the
     AR admissible region.}
}

\references{ 
McLeod, A.I. and Zhang, Y. (2006).  
Partial Autocorrelation Parameterization for Subset Autoregression.  
Journal of Time Series Analysis, 27, 599-612.

McLeod, A.I. and Zhang, Y. (2008a). Faster ARMA Maximum Likelihood Estimation, 
Computational Statistics and Data Analysis
52-4, 2166-2176.  
DOI link: \url{http://dx.doi.org/10.1016/j.csda.2007.07.020}. 

McLeod, A.I. and Zhang, Y. (2008b, Submitted).
Improved Subset Autoregression: With R Package.
Journal of Statistical Software.
}


\author{ A.I. McLeod }

\seealso{ 
\code{\link{FitAR}}, 
\code{\link{FitARz}}, 
\code{\link{GetFitARz}}, 
\code{\link{FitARp}}, 
\code{\link{RacfPlot}} 
}
\examples{
#Example 1. MLE works
z<-log(lynx)
p<-c(1,2,4,7,10,11)
GetFitARpMLE(z, p)
#
#Example 2. MLE fails with error.
p<-c(1,2,9,12)
\dontrun{GetFitARpMLE(z, p)}
}

\keyword{ ts }
