\name{GetFitAR}
\alias{GetFitAR}
\title{MLE for AR, ARp and ARz}
\description{
Obtains the exact MLE for AR(p) or subset AR models ARp or ARz.
This function is used by FitAR.
One might prefer to use GetFitAR for applications such
as bootstrapping since it is faster than FitAR.
}

\usage{
GetFitAR(z, p, ARModel = "ARz", ...)
}

\arguments{
  \item{z}{ time series }
  \item{p}{ model order or subset lags }
  \item{ARModel}{ either "ARp" or "ARz"  corresponding to \code{GetFitARp} or \code{GetFitARz}}
  \item{\dots}{ optional arguments which are passed to \code{GetFitARp} or \code{GetFitARz} }
}

\details{
This is just a shell which simply invokes either
\code{GetFitARp} or \code{GetFitARz}
}

\value{
  \item{loglikelihood}{value of maximized loglikelihood}
  \item{zetaHat}{estimated zeta parameters}
  \item{phiHat}{estimated phi parameters}
  \item{convergence}{result from optim}
}

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

\author{ A.I. McLeod}

\seealso{ \code{\link{FitAR}} }

\examples{
#compare results from GetFitAR and FitAR
z<-log(lynx)
z<-z - mean(z)
GetFitAR(z, c(1,2,8))
out<-FitAR(log(lynx), c(1,2,8))
out
coef(out)
}
\keyword{ ts }
