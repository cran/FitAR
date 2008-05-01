\name{GetFitARpLS}
\alias{GetFitARpLS}

\title{LS for AR(p) and Subset ARp -- Short Version}
\description{
For ARp subset models, the least squares estimates are
computed. The exact loglikelihood is then determined.
The estimated parameters are checked to see if they are in
the AR admissible region. 
}

\usage{
GetFitARpLS(z, pvec)
}

\arguments{
  \item{z}{ vector or ts object, the time series }
  \item{pvec}{ lags included in subset AR. If pvec = 0, white noise assumed. }
}
\details{
The R function \code{\link{lsfit}} is used.
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
  \item{yX}{ the y vector and X matrix used for the regression fitting}
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


\author{ A.I. McLeod }

\note{ 
This is a helper function for \code{\link{FitARp}}
which is invoked by the main package function \code{\link{FitAR}}.
Normally the user would \code{\link{FitAR}} since this
function provides generic print, summary, resid and plot methods but
GetFitARpLS is sometimes useful in iterative computations like
bootstrapping since it is faster.
 }

\seealso{ 
\code{\link{FitAR}}, 
\code{\link{FitARz}}, 
\code{\link{GetFitARz}}, 
\code{\link{FitARp}}, 
\code{\link{GetFitARpMLE}}, 
\code{\link{RacfPlot}} 
}

\examples{
#Fit subset AR using LS
data(SeriesA)
#normally use FitAR
ans<-FitAR(SeriesA, c(1,2,7), ARModel="ARp", MLEQ=FALSE)
#could also use FitARp
ans<-FitARp(SeriesA, c(1,2,7))
#for some applications GetFitARpLS is simpler and faster
ansLS<-GetFitARpLS(SeriesA, c(1,2,7))
ansLS
}

\keyword{ ts }


