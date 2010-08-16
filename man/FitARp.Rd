\name{FitARp}
\alias{FitARp}
\title{ Fit subset ARp Models}
\description{
The subset ARp is defined as an AR(p) in which some
of the ar-coefficients are constrained to zero.
This is the usual type of subset AR.
In contrast the ARz model constrains some of the partial
autocorrelation coefficients to zero.
}
\usage{
FitARp(z, p, lag.max = "default", MLEQ = FALSE)
}
\arguments{
  \item{z}{ time series, vector or ts object }
  \item{p}{ p specifies the model. If length(p) is 1, an AR(p) is assumed and if p
  has length greater than 1, a subset ARp is assumed.
  For example, to fit a subset model with lags 1 and 4 present set p to c(1,4) or
  equivalently c(1,0,0,4).  To fit a subset model with just lag 4, you must use
  p=c(0,0,0,4) since p=4 will fit a full AR(4). }
  \item{lag.max}{ the residual autocorrelations are tabulated for lags 1, \dots, lag.max. Also
  lag.max is used for the Ljung-Box portmanteau test. }
  \item{MLEQ}{TRUE, use MLE. FALSE, use LS}
}

\details{
Subset ARp model is fit using exact MLE. The built-in \code{arima} function is 
used for MLE.  When MLEQ=FALSE, LS is used.  LS is has been widely used
in past for subset ARp fiting.
}


\value{
  A list with class name "FitAR" and components: 
  \item{loglikelihood }{value of the loglikelihood}
  \item{phiHat }{coefficients in AR(p) -- including 0's}
  \item{sigsqHat }{innovation variance estimate}
  \item{muHat }{estimate of the mean}
  \item{covHat }{covariance matrix of the coefficient estimates}
  \item{zetaHat }{transformed parameters, length(zetaHat) = \# coefficients estimated}
  \item{RacfMatrix }{residual autocorrelations and sd for lags 1, \dots, lag.max}
  \item{LjungBox}{table of Ljung-Box portmanteau test statistics}
  \item{SubsetQ }{parameters in AR(p) -- including 0's}
  \item{res}{innovation residuals, same length as z}
  \item{fits}{fitted values, same length as z}
  \item{pvec }{lags used in AR model}
  \item{demean }{TRUE if mean estimated otherwise assumed zero}
  \item{FitMethod }{"MLE" or "LS"}
  \item{IterationCount }{number of iterations in mean mle estimation}
  \item{convergence }{value returned by optim -- should be 0}
  \item{MLEMeanQ }{TRUE if mle for mean algorithm used}
  \item{ARModel}{"ARp" if FitARp used, otherwise "ARz"}
  \item{tsp}{tsp(z)}
  \item{call}{result from match.call() showing how the function was called}
  \item{ModelTitle}{description of model}
  \item{DataTitle}{returns attr(z,"title") }
  \item{z}{time series data input}
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
\code{\link{GetFitARpMLE}}, 
\code{\link{RacfPlot}} 
}

\examples{
#First Example: Fit to AR(4) 
set.seed(3323)
phi<-c(2.7607,-3.8106,2.6535,-0.9238)
z<-SimulateGaussianAR(phi,1000)
#MLE using arima
ans1<-FitARp(z,4,MLEQ=TRUE)
ans1
coef(ans1)
#OLS
ans2<-FitARp(z,4,MLEQ=FALSE)
ans2
coef(ans2)

\dontrun{#save time building package
#Second Example: Fit subset ARp model
z<-log(lynx)
#MLE
FitARp(z, c(1,2,4,7,10,11),MLEQ=TRUE)
#LS
FitARp(z, c(1,2,4,7,10,11),MLEQ=FALSE)

#Third Example: Use UBIC model selection to fit subset models
z<-log(lynx)
p<-SelectModel(z,ARModel="ARp")[[1]]$p
#MLE #error returned by arima
#ans1<-FitARp(z, p, MLEQ=TRUE)
#ans1
#LS
ans2<-FitARp(z, p, MLEQ=FALSE)
ans2
}
}

\keyword{ ts }
