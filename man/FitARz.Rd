\name{FitARz}
\alias{FitARz}
\title{ Subset ARz Model Fitting }
\description{
The subset ARz model, defined by constraining partial autocorrelations to zero,
is fitted using exact MLE.
When length(p)=1, an AR(p) is fit by MLE.
}
\usage{
FitARz(z, p, demean = TRUE, MeanMLEQ = FALSE, lag.max = "default")
}
\arguments{
  \item{z}{ time series, vector or ts object }
  \item{p}{ p specifies the model. If length(p) is 1, an AR(p) is assumed and if p
  has length greater than 1, a subset ARz is assumed.
  For example, to fit a subset model with lags 1 and 4 present set p to c(1,4) or
  equivalently c(1,0,0,4).  To fit a subset model with just lag 4, you must use
  p=c(0,0,0,4) since p=4 will fit a full AR(4).}
  \item{demean}{TRUE, mean estimated. FALSE, mean is zero. }
  \item{MeanMLEQ}{use exact MLE for mean parameter }
  \item{lag.max}{ the residual autocorrelations are tabulated for lags 1, \dots, lag.max. Also
  lag.max is used for the Ljung-Box portmanteau test. }
}
\details{
The model and its properties are discussed in McLeod and Zhang (2006)
and McLeod and Zhang (2008).
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
  \item{z}{time series data input) }
}

\references{ 
McLeod, A.I. and Zhang, Y. (2006).  
Partial Autocorrelation Parameterization for Subset Autoregression.  
Journal of Time Series Analysis, 27, 599-612.

McLeod, A.I. and Zhang, Y. (2008a). Faster ARMA Maximum Likelihood Estimation, 
Computational Statistics and Data Analysis,
52-4, 2166-2176.  
DOI link: \url{http://dx.doi.org/10.1016/j.csda.2007.07.020}. 

McLeod, A.I. and Zhang, Y. (2008b, Submitted).
Improved Subset Autoregression: With R Package.
Journal of Statistical Software.
}


\author{ A.I. McLeod }

\note{
Normally one would use the \code{FitAR} function which
then calls this function for the ARz case.
}
\seealso{ 
\code{\link{FitAR}}, 
\code{\link{FitARp}}, 
\code{\link{GetFitARz}}, 
\code{\link{GetFitARpMLE}}, 
\code{\link{RacfPlot}} 
}

\examples{
#First Example: Fit exact MLE to AR(4) 
set.seed(3323)
phi<-c(2.7607,-3.8106,2.6535,-0.9238)
z<-SimulateGaussianAR(phi,1000)
ans<-FitARz(z,4,MeanMLEQ=TRUE)
ans
coef(ans)

#Second Example: compare with sample mean result
ans<-FitARz(z,4)
coef(ans)

#Third Example: Fit subset ARz 
z<-log(lynx)
FitARz(z, c(1,2,4,7,10,11))
#now obain exact MLE for Mean as well
FitARz(z, c(1,2,4,7,10,11), MeanMLE=TRUE)

#Fourth Example: Fit subset ARz
somePACF<-c(0.5,0,0,0,-0.9)
someAR<-PacfToAR(somePACF)
z<-SimulateGaussianAR(someAR,1000)
ans=FitARz(z, c(1,5),MeanMLEQ=TRUE)
coef(ans)
GetFitARz(z,c(1,5))#assuming a known zero mean
}

\keyword{ ts }
