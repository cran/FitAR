\name{FitAR}
\alias{FitAR}
\title{ Fit AR, ARp and ARz }
\description{
Exact MLE for full AR as well as subset AR. Both subset ARp and subset ARz models
are implemented.  For subset ARp models the R function arima is used.
For full AR and subset ARz models, algorithm of McLeod & Zhang (2006) is implemented.
The LS algorithm for subset ARp is also available as an option.

}
\usage{
FitAR(z, p, lag.max = "default", ARModel = "ARz", ...)
}

\arguments{
  \item{z}{ time series, vector or ts object. }
  \item{p}{ p specifies the model. If length(p) is 1, an AR(p) is assumed and if p
  has length greater than 1, a subset ARp or ARz is assumed - the default is ARz.
  For example, to fit a subset model with lags 1 and 4 present, set p to c(1,4) or
  equivalently c(1,0,0,4).  To fit a subset model with just lag 4, you must use
  p=c(0,0,0,4) since p=4 will fit a full AR(4). }
  \item{lag.max}{ the residual autocorrelations are tabulated for lags 1, \dots, lag.max. Also
  lag.max is used for the Ljung-Box portmanteau test. }
  \item{ARModel}{ which subset model, ARz or ARp }
  \item{\dots}{ optional arguments which are passed to
                \code{FitARz} or \code{FitARp} }
}


\details{
The exact MLE for AR(p) and subset ARz use methods described in
McLeod and Zhang (2006).
In addition the exact MLE for the mean can be computed using an
iterative backfitting approach described in McLeod and Zhang (2008).

The subset ARp model can be fit by exact MLE using the R function \code{arima}
or by least-squares.

The default for lag.max is
min(300, ceiling(length(z)/5))
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
  \item{pvec}{lags used in AR model}
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
Computational Statistics and Data Analysis,
52-4, 2166-2176.  
DOI link: \url{http://dx.doi.org/10.1016/j.csda.2007.07.020}. 

McLeod, A.I. and Zhang, Y. (2008b, Submitted).
Improved Subset Autoregression: With R Package.
Journal of Statistical Software.
}


\author{ A.I. McLeod }

\note{ 
There are generic print, summary, coef and resid functions for class "FitAR".

It is somewhat surprising that in the 'ARp' subset autoregression
quite different subsets may be chosen depending on the choice of 'lag.max'.
For example, with the 'lynx' taking lag.max = 15, 20
produces subsets {1, 2, 4, 10, 11} and {1, 2, 10, 11} using the BIC.
This also occurs even with the AIC.
See sixth example below.
}



\seealso{ 
\code{\link{FitARp}}, 
\code{\link{FitARz}}, 
\code{\link{GetFitARz}}, 
\code{\link{FitARp}}, 
\code{\link{GetFitARpMLE}}, 
\code{\link{RacfPlot}} 
}

\examples{
#First example: fit exact MLE to AR(4) 
set.seed(3323)
phi<-c(2.7607,-3.8106,2.6535,-0.9238)
z<-SimulateGaussianAR(phi,1000)
ans<-FitAR(z,4,MeanMLEQ=TRUE)
ans
coef(ans)

\dontrun{ #save time building package!
#Second example: compare with sample mean result
ans<-FitAR(z,4)
coef(ans)

#Third example: fit subset ARz and ARp models
z<-log(lynx)
FitAR(z, c(1,2,4,7,10,11))
#now obtain exact MLE for Mean as well
FitAR(z, c(1,2,4,7,10,11), MeanMLE=TRUE)
#subset ARp using exact MLE
FitAR(z, c(1,2,4,7,10,11), ARModel="ARp", MLEQ=TRUE)
#subset ARp using LS
FitAR(z, c(1,2,4,7,10,11), ARModel="ARp", MLEQ=FALSE)
#or
FitAR(z, c(1,2,4,7,10,11), ARModel="ARp")


#Fourth example: use UBIC model selection to fit subset models
z<-log(lynx)
#ARz case
p<-SelectModel(z,ARModel="ARz")[[1]]$p
ans1<-FitAR(z, p)
ans1
ans1$ARModel

#ARp case
p<-SelectModel(z,ARModel="ARp")[[1]]$p
ans2<-FitAR(z, p, ARModel="ARp")
ans2
ans2$ARModel

#Fifth example: fit a full AR(p) using AIC/BIC methods
z<-log(lynx)
#BIC
p<-SelectModel(z,ARModel="AR")[1,1]
ans1<-FitAR(z, p)
ans1
#AIC
p<-SelectModel(z, ARModel="AR", Criterion="AIC")[1,1]
ans2<-FitAR(z, p)
ans2
}

#Sixth Example: Subset autoregression depends on lag.max!
#Because least-squares is used, P=lag.max observations are
#  are deleted. This causes different results depending on lag.max.
#This phenomenon does not happen with "ARz" subset models
#ARp models depend on lag.max
SelectModel(z,lag.max=15,ARModel="ARp", Criterion="BIC", Best=1)
SelectModel(z,lag.max=20,ARModel="ARp", Criterion="BIC", Best=1)
#ARz models do NOT depend in this way on lag.max.
#Obviously if some lags beyond the initial value of lag.max are
# found to be important, then there is a dependence but this
# is not a problem!
SelectModel(z,lag.max=15,ARModel="ARz", Criterion="BIC", Best=1)
SelectModel(z,lag.max=20,ARModel="ARz", Criterion="BIC", Best=1)
}

\keyword{ ts }

