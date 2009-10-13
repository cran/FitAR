\name{BoxCox.Arima}
\alias{BoxCox.Arima}
\title{ Box-Cox Analysis for "Arima" Objects }
\description{
Implements Box-Cox analysis for "Arima" class objects, the output from \code{arima}, a \R built-in function.
Variance change in time series is an important topic. 
In some cases using a Box-Cox transformation will provide a much simpler
analysis than the much more complex ARMA-GARCH approach.
See US Tobacco series example given below for an example. 
}
\usage{
\method{BoxCox}{Arima}(object, interval = c(-1, 1), type = "BoxCox", InitLambda = "none", ...)
}
\arguments{
  \item{object}{ output from arima, a R built-in function }
  \item{interval}{ interval to be searched for the optimal transformation }
  \item{type}{ ignored unless, InitLambda!="none". Type of transformation, default is "BoxCox". Otherwise a simple power transformation. }
  \item{InitLambda}{ default "none". Otherwise a numerical value giving the transformation parameter. }
  \item{\dots}{ optional arguments passed to optimize }
}
\details{
If no transformation is used on the data, then the original data is used.
But if a transformation has already been used, we need to inverse transform
the data to recover the untransformed data.

For \eqn{\lambda \ne 0}{lambda!=0}, the Box-Cox transformation of x is
\eqn{(x-1)^\lambda/\lambda}{(x-1)^lambda/lambda} whereas the regular power transformation 
is simply \eqn{x^\lambda}{x^lambda}.  
When \eqn{\lambda=0}{lambda=0}, it is log in both cases.

The log of the Jacobian is 
\eqn{(\lambda-1)\sum_{t=D+1}^n \log(z_t)}{%
(lambda-1)*sum(log(z[(D+1):n])},
where \eqn{\lambda}{lambda} is the transformation,
n=length(z), 
z is the vector of data and D = d + ds*s,
where d is the degree of regular differencing, ds is the degree of seasonal
differencing and s is the seasonal period.
The correct expression for the loglikelihood function was first given
in Hipel and McLeod (1977, eqn. 10).
Using the wrong expression for the Jacobian has a disasterous effect in
many situations.  For example with the international airline passenger time
series, the MLE for lambda would be about 1.958 instead of close to zero.

If the minimum data value is <= 0, a small positive constant, equal to
the negative of the minimum plus 0.25, is added to all the
data values.

}

\value{
No value returned. Graphical output is produced as side-effect.
The plot shows relative likelihood function as well as the MLE and a confidence interval.
}

\references{ 
Hipel, K.W. and McLeod, A.I. (1977). Advances in Box-Jenkins Modelling. 
Part 1, Model Construction. Water Resources Research 13, 567-575.
}

\author{ A.I. McLeod and Y. Zhang}

\note{ The MASS package has a similar function \code{boxcox} but this is implemented
only for regression and analysis of variance. }

\seealso{
 \code{\link{arima}},
 \code{\link{BoxCox}}, 
 \code{\link{BoxCox.FitAR}} 
}

\examples{
#Tobacco Production
 data(USTobacco)
 plot(USTobacco)
 USTobacco.arima<-arima(USTobacco,order=c(0,1,1))
 BoxCox(USTobacco.arima)
#
 air.arima<-arima(AirPassengers, c(0,1,1), seasonal=list(order=c(0,1,1), period=12))
 BoxCox(air.arima)
#
#In this example, we fit a model to the square-root of the sunspots and
#back transform in BoxCox.
sqrtsun.arima<-arima(sqrt(sunspot.year),c(2,0,0))
BoxCox(sqrtsun.arima, InitLambda=0.5, type="power")
#
#Back transform with AirPassengers
Garima<-arima(log(AirPassengers), c(0,1,1), seasonal=list(order=c(0,1,1),period=12))
BoxCox(Garima, InitLambda=0)
}

\keyword{ ts }




