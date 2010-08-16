\name{BoxCox.FitAR}
\alias{BoxCox.FitAR}
\title{ Box-Cox Analysis for "FitAR" Objects }
\description{
  This is a methods function to do a Box-Cox analysis for models
  fit using FitAR.  
}

\usage{
\method{BoxCox}{FitAR}(object, interval = c(-1, 1), type = "BoxCox", InitLambda = "none", ...)
}

\arguments{
  \item{object}{ output from FitAR }
  \item{interval}{ interval to be searched for the optimal transformation }
  \item{type}{ Ignored unless, InitLambda!="none".  Type of transformation, default is "BoxCox". Otherwise a simple power transformation. }
  \item{InitLambda}{ default "none". Otherwise a numerical value giving the transformation parameter. }
  \item{\dots}{ optional arguments passed to optimize }
}

\details{
If no transformation is used on the data, then the original data is used.
But if a transformation has already been used, we need to inverse transform
the data to recover the untransformed data.

For \eqn{\lambda \ne 0}{lambda!=0}, the Box-Cox transformation is of x is
\eqn{(x-1)^\lambda/\lambda}{(x-1)^lambda/lambda} whereas the regular power transformation 
is simply \eqn{x^\lambda}{x^lambda}.  
When \eqn{\lambda=0}{lambda=0}, it is log in both cases.

If the minimum data value is <= 0, a small positive constant, equal to
the negative of the minimum plus 0.25, is added to all the
data values.
}


\value{
No value returned. Graphical output produced as side-effect.
The plot shows relative likelihood function as well as the MLE and a confidence interval.
}

\references{  
Box, G. E. P. and Cox, D. R. (1964) An analysis of transformations. 
Journal of Royal Statistical Society, Series B, vol. 26, pp. 211-246.

McLeod, A.I. and Zhang, Y. (2006a).  
Partial autocorrelation parameterization for subset autoregression.  
Journal of Time Series Analysis, 27, 599-612.

McLeod, A.I. and Zhang, Y. (2006b, under review).  
Subset Autoregression Modelling.
Journal of Statistical Sofware.
}

\author{ A.I. McLeod }

\note{ The MASS package has a similar function \code{boxcox} but this is implemented
only for regression and analysis of variance. }


\seealso{
 \code{\link{BoxCox}}, 
 \code{\link{BoxCox.Arima}} 
 }

\examples{
\dontrun{ #takes a few seconds
#lynx time series. ARp subset model.
out<-FitAR(lynx, c(1,2,4,10,11), ARModel="ARp")
BoxCox(out)
#
#lynx time series. ARz subset model.
p<-SelectModel(lynx, ARModel="ARz", lag.max=25, Best=1)
out<-FitAR(lynx, p)
BoxCox(out)
}

}

\keyword{ ts }

