\name{BoxCox.numeric}
\alias{BoxCox.numeric}
\title{ Box-Cox Analysis for a Time Series}
\description{
An AR(p) model is selected using AIC and then the best Box-Cox transformation
is determined.	Requires package FitAR.
}

\usage{
\method{BoxCox}{numeric}(object, interval = c(-1, 1), IIDQ = FALSE, ...)
}

\arguments{
  \item{object}{ a vector of time series values }
  \item{interval}{ interval to be searched }
  \item{IIDQ}{ If true, IID is assumed, ie. p=0. If FALSE, AR(p) is fit with p determined using AIC. }
  \item{\dots}{ optional arguments }
}

\details{
If the minimum data value is <= 0, a small positive constant, equal to
the negative of the minimum plus 0.25, is added to all the
data values.
If length(object) < 20, no AR model is used, that is, p=0.
}

\value{
No value returned. Graphical output produced as side-effect.
The plot shows relative likelihood funciton as well as the MLE and a confidence interval.
}

\references{ 
Box, G. E. P. and Cox, D. R. (1964) An analysis of transformations. 
Journal of Royal Statistical Society, Series B, vol. 26, pp. 211-246.
}

\author{ A.I. McLeod and Y. Zhang}

\note{ The MASS package has a similar function boxcox but this is implemented
only for regression and analysis of variance. }

\seealso{
 \code{\link{BoxCox.FitAR}},
 \code{\link{BoxCox.Arima}},
 \code{\link{BoxCox.ts}}
 }

\examples{
\dontrun{#takes a few seconds
#annual sunspot series
BoxCox(sunspot.year, IIDQ=FALSE)
#
#non-time series example, lengths of rivers
BoxCox(rivers)
}
}
\keyword{ ts }
