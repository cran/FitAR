\name{BoxCox.ts}
\alias{BoxCox.ts}
\title{ Box-Cox Analysis for a Time Series }
\description{
  The time series is converted to a vector and BoxCox.numeric is used.
}
\usage{
\method{BoxCox}{ts}(object, interval = c(-1, 1), ...)
}

\arguments{
  \item{object}{ a vector of time series values }
  \item{interval}{ interval to be searched }
  \item{\dots}{ optional arguments  }
}

\details{
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
}

\author{ A.I. McLeod }

\note{ The MASS package has a similar function boxcox but this is implemented
only for regression and analysis of variance. }

\section{Warning}{
It is important not to transform the data when fitting it with
AR since the optimal transformation would be found for the transformed
data -- not the original data. Normally this would not be a sensible thing
to do.
}

\seealso{
 \code{\link{BoxCox.FitAR}},
 \code{\link{BoxCox.Arima}},
 \code{\link{BoxCox.numeric}}
 }

\examples{
#
\dontrun{#takes a few seconds
BoxCox(sunspot.year)
}
}
\keyword{ ts }

