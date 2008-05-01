\name{USTobacco}
\alias{USTobacco}
\docType{data}
\title{ U.S. Tobacco Production, 1871-1984}
\description{
 Annual U.S. tobacco production, 1871-1984, in millions of pounds.
}
\usage{data(USTobacco)}
\format{
  The format is:
 Time-Series [1:114] from 1871 to 1984: 327 385 382 217 609 466 621 455 472 469 ...
 - attr(*, "title")= chr "Tobacco production,US, 1871-1984"
}
\details{
Wei (2006, p.120, Example 6.6) fits an ARIMA(0,1,1) to the logarithms.
But a more accurate Box-Cox analysis indicates a square-root transformation
should be used.  A more complex ARIMA-GARCH model is also suggested by Wei (2006).
}
\source{
Wei, W.W.S. (2006, Series W6, p.570), \emph{Time Series Analysis: Univariate and Multivariate Methods}.
2nd Ed., New York: Addison-Wesley.

}

\examples{
#From a plot of the series, we see that the variance is increasing with level.
#From the acf of the first differences an ARIMA(0,1,1) is suggested.
 data(USTobacco)
# layout(matrix(c(1,2,1,2),ncol=2))
 plot(USTobacco)
 lines(lowess(time(USTobacco), USTobacco), lwd=2, col="blue")
 acf(diff(USTobacco, differences=1))

}
\keyword{datasets}
