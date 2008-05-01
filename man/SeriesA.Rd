\name{SeriesA}
\alias{SeriesA}
\docType{data}
\title{ Series A, Chemical Process Concentration Readings}
\description{
Chemical process concentration readings for every 2 hours.
}
\usage{data(SeriesA)}
\format{
ts object with attribute "title"
}
\details{
Box and Jenkins (1970) fit an ARMA(1,1) and ARIMA(0,1,1) to this series.
Cleveland(1971) suggested a subset AR(1,2,7).
McLeod and Zhang (2006) fit a subset ARz(1,2,6,7) parameterized
with the partial autocorrelations.
}
\source{
Box and Jenkins (1970). 
Time Series Analysis: Forecasting and Control.
}
\references{
Cleveland, W.S. (1971)
The inverse autocorrelations of a time series and their applications.
Technometrics 14, 277-298.

McLeod, A.I. and Zhang, Y. (2006).  
Partial autocorrelation parameterization for subset autoregression.  
Journal of Time Series Analysis, 27, 599-612.

} 
\examples{
data(SeriesA)
#fit subset models
FitAR(SeriesA, c(1,2,7), ARModel="ARp")
FitAR(SeriesA, c(1,2,6,7), ARModel="ARz")
}
\keyword{datasets}
