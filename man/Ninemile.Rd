\name{Ninemile}
\alias{Ninemile}
\docType{data}
\title{ Douglas Fir Treerings, Nine Mile Canyon, Utah, 1194-1964}
\description{
A treering time series comprises of 771 values showing a periodicity
of around 10 years.
}
\usage{data(Ninemile)}
\format{
ts object with title attribute
}

\source{
Hipel, K.W. and McLeod, A.I. (2006).
Time Series Modelling of Water Resources and Environmental Systems.
}
\references{ 
McLeod, A.I. and Zhang, Y. (2006).  
Partial autocorrelation parameterization for subset autoregression.  
Journal of Time Series Analysis, 27, 599-612.
}


\examples{
ans<-FitAR(Ninemile, c(1,2,9))
summary(ans)

}
\keyword{datasets}




