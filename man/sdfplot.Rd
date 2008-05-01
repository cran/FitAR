\name{sdfplot}
\alias{sdfplot}
\title{ Autoregressive Spectral Density Estimation }
\description{
  Generic function. Methods are available for "FitAR", "ar", "Arima", "ts" and "numeric".
}
\usage{
sdfplot(obj, ...)
}
\arguments{
  \item{obj}{ input object }
  \item{\dots}{ optional arguments }
}
\value{
Plot produced as side-effect. No output.
}
\author{A.I. McLeod }

\seealso{ 
\code{\link{sdfplot}},
\code{\link{FitAR}}
}

\examples{
#Use AIC to select best subset model to fit to lynx data and
#plot spectral density function
data(SeriesA)
pvec<-SelectModel(SeriesA, ARModel="ARp", lag.max=10, Best=1)
ans<-FitAR(SeriesA, pvec)
sdfplot(ans)
#
# fit ARMA and plot sdf
ans<-arima(SeriesA, c(1,0,1))
sdfplot(ans)
}
\keyword{ ts }


