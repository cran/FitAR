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
Plot is produced using plot.
Matrix with 2 columns containing the frequencies and spectral density is returned invisibly.
}
\author{A.I. McLeod }

\seealso{ 
\code{\link{sdfplot}},
\code{\link{FitAR}}
}

\examples{
#Example 1
#Use AIC to select best subset model to fit to lynx data and
#plot spectral density function
pvec<-SelectModel(SeriesA, ARModel="ARp", lag.max=10, Best=1)
ans<-FitAR(SeriesA, pvec)
sdfplot(ans)
#
#Example 2
#Fit ARMA and plot sdf
ans<-arima(SeriesA, c(1,0,1))
sdfplot(ans)
#
#Example 3
#Fit ARz model using AIC to monthly sunspots and plot spectral density
#Warning: this may take 10 minutes or so.
\dontrun{
pvec<-SelectModel(sunspots, lag.max=200, ARModel="ARz", Criterion="AIC", Best=1)
ans<-FitAR(sunspots, pvec)
sdfplot(ans)
}
}
\keyword{ ts }


