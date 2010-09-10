\name{Caffeine}
\alias{Caffeine}
\docType{data}
\title{Caffeine industrial time series}

\description{
Hamilton and Watts (1978) state this series is produced from
a cyclic industrial process with a period of 5.
}

\usage{data(Caffeine)}

\format{
  The format is:
 num [1:178] 0.429 0.443 0.451 0.455 0.44 0.433 0.423 0.412 0.411 0.426 ...
}

\details{
The dataset are from the paper by Hamilton and Watts (1978, Table 1).
The series is used to illustrate how a multiplicative seasonal
ARMA model may be identified using the partial autocorrelations.
Chatfield (1979) argues that the inverse autocorrelations are
more effective for model identification with this example. 
}

\source{
Hamilton, David C. and Watts, Donald G. (1978).
Interpreting Partial Autocorrelation Functions of Seasonal Time Series Models.
Biometrika 65/1, 135-140.
}

\references{
Hamilton, David C. and Watts, Donald G. (1978).
Interpreting Partial Autocorrelation Functions of Seasonal Time Series Models.
Biometrika 65/1, 135-140.

Chatfield, C. (1979). Inverse Autocorrelations.
Journal of the Royal Statistical Society. Series A (General) 142/3, 363--377.
}

\examples{
#Example 1
sdfplot(Caffeine)
TimeSeriesPlot(Caffeine)
#
#Example 2
a<-numeric(3)
names(a)<-c("AIC", "BIC", paste(sep="","BIC(q=", paste(sep="",c(0.85),")")))
z<-Caffeine
lag.max <- ceiling(length(z)/4)
a[1]<-SelectModel(z, lag.max=lag.max, ARModel="AR", Best=1, Criterion="AIC")
a[2]<-SelectModel(z, lag.max=lag.max, ARModel="AR", Best=1, Criterion="BIC")
a[3]<-SelectModel(z, lag.max=lag.max, ARModel="AR", Best=1, Criterion="BICq", Q=0.85)
a
}
\keyword{datasets}
