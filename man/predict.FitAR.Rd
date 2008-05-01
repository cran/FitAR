\name{predict.FitAR}
\alias{predict.FitAR}
\title{ Predict Subset AR Model }
\description{
After fitting we predict at origin times n, n+1, ..., n+m,
where m is the length of the vector newdata and for
lead time series as specified by n.ahead.
}
\usage{
\method{predict}{FitAR}(object, n.ahead = 1, newdata = numeric(0), ...)
}
\arguments{
  \item{object}{ \sQuote{FitAR} object }
  \item{n.ahead}{ lead time }
  \item{newdata}{ new time series values }
  \item{\dots}{ optional arguments }
}
\details{
The prediction algorithm described in McLeod, Yu and Zinovi (2008) is used.  
}

\value{
A list with components
  \item{Forecasts }{matrix with m+1 rows and maxLead columns with the forecasts}
  \item{SDForecasts }{matrix with m+1 rows and maxLead columns with the sd of the forecasts}
}

\references{ 
McLeod AI, Yu H, Zinovi K (2008). 
Linear Time Series Modeling with R Package. 
Journal of Statistical Software, 23/5, 1-26.
}

\author{ A.I. McLeod }

\seealso{ \code{\link[ltsa]{TrenchForecast}}
}

\examples{
#Example 1.  
#Compare the predictions for the monthly sunspots using the ARz
#  fitted using the UBIC and BIC.
#  This computation takes about 3-4 minutes.

`getRMSE` <-
function(obj, zTOT, n.ahead = 12, newdata=numeric(0)){
    ans<-predict(obj, n.ahead=n.ahead, newdata=newdata)
    ansf<-ans$Forecasts
    nL<-as.numeric(colnames(ansf))
    nO<-as.numeric(rownames(ansf))
    err<-ansf-zTOT[-1+outer(nO,nL,FUN="+")]
    s<-apply(err, MARGIN=2, FUN=rmse)
    s
}

`rmse` <-
function(x){
y<-x[!is.na(x)]
sqrt(sum(y^2)/length(y))
}

zTOT <- sqrt(sunspots)
nTOT <- length(zTOT)
nOUT <- 12*3 #using last 3 years for out-of-sample forecasts
ind<- (1:nTOT)<(nTOT-nOUT+1)
newdata<-zTOT[!ind]
z<-zTOT[ind]
lag.max<-12*11 #using lags up to last 11 years in subset model
nahead<-4 #forecasts for 1 to 4 months ahead
pUBIC <- SelectModel(z, ARModel="ARz", lag.max=lag.max, Best=1)
zUBIC <- FitAR(z, pUBIC, ARModel="ARz")
pBIC <- SelectModel(z, ARModel="ARz",lag.max=lag.max,Best=1,Criterion="BIC")
zBIC <- FitAR(z, pBIC, ARModel="ARz")
fubic<-getRMSE(zUBIC, zTOT, n.ahead=nahead, newdata=newdata)
fbic<-getRMSE(zBIC, zTOT, n.ahead=nahead, newdata=newdata)
m<-matrix(c(fubic,fbic), ncol=2)
dimnames(m)<-list(1:nahead, c("fubic","fbic")) 
m
#
#Example 2.
#Compute predictions and plot observed - predicted
z <- sqrt(sunspots)
pUBIC <- SelectModel(z, ARModel="ARz", lag.max=240, Best=1)
zUBIC <- FitAR(z, pUBIC, ARModel="ARz")
out<-predict(zUBIC, n.ahead=24)
zf<-out$Forecasts
zsd<-out$SDForecasts
y<-cts(z, zf)
plot(window(y,start=1980), type="n", ylab="sqrt sunspot number")
y1<-window(y, start=1980, end=1983)
lines.ts(y1,col="blue",type="o", lwd=2, pch=16)
y2<-window(y, start=c(1983,1))
lines.ts(y2,col="red",type="o",lwd=2, pch=16)
legend(1984,12, legend=c("observed", "forecast"),col=c("red","blue"),lwd=c(2,2),pch=c(16,16))
}

\keyword{ ts }
