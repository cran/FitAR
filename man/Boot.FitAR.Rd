\name{Boot.FitAR}
\alias{Boot.FitAR}
\title{ Simulate a Fitted AR }
\description{
Simulate a realization from a fitted AR model. This is useful in the
parametric bootstrap.  Generic function for "Boot" method.
}
\usage{
\method{Boot}{FitAR}(obj, R=1, ...)
}
\arguments{
  \item{obj}{ the output from FitAR }
  \item{R}{ number of bootstrap replications}
  \item{...}{ optional arguments }
}
\value{
A simulated time series with the same length as the original fitted time series is produced when R=1.
When R>1, a matrix with R columns is produced with each column a separate bootstrap realization.
}
\author{ A.I. McLeod and Y. Zhang}

\seealso{  
\code{\link{Boot}} 
\code{\link{SimulateGaussianAR}} 
}
\examples{
#Plot log(lynx) time series and simulation
#
data(lynx)
ans <- FitAR(log(lynx), 8)
z<-Boot.FitAR(ans)
par(mfrow=c(2,1))
TimeSeriesPlot(log(lynx))
title(main="log(lynx) time series")
TimeSeriesPlot(z)
title(main="Simulated AR(8), fitted to log lynx")
#par(mfrow=c(1,1)
#
#Use bootstrap to compute standard errors of parameters
#takes about 18 seconds on a 3.6 GHz PC
ptm <- proc.time() #user time
R<-100  #number of bootstrap iterations
p<-c(1,2,4,7,10,11)
ans<-FitAR(log(lynx),p)
out<-Boot(ans, R)
fn<-function(z) GetFitARz(z,p)$zetaHat
sdBoot<-sqrt(diag(var(t(apply(out,fn,MARGIN=2)))))
sdLargeSample<-coef(ans)[,2][1:6]
sd<-matrix(c(sdBoot,sdLargeSample),ncol=2)
dimnames(sd)<-list(names(sdLargeSample),c("Bootstrap","LargeSample"))
ptm<-(proc.time()-ptm)[1]
sd
}
\keyword{ts}
