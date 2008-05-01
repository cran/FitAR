\name{ChampernowneD}
\alias{ChampernowneD}
\title{ Champernowne Matrix }
\description{
Computes sufficient statistics for AR.
}
\usage{
ChampernowneD(z, p, MeanZero = FALSE)
}
\arguments{
  \item{z}{ time series data}
  \item{p}{ order of the AR }
  \item{MeanZero}{ Assume mean is zero. Default is FALSE so the sample mean
  is subtracted from the data first. Otherwise no sample mean correction is made. }
}
\details{
This matrix is defined in McLeod & Zhang (2006).
}
\value{
The matrix D defined following eqn. (3) of McLeod & Zhang (2006) is computed.
}
\references{ 
McLeod, A.I. and Zhang, Y. (2006).  
Partial autocorrelation parameterization for subset autoregression.  
Journal of Time Series Analysis, 27, 599-612.
}
\author{ A.I. McLeod and Y. Zhang}
\note{ This function is used by GetFitAR. It may be used
to compute the exact loglikelihood for an AR. }

\seealso{ 
 \code{\link{GetFitARz}}, 
 \code{\link{FastLoglikelihoodAR}}, 
 \code{\link{FitAR}} 
 }
\examples{
#compute the exact concentrated loglikelihood function, (McLeod & Zhang, 2006, eq.(6)),
# for AR(p) fitted by Yule-Walker to logged lynx data
#
p<-8
CD<-ChampernowneD(log(lynx), p)
n<-length(lynx)
phi<-ar(log(lynx), order.max=p, aic=FALSE, method="yule-walker")$ar
LoglYW<-FastLoglikelihoodAR(phi,n,CD)
phi<-ar(log(lynx), order.max=p, aic=FALSE, method="burg")$ar
LoglBurg<-FastLoglikelihoodAR(phi,n,CD)
phi<-ar(log(lynx), order.max=p, aic=FALSE, method="ols")$ar
LoglOLS<-FastLoglikelihoodAR(phi,n,CD)
phi<-ar(log(lynx), order.max=p, aic=FALSE, method="mle")$ar
LoglMLE<-FastLoglikelihoodAR(phi,n,CD)
ans<-c(LoglYW,LoglBurg,LoglOLS,LoglMLE)
names(ans)<-c("YW","Burg","OLS","MLE")
ans
#compare the MLE result given by ar with that given by FitAR
FitAR(log(lynx),p)
}
\keyword{ ts }
