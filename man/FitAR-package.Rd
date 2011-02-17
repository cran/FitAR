\name{FitAR-package}
\alias{FitAR-package}
\docType{package}
\title{
Fits AR and subset AR models and provides complete model building 
capabilities.
FitAR
}
\description{
For model estimation the main function is FitAR for which generic methods
print, summary, coef, plot and predict are implemented.
For model identification, there is a new PacfPlot for subset ARz idenfication.
Subset models may also be selected using AIC, BIC and UBIC criteria
with the function SelectModel.
SelectModel produces a S3 class object, "SelectModel",
for which their is a plot method.
The main fitting function is FitAR.
New methods and generic functions, BoxCox, Boot and sdfplot are given.
Methods for print, summary, coef, residuals, fitted and predict implemented.

}

\details{

\tabular{ll}{
Package: \tab FitAR\cr
Type: \tab Package\cr
Version: \tab 1.91\cr
Date: \tab 2011-02-17\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
LazyData: \tab yes\cr
}

To get started please see the documentation and examples given in the
functions PacfPlot, SelectModel and FitAR.

R functions for model diagnostic checking, simulation and forecasting are also available. 
The function plot provides many graphical diagnostic plots.


Model Selection: 
\code{\link{TimeSeriesPlot}}, 
\code{\link{PacfPlot}}, 
\code{\link{SelectModel}} 

Model Estimation:   
\code{\link{FitAR}}, 
\code{\link{AR1Est}}   

Model Checking:   
\code{\link{plot.FitAR}}, 
\code{\link{BoxCox}}, 
\code{\link{LBQPlot}},
\code{\link{RacfPlot}},
\code{\link{JarqueBeraTest}},

Model Applications:   
\code{\link{Boot}},
\code{\link{SimulateGaussianAR}}   

Methods Functions:  
\code{\link{coef}},
\code{\link{fitted}},
\code{\link{predict}},
\code{\link{print}},
\code{\link{summary}},
\code{\link{residuals}}  

Useful Utility Functions:  
\code{\link{Readts}}, 
\code{\link{cts}}   

New Generic and Methods Functions:
\code{\link{Boot}}, 
\code{\link{BoxCox}}, 
\code{\link{sdfplot}} 
}

\author{
A. I. McLeod and Ying Zhang

Maintainer: aimcleod@uwo.ca
}

\references{ 
McLeod, A.I. and Zhang, Y. (2006).  
Partial autocorrelation parameterization for subset autoregression.  
Journal of Time Series Analysis, 27, 599-612.

McLeod, A.I. and Zhang, Y. (2008a). Faster ARMA Maximum Likelihood Estimation, 
Computational Statistics and Data Analysis
52-4, 2166-2176.  
DOI link: \url{http://dx.doi.org/10.1016/j.csda.2007.07.020}. 

McLeod, A.I. and Zhang, Y. (2008b).
Improved Subset Autoregression: With R Package.
Journal of Statistical Software.

Changjiang Xu and A. I. McLeod (2010).
Bayesian information criterion with Bernoulli prior.
Submitted for publication.

Changjiang Xu and A. I. McLeod (2010).
Model selection using generalized information criterion.
Submitted for publication.

}

\keyword{ ts }

\keyword{ package }

\examples{
#Scripts are given below for all Figures and Tables in McLeod and Zhang (2008b).
#

#Figure 1. Plot of lynx time series using plot.ts
plot(lynx)

#Figure 2. Plot of lynx series using TimeSeriesPlot
TimeSeriesPlot(lynx, type="o", pch=16, ylab="# pelts", main="Lynx Trappings")

#Figure 3. Trellis plot for Ninemile series
graphics.off() #clear previous graphics
data(Ninemile)
print(TimeSeriesPlot(Ninemile, SubLength=200))

#Figure 4. Partial autocorrelation plot of lynx series 
graphics.off() #clear previous graphics
PacfPlot(log(lynx))

\dontrun{#takes some time for all these examples
#Figure 5. Using SelectModel to select the best subset ARz or ARp and
#          comparing BIC and UBIC subset selection.
#
graphics.off() #clear previous graphics
layout(matrix(1:4,ncol=2),respect=TRUE)
ansBICp<-SelectModel(log(lynx),lag.max=15,Criterion="BIC", ARModel="ARp", Best=3)
ansUBICp<-SelectModel(log(lynx),lag.max=15, ARModel="ARp", Best=3)
ansBICz<-SelectModel(log(lynx),lag.max=15,Criterion="BIC", ARModel="ARz", Best=3)
ansUBICz<-SelectModel(log(lynx),lag.max=15, ARModel="ARz", Best=3)
par(mfg=c(1,1))
plot(ansBICp)
par(mfg=c(2,1))
plot(ansUBICp)
par(mfg=c(1,2))
plot(ansBICz)
par(mfg=c(2,2))
plot(ansUBICz)

#Figure 6. Logged spectral density function fitted to square-root of monthly
#          sunspot series using the non-subset AR and subset ARz.
#          AIC and BIC are used for the AR while BIC and UBIC are used
#          for the ARz. Takes about 115 seconds on 3.6 GHz Pentium PC.
graphics.off() #clear previous graphics
layout(matrix(1:4,ncol=2),respect=TRUE)
z<-sqrt(sunspots)
P<-200
pAIC<-SelectModel(z, lag.max=P, ARModel="AR", Best=1, Criterion="AIC")
ARAIC<-FitAR(z, pAIC)
par(mfg=c(1,1))
sdfplot(ARAIC)
title(main="AIC Order Selection")
pBIC<-SelectModel(z, lag.max=P, ARModel="AR", Best=1, Criterion="BIC")
ARBIC<-FitAR(z, pBIC)
par(mfg=c(1,2))
sdfplot(ARBIC)
title(main="BIC Order Selection")
SunspotMonthARzBIC<-SelectModel(z,lag.max=P, ARModel="ARz", Best=1, Criterion="BIC")
ARzBIC<-FitAR(z, SunspotMonthARzBIC)
par(mfg=c(2,1))
sdfplot(ARzBIC)
title(main="BIC Subset Selection")
SunspotMonthARzUBIC<-SelectModel(z,lag.max=P, ARModel="ARz", Best=1)
ARzUBIC<-FitAR(z, SunspotMonthARzUBIC)
par(mfg=c(2,2))
sdfplot(ARzUBIC)
title(main="UBIC Subset Selection")

#Table 3.
#First part of table: AR(1) and AR(2).
#Only timings for GetFitAR and FitAR since the R function ar produces too many
#  warnings and an error message as noted in McLeod and Zhang (2008b, p.12). 
#The ar function with mle option is not recommended.
start.time<-proc.time()
set.seed(661177723)
NREP<-100 #takes about 156 sec
NREP<-10 #takes about 16 sec
ns<-c(50,100,200,500,1000)
ps<-c(1,2) #AR(p), p=1,2
tmsA<-matrix(numeric(4*length(ns)*length(ps)),ncol=4)
ICOUNT<-0
for (IP in 1:length(ps)){
p<-ps[IP]
for (ISIM in 1:length(ns)){
    ICOUNT<-ICOUNT+1
    n<-ns[ISIM]
    ptm <- proc.time()
    for (i in 1:NREP){
        phi<-PacfToAR(runif(p, min=-1, max =1))
        z<-SimulateGaussianAR(phi,n)
        phiHat<-try(GetFitAR(z,p,MeanValue=mean(z))$phiHat)
        }
    t1<-(proc.time() - ptm)[1]
#
    ptm <- proc.time()
    for (i in 1:NREP){
        phi<-PacfToAR(runif(p, min=-1, max =1))
        z<-SimulateGaussianAR(phi,n)
        phiHat<-try(FitAR(z,p,MeanMLEQ=TRUE)$phiHat)
        }
    t2<-(proc.time() - ptm)[1]
#
    ptm <- proc.time()
    for (i in 1:NREP){
        phi<-PacfToAR(runif(p, min=-1, max =1))
        z<-SimulateGaussianAR(phi,n)
        #uncomment this line and next two lines for ar timings -- expect lots of
        #   warnings and an error message!!
        #phiHat<-try(ar(z,aic=FALSE,order.max=p,method="mle")$ar)
        #delete this line and the next one
        phiHat<-NA
    }
    #uncomment this line for ar timings
         #t3<-(proc.time() - ptm)[1]
         t3<-NA #delete this line for ar timings

        tmsA[ICOUNT,]<-c(n,t1,t2,t3)
  } 
} 
rnames<-c(rep("AR(1)", length(ns)),rep("AR(2)", length(ns)) )
cnames<-c("n", "GetFitAR", "FitAR", "ar")  
dimnames(tmsA)<-list(rnames,cnames)
tmsA[,-1]<-round(tmsA[,-1]/NREP,2)
end.time<-proc.time()
total.time<-(end.time-start.time)[1]

#Second part of table: AR(20) and AR(40).
#NOTE: ar is not recommended with method="mle" produces numerous warnings
#      and also takes a long time!
        start.time<-proc.time()
        set.seed(661177723)
        NREP<-100 #takes 7.5 hours
        NREP<-10 #takes 45 minutes
        ns<-c(1000,2000,5000)
        ps<-c(20,40)
        tmsB<-matrix(numeric(4*length(ns)*length(ps)),ncol=4)
        ICOUNT<-0
        for (IP in 1:length(ps)){
        p<-ps[IP]
        phi<-PacfToAR(0.8/(1:p))
        for (ISIM in 1:length(ns)){
            ICOUNT<-ICOUNT+1
            n<-ns[ISIM]
            ptm <- proc.time()
            for (i in 1:NREP){
                z<-SimulateGaussianAR(phi,n)
                phiHat<-try(GetFitAR(z,p,MeanValue=mean(z))$phiHat)
                }
            t1<-(proc.time() - ptm)[1]
            ptm <- proc.time()
            for (i in 1:NREP){
                z<-SimulateGaussianAR(phi,n)
                phiHat<-try(FitAR(z,p,MeanMLEQ=TRUE)$phiHat)
                }
            t2<-(proc.time() - ptm)[1]
            ptm <- proc.time()
            for (i in 1:NREP){
                z<-SimulateGaussianAR(phi,n)
                phiHat<-try(ar(z,aic=FALSE,order.max=p,method="mle")$ar)
                }
            t3<-(proc.time() - ptm)[1]
            tmsB[ICOUNT,]<-c(n,t1,t2,t3)
        } 
        }   
        rnames<-c( rep("AR(20)", length(ns)), rep("AR(40)", length(ns)) )
        cnames<-c("n", "GetFitAR", "FitAR", "ar")  
        dimnames(tmsB)<-list(rnames,cnames)
        tmsB[,-1] <- round(tmsB[,-1]/NREP,2)
        end.time<-proc.time()
        total.time<-(end.time-start.time)[1]

#Figure 7. Comparing Box-Cox analyses using FitAR and MASS
library(MASS)
graphics.off() #clear previous graphics
layout(matrix(c(1,2,1,2),ncol=2))
pvec<-c(1,2,4,10,11)
out<-FitAR(lynx, ARModel="ARp", pvec)
BoxCox(out)
PMAX<-max(pvec)
Xy <- embed(lynx, PMAX + 1)
y <- Xy[, 1]
X <- (Xy[, -1])[, pvec] #pvec != 1
outlm<-lm(y~X)
boxcox(outlm,lambda=seq(0.0,0.6,0.05))

#Figure 8
graphics.off() #clear previous graphics
BoxCox(AirPassengers) #takes about 30 sec

#Figure 9
graphics.off() #clear previous graphics
data(rivers)
BoxCox(rivers)
title(sub="Length of 141 North American Rivers")

#Figure 10
graphics.off() #clear previous graphics
data(USTobacco)
TimeSeriesPlot(USTobacco, aspect=1)

#Figure 11
graphics.off() #clear previous graphics
data(USTobacco)
outUST<-arima(USTobacco, c(0,1,1))
BoxCox(outUST)

#Figure 12. Basic diagnostic plots for ARp fitted to the log lynx series
graphics.off() #clear previous graphics
out<-FitAR(log(lynx), ARModel="ARp", c(1,2,4,10,11))
plot(out, terse=TRUE)

#Figure 13. RSF plot for ARp fitted to log lynx series
graphics.off() #clear previous graphics
out<-FitAR(log(lynx), ARModel="ARp", c(1,2,4,10,11))
rfs(out)

#Table 6. Comparison of bootstrap and large-sample sd
#Use bootstrap to compute standard errors of parameters
#takes about 34 seconds on a 3.6 GHz PC
ptm <- proc.time() #user time
set.seed(2491781) #for reproducibility
R<-100  #number of bootstrap iterations
p<-c(1,2,4,7,10,11)
ans<-FitAR(log(lynx),p)
out<-Boot(ans, R)
fn<-function(z) FitAR(z,p)$zetaHat
sdBoot<-sqrt(diag(var(t(apply(out,fn,MARGIN=2)))))
sdLargeSample<-coef(ans)[,2][1:6]
sd<-matrix(c(sdBoot,sdLargeSample),ncol=2)
dimnames(sd)<-list(names(sdLargeSample),c("Bootstrap","LargeSample"))
ptm<-(proc.time()-ptm)[1]
sd

}


}
