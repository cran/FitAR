`LjungBoxTest` <-
function(res, k=0, lag.max=30, StartLag=10, SquaredQ=FALSE){
n<-length(res)
L0<-StartLag
if (lag.max<=k)
    stop("invalid parameter setting: lag.max must be > k")
if (StartLag>lag.max || StartLag<=k) 
    L0<-lag.max
if (SquaredQ) {
    z<-(res-mean(res))^2
    kpar<-0
    }
else {
    z<-res
    kpar<-k
}
ra<-(acf(z, lag.max=lag.max, plot=FALSE)$acf)[-1]
lags<-L0:lag.max
QQ<-n*(n+2)*cumsum((ra^2)/(n-(1:lag.max)))[lags]
pv<-1-pchisq(QQ,lags-kpar)
QQ<-round(QQ,2)
a<-matrix(c(lags,QQ,pv),ncol=3)
dimnames(a)<-list(rep("",length(QQ)),c("m","Qm", "pvalue"))
a
}

