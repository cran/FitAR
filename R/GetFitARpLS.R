`GetFitARpLS` <-
function(z, pvec){
stopifnot(length(z)>0, length(z)>max(pvec), length(pvec)>0)
PMAX<-max(pvec)
if (PMAX == 0){
    phiHat<-numeric(0)
    constantTerm<-mean(z)
    res<-z-constantTerm
    yX<-matrix(z,ncol=1)
    Iq<-TRUE
    LL<-LoglikelihoodAR(0,z, MeanValue=constantTerm)
    }
else {
    Xy <- embed(z, PMAX+1)
    y <- Xy[,1]
    if (length(pvec)==1 && pvec == 1)
        X <- matrix(Xy[,-1], ncol=1)
    else    
        X <- (Xy[,-1])[,pvec]
    yX<-cbind(y,X)
    ans<-lsfit(X,y)
    res<-ans$residuals
    estimates<-lsfit(X,y)$coefficients
    betaHat<-as.vector(estimates[-1])
    constantTerm<-as.vector(estimates[1])
    phiHat<-numeric(PMAX)
    phiHat[pvec]<-betaHat
    Iq<-InvertibleQ(phiHat)
    if (Iq)
        LL<-LoglikelihoodAR(phiHat,z, MeanValue=mean(z))
    else
        LL<- -1e30
    }
list(loglikelihood=LL, phiHat=phiHat, constantTerm=constantTerm,res=res,pvec=pvec,
     InvertibleQ=Iq,  yX=yX)
}
