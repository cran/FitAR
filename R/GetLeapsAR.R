GetLeapsAR <-
function(z, lag.max=15, Criterion="UBIC", Best=3, Candidates=5, G=0.5, Q=0.25, level=0.99, ExactQ=FALSE){
stopifnot(length(z)>0, length(z)>lag.max, lag.max>1, Best>0)
is.wholenumber <-
    function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
stopifnot(is.wholenumber(lag.max))
method<-Criterion
if (is.na(pmatch(method,c("UBIC","AIC","BIC","EBIC","BICq","GIC"))))
    method<-"UBIC"
if (ExactQ){
    stop("Sorry this option is not available yet")
    M<-2^lag.max
    logL <- numeric(M)
    logL[1] <- GetFitARpMLE(z, p=0)$loglikelihood
    for (i in 1:(M-1)){
        ind<-as.logical(rev(to.binary(i, lag.max)))
        pvec <- (1:lag.max)[ind]
        out<-GetFitARpMLE(z, p=pvec)
        logL[i+1] <- out$loglikelihood
    }
}
#GIC/BICq treated as a special case
if (method=="GIC"){
    pvec <- 1:lag.max
    n <- length(z)-lag.max
    ind <- (lag.max+1):length(z)
    y<-z[ind]
    X<-matrix(rep(0,n*lag.max), nrow=n, ncol=lag.max)
    for (i in 1:lag.max)
        X[,i] <- z[ind-pvec[i]]
    outLeaps <- leaps(y=y,x=X,nbest=1,method="r2",strictly.compatible=FALSE)
    k <- outLeaps$size #k=(2,...,lag.max+1) since k=1 is mean only.
# approximate likelihood approach
    TotSS <- sum((y-mean(y))^2)
    RSS <- TotSS*(1-outLeaps$r2)
    LogL <- (-n/2)*log(c(TotSS/n, RSS/n))
    ans <- BICqLL(logL=LogL, n=n, level=level)
    kHat <- (ans$khat)[[1]]
    kHat <- kHat - 1 #null model corresponds to kHat=0
    pvec <-0
    if (kHat > 0)
        pvec <- (1:lag.max)[(outLeaps$which)[kHat,]]
    return(pvec)
}
pvec <- 1:lag.max
LagRange<-1:lag.max
n <- length(z)-lag.max
ind <- (lag.max+1):length(z)
y<-z[ind]
X<-matrix(rep(0,n*lag.max), nrow=n, ncol=lag.max)
for (i in 1:lag.max)
    X[,i] <- z[ind-pvec[i]]
outLeaps <- leaps(y=y,x=X,nbest=1,method="r2", strictly.compatible=FALSE)
k <- outLeaps$size
#
#this defines an approximate likelihood approach
TotSS <- sum((y-mean(y))^2)
RSS <- TotSS*(1-outLeaps$r2)
LogL <- (-n/2)*log(RSS/n)
if (method=="AIC")
    ic<- -2*LogL + 2*k
if (method=="BIC")
    ic<- -2*LogL + log(n)*k
if (method=="UBIC")
    ic<- -2*LogL + log(n)*k + 2*lchoose(lag.max+1, k)
if (method=="EBIC")
    ic<- -2*LogL + log(n)*(1+LagRange)+2*G*lchoose(lag.max, LagRange)
if (method=="BICq")
    ic<- -2*LogL + log(n)*(1+LagRange)-2*(LagRange*log(Q)+(lag.max+1-LagRange)*log(1-Q)) 
indBest<-order(ic)
#extra step needed because leaps does not include null model
#Very important: refit with exact MLE
LogL<-numeric(Candidates+1)
#LogL[1]<-GetFitARpLS(z, 0)$loglikelihood #null model included here
LogL <-GetFitAR(z, 0)$loglikelihood #null model included here
for (i in 1:Candidates)
       LogL[i+1]<-GetFitAR(z,pvec[outLeaps$which[indBest[i],]])$loglikelihood
       #LogL[i+1]<-GetFitARpLS(z,pvec[outLeaps$which[indBest[i],]])$loglikelihood
k<-c(1,k[indBest[1:Candidates]])
if (method=="AIC")
    ic<- -2*LogL + 2*k
if (method=="BIC")
    ic<- -2*LogL + log(n)*k
if (method=="UBIC")
    ic<- -2*LogL + log(n)*k + 2*lchoose(lag.max+1, k)
if (method=="BICq")
    ic<- -2*LogL + log(n)*k-2*(k*log(Q)+(lag.max+1-k)*log(1-Q)) 
icBest<-order(ic)[1:Best] #best models exact
m<-as.list(numeric(Best))
for (i in 1:Best){
    ind<-icBest[i]
    if (indBest[ind] ==1)
        p<-0
    else
        p<-pvec[outLeaps$which[indBest[ind-1],]]
    if (method == "AIC")
        m[[i]] <-list(p=p, AIC=ic[ind])
    if (method == "BIC")
        m[[i]] <-list(p=p, BIC=ic[ind])
    if (method == "UBIC")
        m[[i]] <-list(p=p, UBIC=ic[ind])
    if (method == "BICq")
        m[[i]] <-list(p=p, BICq=ic[ind])
    }
class(m)<-"Selectmodel"
attr(m,"model")<-"ARp"
if (Best>1)
    ans<-m
else {
    ans <- m[[1]]$p
	if (length(ans)==0)
		ans<-0
	}
ans
}
