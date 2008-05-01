`GetLeapsAR` <-
function(z, lag.max=15, Criterion="UBIC", Best=3, Candidates=5){
stopifnot(length(z)>0, length(z)>lag.max, lag.max>1, Best>0)
method<-Criterion
if (is.na(pmatch(method,c("UBIC","AIC","BIC"))))
    method<-"UBIC"
pvec <- 1:lag.max
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
indBest<-order(ic)
#extra step needed because leaps does not include null model
LogL<-numeric(Candidates+1)
LogL[1]<-GetFitARpLS(z, 0)$loglikelihood #null model included here
for (i in 1:Candidates)
       LogL[i+1]<-GetFitARpLS(z,pvec[outLeaps$which[indBest[i],]])$loglikelihood
k<-c(1,k[indBest[1:Candidates]])
if (method=="AIC")
    ic<- -2*LogL + 2*k
if (method=="BIC")
    ic<- -2*LogL + log(n)*k
if (method=="UBIC")
    ic<- -2*LogL + log(n)*k + 2*lchoose(lag.max+1, k)
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
    }
class(m)<-"Selectmodel"
attr(m,"model")<-"ARp"
if (Best>1)
    m
else
    m[[1]]$p
}

