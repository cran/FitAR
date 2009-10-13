`SelectModel` <-
function(z, lag.max=15, ARModel=c("AR","ARz","ARp"), Criterion="default", Best=3, Candidates=5, G=0.5, Q=0.25){
stopifnot(length(z)>0, length(z)>lag.max, lag.max>1, Best>0, Candidates>0)
BestCandidates<-Candidates
IsValidCriterionQ <- Criterion %in% c("default", "AIC", "BIC", "UBIC", "EBIC", "QBIC")
if (!IsValidCriterionQ)
    stop("Criterion = ", Criterion, " not known.")
ARModel <- match.arg(ARModel)
if (Best > BestCandidates)
    BestCandidates<-Best
if (ARModel=="ARp") #subset ARp
    return(GetLeapsAR(z, lag.max=lag.max, Criterion=Criterion, Best=Best, BestCandidates))
if (ARModel=="ARz")
    SubsetQ <- TRUE
else
    SubsetQ <- FALSE
method<-Criterion
if (Criterion == "default")
    if (SubsetQ)
        method <- "QBIC"
    else
        method <- "BIC"
if (!SubsetQ && Criterion=="UBIC")
    method <- "BIC"
zta<-ARToPacf(ar.burg(z,aic=FALSE,order.max=lag.max)$ar)
n<-length(z)
LagRange<-1:lag.max
if (method=="UBIC")
    penalty<-log(n)*(1+LagRange)+2*lchoose(lag.max, LagRange)
if (method=="EBIC")
    penalty<-log(n)*(1+LagRange)+2*G*lchoose(lag.max, LagRange)
if (method=="QBIC")
    penalty<-log(n)*(1+LagRange)-2*(LagRange*log(Q)+(lag.max+1-LagRange)*log(1-Q)) 
if (method=="BIC")
    penalty<-log(n)*(1+LagRange)
if (method=="AIC")
    penalty<-2*(1+LagRange)
if (SubsetQ)
    LagsEntering<-order(abs(zta),decreasing=TRUE)
else  
    LagsEntering<-1:lag.max
LLapprox <- n*log(cumprod(1-zta[LagsEntering]^2))
AnIC <- LLapprox + penalty
IndCandidates<-order(AnIC)[1:BestCandidates]
AnICexact<-numeric(BestCandidates+1)
if (SubsetQ){ #subset. AR model subset selection.
    m<-as.list(numeric(BestCandidates+1))
    for (isub in 1:BestCandidates){
        ModelLags<-sort(LagsEntering[1:IndCandidates[isub]])
        LL<-GetFitAR(z-mean(z), ModelLags)$loglikelihood
        k<-length(ModelLags)+1 #mean is included and k>=2 here
        if (method=="UBIC") {
            UBIC <- -2*LL + log(n)*k + 2*lchoose(lag.max+1, k) 
            AnICexact[isub]<-UBIC
            m[[isub]] <- list(p=ModelLags, UBIC=UBIC)
            }
        if (method=="EBIC") {
            EBIC <- -2*LL + log(n)*k + 2*G*lchoose(lag.max+1, k) 
            AnICexact[isub]<-EBIC
            m[[isub]] <- list(p=ModelLags, EBIC=EBIC)
            }
        if (method=="QBIC") {
            QBIC <- -2*LL + log(n)*k -2*(k*log(Q)+(lag.max+1-k)*log(1-Q))
            AnICexact[isub]<-QBIC
            m[[isub]] <- list(p=ModelLags, QBIC=QBIC)
            }
        if (method=="AIC"){
            AIC <- -2*LL+2*k
            AnICexact[isub]<-AIC
            m[[isub]] <- list(p=ModelLags, AIC=AIC)
            }
        if (method=="BIC") {
            BIC <- -2*LL+log(n)*k
            AnICexact[isub]<-BIC
            m[[isub]] <- list(p=ModelLags, BIC=BIC)
            }
    }
        #null model
        LL<-GetFitAR(z-mean(z), 0)$loglikelihood
        if (method=="UBIC") {#parameters=1, just mean
            UBIC <- -2*LL + log(n) + 2*lchoose(lag.max+1, 1) 
            AnICexact[BestCandidates+1]<-UBIC
            m[[BestCandidates+1]] <- list(p=0, UBIC=UBIC)
        }
       if (method=="EBIC") {
            EBIC <- -2*LL + log(n) + 2*G*lchoose(lag.max+1, 1) 
            AnICexact[BestCandidates+1]<-EBIC
            m[[BestCandidates+1]] <- list(p=0, EBIC=EBIC)
        }
        if (method=="QBIC") {
            QBIC <- -2*LL + log(n) -2*(log(Q)+(lag.max)*log(1-Q)) 
            AnICexact[BestCandidates+1]<-QBIC
            m[[BestCandidates+1]] <- list(p=0, QBIC=QBIC)
        }
        if (method=="AIC"){
            AIC <- -2*LL+2
            AnICexact[BestCandidates+1] <- AIC
            m[[BestCandidates+1]]<-list(p=0,AIC=AIC)
            }
        if (method=="BIC") {
            BIC <- -2*LL+log(n)
            AnICexact[BestCandidates+1] <- BIC
            m[[BestCandidates+1]]<-list(p=0,BIC=BIC)
            }
        #final model select based on exact likelihood
        i<-order(AnICexact)
        m<-m[i]
        m<-m[1:Best] 
        attr(m, "model")<-ARModel              
    }
else  { #non-subset. AR model order selection.
    AnICexact<-numeric(BestCandidates)
    AnICApprox<-numeric(BestCandidates)
    for (i in 1:BestCandidates){
        p<-LagsEntering[IndCandidates[i]]-1
        AnICApprox[i]<-AnIC[p+1]
        ans<-GetFitAR(z-mean(z), 0:p)
        LL<-ans$loglikelihood
#mean included in all models
        if (method=="AIC")
                penalty<-2*(p+1)
        else
                penalty<-log(n)*(p+1)
        AnICexact[i]<- -2*LL+penalty
    }
    m<-c(LagsEntering[IndCandidates]-1,AnICexact,AnICApprox)
    m<-matrix(m,ncol=3)
    m<-m[order(AnICexact),]
    m<-m[1:Best,]
    if (Best > 1)
        if (method=="AIC")
            dimnames(m)<-c(list(1:Best),list(c("p", "AIC-Exact", "AIC-Approx")))
        else
            dimnames(m)<-c(list(1:Best),list(c("p", "BIC-Exact", "BIC-Approx")))
    }
if (SubsetQ) class(m)<-"Selectmodel"
if (Best > 1)
    m
else 
    if (is.list(m))
        m[[1]]$p
    else
        as.vector(m[1])
}
