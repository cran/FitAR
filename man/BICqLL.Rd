\name{BICqLL}
\alias{BICqLL}
\title{
Select best model using BICq
}

\description{
Given the loglikelihoods for a set of models arranged in ascending order of size,
the best model is selected using the BICq criterion for a specified size.
}
\usage{
BICqLL(logL, n, level = 0.99, mSize = 1:length(logL), mComplex = function(k) k)
}
\arguments{
  \item{logL}{vector of loglikelihoods}
  \item{n}{sample size}
  \item{level}{probability}
  \item{mSize}{model sizes}
  \item{mComplex}{a complexity function}
}

\details{See reference}

\value{
\item{khat}{dataframe with columns: k, a.1, a.2 q.1, q.2, level=level,
where k is the optimal model, (a.1,a.2) is the interval for alpha in the
GIC, (q.1, q.2) is the interval for q and level is the probability. 
Each row corresponds to an entry in 'level'. }
\item{table }{This table indicates which models can be selected
for some values of alpha or q.}
}

\references{
Changjiang Xu and A. Ian McLeod (2010).
Bayesian information criterion with Bernoulli prior. 
Submitted for publication.
}

\author{Changjiang Xu and A. Ian McLeod}

\note{
AIC corresponds to setting level=0.84.
BIC corresponds to setting level=pchisq(log(n), 1).
So for n=100, 1000; BIC=0.96, 0.97
}

\seealso{
\code{\link{SelectModel}}
}
\examples{
#Example 1.
#AR(p) Order selection for 'lynx' series
z <- log(lynx)
n <- length(z)
lag.max <- 20
zta<-ARToPacf(ar.burg(z,aic=FALSE,order.max=lag.max)$ar)
LagsEntering<-1:lag.max
LLapprox <- (-n)*log(cumprod(1-zta[LagsEntering]^2))
ans<-BICqLL(logL=LLapprox, n=n, level=c(0.9, 0.95, 0.99))
ans$khat
ans$table
#if we just want the best model for level=0.99 then,
(BICqLL(logL=LLapprox, n=n, level=0.99)$khat)[[1]]
#aic for comparison
aic<-(-2*LLapprox)+2*LagsEntering
which.min(aic)
plot(LagsEntering, aic)
#
#Example 2. AR(p) Order Selection
#white noise. We do NumRep simulations and 
#  count the number of overfit models.
set.seed(231789) #make reproducible
n <- 100
lag.max <- 30
LagsEntering<-0:lag.max
NumRep<-25
level<-c(0.99, 0.95, 0.9)
k<-numeric(length(level))
for (i in 1:NumRep){
    z <- rnorm(n)
    zta<-ARToPacf(ar.burg(z,aic=FALSE,order.max=lag.max)$ar)
    LLapprox <- c(0,(-n)*log(cumprod(1-zta[LagsEntering]^2)))
    k<-k+as.numeric(0<(BICqLL(logL=LLapprox, n=n, level=level,mSize=LagsEntering)$khat)[,1])
    }
ans<-k
names(ans)<-level
ans
#
#Example 3. AR(p) best subset. ARz Family.
z <- log(lynx)
n <- length(z)
lag.max <- 20
zta <- ARToPacf(ar.burg(z,aic=FALSE,order.max=lag.max)$ar)
LagsEntering <- order(abs(zta),decreasing=TRUE)
LLapprox <- c(0, -n*log(cumprod(1-zta[LagsEntering]^2)))
kHat <- (BICqLL(logL=LLapprox, n=n, level=0.99)$khat)[[1]]
pvec<-LagsEntering[1:kHat]
pvec
#pvec above shows the lags in order of importance
#
#Example 4. AR(p) best subset. ARp Family.
#could also try z <- sunspot.year
z <- log(lynx)
lag.max <- 15
pvec <- 1:lag.max
n <- length(z)-lag.max
ind <- (lag.max+1):length(z)
y<-z[ind]
X<-matrix(rep(0,n*lag.max), nrow=n, ncol=lag.max)
for (i in 1:lag.max)
    X[,i] <- z[ind-pvec[i]]
outLeaps <- leaps(y=y,x=X,nbest=1,method="r2",strictly.compatible=FALSE)
# approximate likelihood approach
TotSS <- sum((y-mean(y))^2)
RSS <- TotSS*(1-outLeaps$r2)
LogL <- (-n/2)*log(c(TotSS/n, RSS/n))#null model included
ans<-BICqLL(logL=LogL, n=n, level=0.99)
kHat <- (ans$khat)[[1]]-1 #kHat=0 is null model
pvec <-0
if (kHat > 0)
    pvec <- (1:lag.max)[(outLeaps$which)[kHat,]]
pvec

}
\keyword{ ts }
