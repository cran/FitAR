\name{FastLoglikelihoodAR}
\alias{FastLoglikelihoodAR}
\title{ Fast Computation of the Loglikelihood Function in AR }
\description{
Computation of the loglikelihood is O(1) flops in repeated
evaluations of the loglikelihood holding the data fixed and varying
the parameters. This is useful in exact MLE estimation.
}

\usage{
FastLoglikelihoodAR(phi, n, CD)
}

\arguments{
  \item{phi}{ AR coefficients }
  \item{n}{ length of series }
  \item{CD}{ Champernowne matrix }
}
\details{
The details of this computation are described in McLeod and Zhang (2006). 
}

\value{
Loglikelihood
}

\references{ 
McLeod, A.I. and Zhang, Y. (2006).  
Partial autocorrelation parameterization for subset autoregression.  
Journal of Time Series Analysis, 27, 599-612.
}

\author{ A.I. McLeod and Y. Zhang}

\seealso{
 \code{\link{ChampernowneD}}, 
 \code{\link{LoglikelihoodAR}} 
}

\examples{
#Compute the loglikelihood using the direct method as implemented
# in LoglikelihoodAR and using the fast method
data(SeriesA)
phi<-PacfToAR(rep(0.5,10))
p<-length(phi)
z<-SeriesA-mean(SeriesA)
n<-length(z)
L1<-LoglikelihoodAR(phi, z)
cd<-ChampernowneD(z,p,MeanZero=TRUE)
L2<-FastLoglikelihoodAR(phi,n,cd)
out<-c(L1,L2)
names(out)<-c("direct","fast")
out
}
\keyword{ ts }
