\name{Get1G}
\alias{Get1G}
\title{ Internal Utility Function: BLUE Mean }
\description{
  This function is not normally used directly by the user.
  It is used in the exact mle for mean.
}
\usage{
Get1G(phi, n)
}
\arguments{
  \item{phi}{ a vector of AR coefficients }
  \item{n}{ length of series }
}

\value{
A vector used in the mle computation of the mean.
}
\author{ A.I. McLeod }

\seealso{ \code{\link{GetARMeanMLE}} }
\examples{
#Simulate an AR(2) and compute the exact mle for mean
set.seed(7771111)
n<-50
phi<-c(1.8,-0.9)
z<-SimulateGaussianAR(phi, n)
g1<-Get1G(phi, length(z))
sum(g1*z)/sum(g1)
#sample mean
mean(z)
#more directly with getArMu
GetARMeanMLE(z,phi)
}
\keyword{ ts }
