\name{SiddiquiMatrix}
\alias{SiddiquiMatrix}
\title{ Covariance Matrix of MLE Parameters in an AR(p) }
\description{
A direct method of computing the inverse of the covariance
matrix of p successive observations in an AR(p) with unit
innovation variance given by Siddiqui (1958) is implemented.
This matrix, divided by n = length of series, 
is the covariance matrix for the MLE estimates
in a regular AR(p).
}
\usage{
SiddiquiMatrix(phi)
}
\arguments{
  \item{phi}{ coefficients in a regular AR(p) }
}
\value{
Matrix, covariance matrix of MLE estimates
}

\references{ 
Siddiqui, M.M. (1958)
On the inversion of the sample covariance matrix in a 
stationary autoregressive process.
Annals of Mathematical Statistics 29, 585-588.

Pagano, M. (1973), When is an autoregressive scheme stationary? 
Communications in Statistics A 1, 533-544.
}
\author{ A.I. McLeod }
\note{ 
No check on whether the parameters are in the stationary region is done.
It has been shown a necessary and sufficient condition for the
parameters to be in the stationary region is that this matrix should
be positive-definite (Pagano, 1973).
But computationally it is probably better to test for stationarity
by using \code{\link{ARToPacf}} to transform to the PACF and then
check that the absolute value of all partial autocorrelations are
less than 1.
}

\seealso{ \code{\link{FitAR}} }
\examples{
#compute the inverse directly and by Siddiqui's method and compare:
phi<-PacfToAR(rep(0.8,5))
A<-SiddiquiMatrix(phi)
B<-solve(toeplitz(TacvfAR(phi, lag.max=length(phi)-1)))
max(abs(A-B))

}
\keyword{ ts }
