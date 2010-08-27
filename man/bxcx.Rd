\name{bxcx}
\alias{bxcx}
\title{ Box-Cox Transformation and its Inverse }
\description{
Box-Cox or power transformation or its inverse.
For \eqn{lambda \ne 0}{lambda!=0}, the Box-Cox transformation of x is
\eqn{(x^\lambda - 1)/\lambda}{(x^lambda-1)/lambda},
whereas the regular power transformation 
is simply \eqn{x^\lambda}{x^lambda}.  
When \eqn{\lambda=0}{lambda=0}, it is log in both cases.
The inverse of the Box-Cox and the power transform can
also be obtained.
}

\usage{
bxcx(x, lambda, InverseQ = FALSE, type = "BoxCox")
}

\arguments{
  \item{x}{ a vector or time series }
  \item{lambda}{ power transformation parameter }
  \item{InverseQ}{ if TRUE, the inverse transformation is done}
  \item{type}{ either "BoxCox" or "power"}
}


\value{
A vector or time series of the transformed data
}

\references{ 
Box, G. E. P. and Cox, D. R. (1964) An analysis of transformations. 
Journal of Royal Statistical Society, Series B, vol. 26, pp. 211-246.
 }
\author{ A.I. McLeod }


\seealso{ \code{\link{BoxCox}} }

\examples{
#lambda=0.5
z<-AirPassengers; lambda<-0.5
y<-bxcx(z, lambda)
z2<-bxcx(y, lambda, InverseQ=TRUE)
sum(abs(z2-z))
#
z<-AirPassengers; lambda<-0.0
y<-bxcx(z, lambda)
z2<-bxcx(y, lambda, InverseQ=TRUE)
sum(abs(z2-z))
}
\keyword{ ts }
