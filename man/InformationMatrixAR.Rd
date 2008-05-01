\name{InformationMatrixAR}
\alias{InformationMatrixAR}
\title{ Information Matrix for AR(p) }
\description{
The Fisher large-sample information matrix per observation for
the p coefficients in an AR(p) is computed.
}
\usage{
InformationMatrixAR(phi)
}
\arguments{
  \item{phi}{ vector of length p corresponding to the AR(p) coefficients }
}
\details{
The Fisher information matrix is computed as the covariance matrix
of an AR(p) process with coefficients given in the argument
\code{phi} and with unit innovation variance.
The \code{TacvfAR} function is used to compute the necessary
autocovariances.
\code{FitAR} uses \code{InformationMatrixAR} to obtain estimates
of the standard errors for the estimated parameters in the case
of the full AR(p) model. 
}
\value{
a p-by-p Toeplitz matrix, p = length(phi)
}
\references{ 
McLeod, A.I. and Zhang, Y. (2006).  
Partial autocorrelation parameterization for subset autoregression.  
Journal of Time Series Analysis, 27, 599-612.
}

\author{ A.I. McLeod and Y. Zhang}

\seealso{
\code{\link{FitAR}}, 
\code{\link{InformationMatrixARp}}, 
\code{\link{TacvfAR}}, 
\code{\link{InformationMatrixARz}} 
}

\examples{
InformationMatrixAR(c(1.8,-0.6))
}
\keyword{ ts }
