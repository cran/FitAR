\name{ARSdf}
\alias{ARSdf}
\title{ Autoregressive Spectral Density Function }
\description{
Spectral density function of AR(p) is computed.
}
\usage{
ARSdf(phi, pFFT = 8)
}
\arguments{
  \item{phi}{AR Coefficient vector }
  \item{pFFT}{ FFT with \eqn{2^pFFF} frequencies, default 8 }
}
\details{
The Fast Fourier Transform (FFT) is used to compute the spectral density function.
}
\value{
A vector of the density function values, \eqn{(f(1),..., f(2^pFFF))}
}
\author{ A.I. McLeod and Y. Zhang}

\seealso{ 
\code{\link{spectrum}}, 
\code{\link{spec.pgram}}, 
\code{\link{spec.ar}} 
}
\examples{
ARSdf(0.8)
ARSdf(c(0.1,0.2))
}
\keyword{ ts }
