\name{AR1Est}
\alias{AR1Est}
\title{ Exact MLE Mean-Zero AR(1) }
\description{
This function is used by GetFitAR in the AR(1) case. It is a fast exact solution using
the root of a cubic equation.
}
\usage{
AR1Est(z, MeanValue = 0)
}
\arguments{
  \item{z}{ time series or vector }
  \item{MeanValue}{ known mean }
}
\details{
The exact MLE for mean-zero AR(1) satisfies a cubic equation.  The solution of this
equation for the MLE given by Zhang (2002) is used.  This approach
is more reliable as well as faster than the usual approach to the
exact MLE using a numerical optimization technique which can occasionally
have convergence problems.
}
\value{
MLE for the parameter
}

\references{
Zhang, Y. (2002). 
Topics in Autoregression,
Ph.D. Thesis, University of Western Ontario.
}

\author{ A.I. McLeod and Y. Zhang }

\seealso{ \code{\link{GetFitARz}} }

\examples{
AR1Est(lynx-mean(lynx))
}

\keyword{ ts }
