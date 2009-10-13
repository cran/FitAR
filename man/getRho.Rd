\name{getRho}
\alias{getRho}
\title{ Normalized rho unit root test statistic }
\description{
Utility function used by UnitRootTest
}
\usage{
getRho(ans)
}
\arguments{
  \item{ans}{ output from FitAR }
}
\value{
Value of the test statistic
}
\author{ A.I. McLeod }

\seealso{ 
\code{\link{getT}} 
\code{\link{UnitRootTest}} 
}
\examples{
z <- cumsum(rnorm(100))
ans <- FitAR(z, p=1)
getRho(ans)
}

\keyword{ ts }
\keyword{ htest }
