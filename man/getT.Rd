\name{getT}
\alias{getT}
\title{ t-statistic for unit root test}
\description{
Utility function used by UnitRootTest
}
\usage{
getT(ans)
}

\arguments{
  \item{ans}{ output from FitAR }
}
\value{
Value of the test statistic
}

\author{ A.I. McLeod }

\seealso{ 
\code{\link{getRho}} 
\code{\link{UnitRootTest}} 
}

\examples{
z <- cumsum(rnorm(100))
ans <- FitAR(z, p=1)
getT(ans)
}

\keyword{ ts }
\keyword{ htest }

