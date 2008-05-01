\name{sdfplot.ts}
\alias{sdfplot.ts}
\title{ Autoregressive Spectral Density Estimation for "ts" Object }
\description{
  Methods function for "ts". 
}
\usage{
\method{sdfplot}{ts}(obj, ...)
}
\arguments{
  \item{obj}{ object, class"ts" }
  \item{\dots}{optional arguments }
}
\value{
Plot produced as side-effect. No output.
}
\author{A.I. McLeod }

\seealso{ \code{\link{sdfplot}} }
\examples{
data(SeriesA)
sdfplot(SeriesA)
}
\keyword{ ts }
