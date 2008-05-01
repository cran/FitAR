\name{sdfplot.numeric}
\alias{sdfplot.numeric}
\title{ Autoregressive Spectral Density Estimation for "numeric" }
\description{
 Method function for vectors, class "numeric"
}

\usage{
\method{sdfplot}{numeric}(obj, ...)
}


\arguments{
  \item{obj}{ object, class"numeric", a vector }
  \item{\dots}{optional arguments }
}
\value{
Plot produced as side-effect. No output.
}
\author{A.I. McLeod }

\seealso{ \code{\link{sdfplot}} }
\examples{
sdfplot(lynx)
}
\keyword{ ts }

