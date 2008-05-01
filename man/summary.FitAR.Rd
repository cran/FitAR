\name{summary.FitAR}
\alias{summary.FitAR}
\title{ Summary Method for "FitAR" Object }
\description{
summary for "FitAR" object.
}
\usage{
\method{summary}{FitAR}(object, ...)
}
\arguments{
  \item{object}{ "FitAR" object }
  \item{...}{ optional arguments }
}
\value{
A printed summary is given
}
\author{ A.I. McLeod }


\seealso{ 
\code{\link{print.FitAR}}, 
\code{\link{FitAR}} 
}
\examples{
data(SeriesA)
out<-FitAR(SeriesA, c(1,2,6,7))
summary(out)
}
\keyword{ ts }
