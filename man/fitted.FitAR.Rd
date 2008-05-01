\name{fitted.FitAR}
\alias{fitted.FitAR}
\title{ Fitted Values from "FitAR" Object }
\description{
Method function, extracts fitted values from \code{FitAR} object.
}
\usage{
\method{fitted}{FitAR}(object, ...)
}
\arguments{
  \item{object}{ object of class "FitAR" }
  \item{...}{ optional arguments }
}
\value{
Vector of fitted values
}
\author{ A.I. McLeod and Y. Zhang}

\seealso{ 
\code{\link{FitAR}}
}
\examples{
data(SeriesA)
out<-FitAR(SeriesA, c(1,2,6,7))
fitted(out)
}
\keyword{ ts }

