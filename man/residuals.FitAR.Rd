\name{residuals.FitAR}
\alias{residuals.FitAR}
\title{ Extract Residuals from "FitAR" Object }
\description{
Method function.
}
\usage{
\method{residuals}{FitAR}(object, \dots)
}
\arguments{
  \item{object}{ object of class "FitAR" }
  \item{\dots}{ optional arguments }
}
\value{
Vector of residuals
}
\author{ A.I. McLeod }

\seealso{ 
\code{\link{FitAR}},
}
\examples{
out<-FitAR(SeriesA, c(1,2,6,7))
resid(out)
}
\keyword{ ts }
