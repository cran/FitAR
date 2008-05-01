\name{Boot}
\alias{Boot}
\title{Generic Bootstrap Function  }
\description{
Generic function to bootstrap a fitted model.
}
\usage{
Boot(obj, R=1, ...)
}
\arguments{
  \item{obj}{fitted object }
  \item{R}{number of bootstrap replicates }
  \item{...}{optional arguments }
}
\details{
At present, the only function implemented is
\code{\link{Boot.FitAR}}.
}
\value{
Parametric bootstrap simulation
}
\author{ A.I. McLeod and Y. Zhang}


\seealso{ \code{\link{Boot.FitAR}} }
\examples{
data(SeriesA)
out<-FitAR(SeriesA, c(1,2,7), ARModel="ARp")
Boot(out)

}
\keyword{ ts }
