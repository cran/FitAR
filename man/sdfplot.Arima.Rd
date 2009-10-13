\name{sdfplot.Arima}
\alias{sdfplot.Arima}
\title{ Spectral Density of Fitted ARIMA Model }
\description{
Method for class "Arima" for sdfplot.  
}

\usage{
\method{sdfplot}{Arima}(obj, ...)
}
\arguments{
  \item{obj}{ object class "Arima" }
  \item{\dots}{ optional arguments }
}

\value{
Plot is produced using plot.
Matrix with 2 columns containing the frequencies and spectral density is returned invisibly.
}

\author{ A.I. McLeod }

\seealso{ 
\code{\link{ARSdf}},
\code{\link{sdfplot}},
\code{\link{sdfplot.FitAR}},
\code{\link{sdfplot.ts}}
 }

\examples{
data(SeriesA)
sdfplot(SeriesA,c(1,0,1))
}

\keyword{ ts }
