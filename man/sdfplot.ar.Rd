\name{sdfplot.ar}
\alias{sdfplot.ar}
\title{ Autoregressive Spectral Density Estimation for "ar" }
\description{
Method for class "ar" for sdfplot. 
}
\usage{
\method{sdfplot}{ar}(obj, ...)
}
\arguments{
  \item{obj}{class "ar" object, output from \code{ar}  }
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
#Fit AR(p) using AIC model selection and Burg estimates. Plot spectral density
#estimate
ans<-ar(lynx, lag.max=20)
sdfplot(ans)
}
\keyword{ ts}
