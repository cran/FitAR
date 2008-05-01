\name{sdfplot.FitAR}
\alias{sdfplot.FitAR}
\title{ Autoregressive Spectral Density Estimation for "FitAR" }
\description{
Methods function for sdfplot. Autoregressive spectral density function
estimation using the result output from FitAR.
}
\usage{
\method{sdfplot}{FitAR}(obj, ...)
}
\arguments{
  \item{obj}{ object, class"FitAR" }
  \item{\dots}{optional arguments }
}
\value{
Plot produced as side-effect. No output.
}
\author{A.I. McLeod }

\seealso{ 
\code{\link{sdfplot}},
\code{\link{FitAR}}
 
}

\examples{
#Use AIC to select best subset model to fit to lynx data and
#plot spectral density function
data(SeriesA)
pvec<-SelectModel(SeriesA, ARModel="ARp", lag.max=10, Best=1)
ans<-FitAR(SeriesA, pvec)
sdfplot(ans)
#
#plot sdf and put your own title
z<-c(SeriesA)
pvec<-SelectModel(z, ARModel="ARp", lag.max=10, Best=1)
ans<-FitAR(z, pvec)
sdfplot(ans)
title(main="Example SDF")
}
\keyword{ ts }

