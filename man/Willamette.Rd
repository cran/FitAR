\name{Willamette}
\alias{Willamette}
\docType{data}
\title{ Willamette Riverflow Time Series}

\description{
Monthly flows of the Willamette River, Salem, Oregon, Oct. 1950 - Aug. 1983.
Percival and Walden (1993, Ch. 10.15) fit high-order AR models to estimate the
spectral density function.
}

\usage{data(Willamette)}

\format{
  The format is:
 Time-Series [1:395] from 1951 to 1984:  8.95  9.49 10.19 10.96 11.08 ...
 - attr(*, "title")= chr "Willamette river, Monthly, Salem, Oregon, Oct. 1950 - Aug. 1983"
}

\source{
\url{http://faculty.washington.edu/dbp/sapabook.html}  
}

\references{
Percival, D.B. and Walden, A.T. (1993),
\emph{Spectral Analysis for Physical Applications}.
Cambridge University Press.   
}

\examples{
#Percival and Walden (1993) fit an AR(27).
#Compare spectral densities with subset AR's.
data(Willamette)
pmax<-27
sdfplot(FitAR(log(Willamette), pmax))
p<-SelectModel(log(Willamette), ARModel="ARz", lag.max=pmax, Best=1)
sdfplot(FitAR(log(Willamette), p))
p<-SelectModel(log(Willamette), ARModel="ARp", lag.max=pmax, Best=1)
sdfplot(FitAR(log(Willamette), p), ARModel="ARp", MLEQ=FALSE)
}
\keyword{datasets}
