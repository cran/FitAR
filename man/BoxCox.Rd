\name{BoxCox}
\alias{BoxCox}
\title{ Generic Box-Cox Analysis Function }
\description{
The function is implemented as a generic function with methods for classes "FitAR", "Arima", "ts" and "numeric".
}

\usage{
BoxCox(object, ...)
}

\arguments{
  \item{object}{ model object }
  \item{\dots}{ optional arguments}
}

\value{
No value returned. Graphical output is produced as side-effect.
The plot shows relative likelihood function as well as the MLE and a confidence interval.
}

\author{ A.I. McLeod and Y. Zhang}

\note{ The MASS package has a similar function boxcox but this is implemented
only for regression and analysis of variance. }

\seealso{ 
\code{\link{BoxCox.Arima}}, 
\code{\link{BoxCox.FitAR}}, 
\code{\link{BoxCox.ts}}, 
\code{\link{BoxCox.numeric}} 
 }

\examples{
BoxCox(lynx)
out<-FitAR(lynx, c(1,2,4,10,11), ARModel="ARp", MLEQ=FALSE)
BoxCox(out)
out<-FitAR(lynx, c(1,2,4,5,7,8,10,11,12))
BoxCox(out)
}

\keyword{ ts }

