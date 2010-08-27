\name{BoxCox}
\alias{BoxCox}
\title{ Generic Box-Cox Analysis Function }
\description{
The function is implemented as a generic function with methods for 
classes "FitAR", "Arima", "ts" and "numeric".


For \eqn{\lambda \ne 0}{lambda!=0}, the Box-Cox transformation is of x is
\eqn{(x^\lambda - 1)/\lambda}{(x^lambda-1)/lambda}.
If the minimum data value is <= 0, a small positive constant, equal to
the negative of the minimum plus 0.25, is added to all the
data values.

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
\dontrun{#takes a few seconds
BoxCox(lynx)
out<-FitAR(lynx, c(1,2,4,10,11), ARModel="ARp", MLEQ=FALSE)
BoxCox(out)
out<-FitAR(lynx, c(1,2,4,5,7,8,10,11,12))
BoxCox(out)
}
}

\keyword{ ts }

