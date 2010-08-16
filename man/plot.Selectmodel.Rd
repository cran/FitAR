\name{plot.Selectmodel}
\alias{plot.Selectmodel}
\title{ Subset AR Graph for "Selectmodel" Object }
\description{
  A graphical depiction is given of the output from \code{SelectModel}.
}
\usage{
\method{plot}{Selectmodel}(x, ...)
}
\arguments{
  \item{x}{ out from \code{SelectModel} }
  \item{...}{ optional arguments }
}
\details{
  The relative plausibility of Model A vs the best Model B, is defined as
  \eqn{R = e^{(AIC_B-AIC_A)/2}}.
  Values of R less than 1% indicate an implausible model.
  R is defined similarly if the BIC/UBIC criterion is used.
}
\value{
No value. Plot produced as side-effect.
}

\author{ A.I. McLeod }

\seealso{ \code{\link{SelectModel}} }
\examples{
#takes about 10 seconds
\dontrun{
out<-SelectModel(log(Willamette),lag.max=150,ARModel="AR",Best=5,Criterion="AIC")
plot(out)
}
}
\keyword{ ts }
