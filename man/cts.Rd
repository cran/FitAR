\name{cts}
\alias{cts}
\title{ Concantenate Time Series }
\description{
Concatenate two ts objects
}
\usage{
cts(x, y)
}
\arguments{
  \item{x}{ one time series, a ts object }
  \item{y}{ a vector or ts object which is concatenated on to x }
}
\details{
If y is a ts object, it is first converted to a vector.
Then the vector is concantenated on to the time series x.
}
\value{
A time series which starts at start(x) and has length equal
to length(x)+length(y).
}
\author{ A.I. McLeod }

\seealso{ 
\code{\link{ts}}, 
\code{\link{start}}, 
\code{\link{window}}
 }

\examples{
#Example 1
#Compare cts and c
#In the current version of R (2.6), they produce different results
z1<-window(lynx,end=1900)
z2<-window(lynx,start=1901)
z<-cts(z1,z2)
y<-c(z1,z2)

#See also Example 2 in predict.FitAR documentation
}

\keyword{ ts }
