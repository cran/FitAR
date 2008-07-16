\name{cts}
\alias{cts}
\title{ Concantenate Time Series }
\description{
Creating a ts object by concatenating y onto x, where x is a ts object and y is a vector.
If y is a ts object, it is simply converted to a vector and then concatenated on to x and
the tsp attribute of y is ignored.
}
\usage{
cts(x, y)
}
\arguments{
  \item{x}{ a time series, a ts object }
  \item{y}{ a vector which is concatenated on to x }
}
\details{
If y is a ts object, it is first converted to a vector.
Then the vector is concantenated on to the time series x.
An error is given if x is not a ts object.
}

\value{
A time series which starts at start(x) and has length equal
to length(x)+length(y).
}
\author{ A.I. McLeod }

\note{ 
The package zoo may also be used to concatenate time series, as in this example,

     x <- ts(1:3)
     y <- ts(4:5, start = 4)
     z <- ts(6:7, start = 7)
     library("zoo")
     as.ts(c(as.zooreg(x), y, z))
}

\section{Warning }{
Only two arguments are allowed, otherwise an error message will be given.

}

\seealso{ 
\code{\link{ts}}, 
\code{\link{start}}, 
\code{\link{window}}
\code{\link[zoo]{as.zooreg}} 
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

#Example 3. 
#Note tsp attribute of second argument is ignored but a warning is given if it is present
# and not aligned with first argument's attribute.
x <- ts(1:3)
z <- ts(6:7, start = 7)
cts(x,z) #warning given
y <- ts(4:5, start = 4)
cts(x,y) #no warning needed in this example.

}

\keyword{ ts }
