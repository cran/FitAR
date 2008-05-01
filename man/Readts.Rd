\name{Readts}
\alias{Readts}
\title{ Input a Time Series }
\description{
This function inputs time series stored in ASCII in a format that
the first line in the file is a title, next few lines, beginning with a \#, are comments,
and the remaining lines contain the data.  Here is an example: \cr

\code{
Changes In Global Temperature, Annual, 1880-1985  \cr
\#Surface air "temperature change" for the globe, 1880-1985.  \cr
\#Degrees Celsius. "Temperature change" actually means temperature  \cr
\#Surface Air Temperature", `Journal of Geophysical Research`, Vol. 92, \cr
-.40    -.37    -.43    ....	 \cr
...............			\cr
 .27     .42     .02     .30     .09     .05  \cr
}

}
\usage{
Readts(file = "", freq = 1, start = 1, VerboseQ=TRUE)
}
\arguments{
  \item{file}{ location for input file }
  \item{freq}{ tsp parameter, =1, annual, =12 monthly etc }
  \item{start}{ tsp parameter }
  \item{VerboseQ}{ normally prompt for arguments but set
  VerboseQ=FALSE to automate }
}
\value{
ts object with attribute 'title'
}
\author{ A.I. McLeod }

\seealso{ 
\code{\link{scan}}, 
\code{\link{ts}}, 

 }
\examples{
#You will need to change save the data given above in a file
#and change the directory as appropriate
#z<-Readts(file="d:/datasets/mhsets/annual/globtp.1", start=1880, VerboseQ=FALSE)
}
\keyword{ ts }
