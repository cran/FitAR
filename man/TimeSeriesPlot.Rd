\name{TimeSeriesPlot}
\alias{TimeSeriesPlot}
\title{ Multi-Panel or Single-Panel Time Series Plot with Aspect-Ratio Control }
\description{
Cleveland (1993) pointed out that the aspect-ratio is important in graphically
showing the rate-of-change or shape information. For many time series,
it is preferably to set this ratio to 0.25 than the default.  In general,
Cleveland (1993) shows that the best choice of aspect-ratio is often obtained
by if the average apparent absolute slope in the graph is about 45 deg.
But for many stationary time series, this would result in an aspect-ratio
which would be too small. 
As a comprise we have chosen a default of 0.25 but the user can select 
other choices.
}
\usage{
TimeSeriesPlot(z, SubLength = Inf, aspect = 0.25, type="l", ylab=NULL, main=NULL, ...)
}
\arguments{
  \item{z}{ ts object or vector, time series data }
  \item{SubLength}{ maximum number of data points per panel. Default SubLength=Inf and regular graphics. For
  trellis graphics, set SubLength to a finite value. }
  \item{aspect}{ optional setting for the aspect-ratio }
  \item{type}{ plot type, default type="l" join points with lines}
  \item{ylab}{optional label for vertical axis }
  \item{main}{ optional title}
  \item{...}{ optional arguments passed to \code{xyplot} }
}
\details{
If z has attribute "title" containing a character string, this is used
on the plot.  
Time series input using the function \code{\link{Readts}} always
have this attribute set.
}

\value{
If \code{SubLength} is finite, the lattice package is used and a graphic object of class trellis is produced.
Otherwise, the standard R graphics system is used and the plot is produced as a side-effect and there is no output.
}

\references{ W.S. Cleveland (1993), \emph{Visualizing Data}. }
\author{ A.I. McLeod }

\note{ Requires \code{lattice} library }


\seealso{ 
\code{\link{plot.ts}},
\code{\link{Readts}}
 }

\examples{
#from built-in datasets
 TimeSeriesPlot(AirPassengers)
 title(main="Monthly number of trans-Atlantic airline passengers")
#
#compare plots for lynx series
plot(lynx)
TimeSeriesPlot(lynx, type="o", pch=16, ylab="# pelts", main="Lynx Trappings")
#
#lattice style plot
data(Ninemile)
TimeSeriesPlot(Ninemile, SubLength=200)
}

\keyword{ ts }
