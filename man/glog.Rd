\name{glog}
\alias{glog}
\title{ glog transformation }
\description{
The glog is a better behaved log transformation when some data 
values are zero or just near zero.
}
\usage{
glog(x, a = 1, InverseQ = FALSE)
}

\arguments{
  \item{x}{ numeric vector of data }
  \item{a}{ additive constant, often 1}
  \item{InverseQ}{ inverse glog }
}
\details{
Basic properties of the glog transformation are illustrated
in the Mathematica notebook glog.nb and its pdf version
glog.pdf which are available in the package directory doc.
}

\value{
transformated data
}

\references{  
W. Huber, A. von Heydebreck, H. Sultmann, A. Poustka, and M. Vingron. Variance stablization
applied to microarray data calibration and to quantification of differential expression.
Bioinformatics, 18: S96-S10 2002.
 }

\author{ A.I. McLeod }

\seealso{ 
\code{\link{bxcx}}
}

\examples{
#usual log transformation doesn't work
all(is.finite(log(sunspot.month)))
#either shifted log
all(is.finite(log(sunspot.month+1)))
#or glog works
all(is.finite(glog(sunspot.month)))
#but glog may be better, especially for values <1 but >=0
}

\keyword{ manip }
