\name{JarqueBeraTest}
\alias{JarqueBeraTest}
\title{ Jarque-Bera Normality Test }
\description{
 A powerful omnibus test for normality.
}
\usage{
JarqueBeraTest(z)
}
\arguments{
  \item{z}{ vector of data }
}
\details{
  This test is derived as a Lagrange multiplier test for normal
  distribution in the family of Pearson distributions (Jarque and Bera, 1987 ).
}  
\value{
  \item{LM }{value of the LM statistic}
  \item{pvalue }{p-value}
}
\references{ 
Jarque, C.M. and Bera, A.K. (1987).  
A Test for Normality of Observations and Regression Residuals.   
International Statistical Review 55, 163-172}

\author{ A.I. McLeod }


\examples{
#some normal data
z<-rnorm(100)
JarqueBeraTest(z)
#some skewed data
z<-rexp(100)
JarqueBeraTest(z)
#some thick tailed data
z<-rt(100,5)
JarqueBeraTest(z)
}
\keyword{ htest }
