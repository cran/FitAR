\name{UnitRootTest}
\alias{UnitRootTest}
\title{ Unit Root Test }
\description{
Unit root test. Test H0: rho=1 vs. H1: rho<1 in the model
model with intercept z[t] = const + rho*z[t-1] + a[t].
}
\usage{
UnitRootTest(z, method = c("MLE", "ExactMLE", "LS", "All"), statistic = c("Z", "T"), NumBoot = 1000, PValueMethod = c("DH", "ET"))
}
\arguments{
  \item{z}{ time series }
  \item{method}{ estimation methods }
  \item{statistic}{ normalized rho or t-statistic }
  \item{NumBoot}{ number of bootstrap iterations }
  \item{PValueMethod}{ p-value can be estimated either as (k+1)/(N+1)
  as recommended by Davison and Hinkley (p. 148) or as k/N as in
  Efron and Tibsharini (p. 221, Algorithm 16.1 ).  }
}

\details{
Bootstrap unit root tests
}

\value{
one-sided P-value
}

\references{ 
Davison, A.C. and Hinkley, D.V. (1997). Bootstrap Methods and their Application. Cambridge.

Efron, B. and Tibshirani, R. (1993). An Introduction to the Bootstrap. Chapman/Hall.

Yu, H., Zhang, Y. and McLeod, A.I. (2009).  
Unit Root Bootstrap Tests with Exact Maximum Likelihood.
}

\author{ A.I. McLeod }

\seealso{ 
\code{\link{PP.test}}
}

\examples{
\dontrun{#takes about 10 seconds
z<-cumsum(rnorm(100))
UnitRootTest(z)
}
}

\keyword{ ts }
\keyword{ htest }
