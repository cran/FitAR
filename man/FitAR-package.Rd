\name{FitAR-package}
\alias{FitAR-package}
\docType{package}
\title{
Fits AR and subset AR models and provides complete model building 
capabilities.
FitAR
}
\description{

}
\details{
\tabular{ll}{
Package: \tab FitAR\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2008-01-08\cr
License: \tab GPL (>= 2)\cr
}
For model estimation the main function is FitAR for which generic methods
print, summary, coef and plot are implemented.
For model identification, there is a new PacfPlot for subset ARz idenfication.
Subset models may also be selected using AIC, BIC and UBIC criteria
with the function SelectModel.
Complete model diagnostic checking is available. 
The function plot provides many graphical diagnostic plots.

}
\author{
A. I. McLeod and Ying Zhang

Maintainer: aimcleod@uwo.ca
}

\references{ 
McLeod, A.I. and Zhang, Y. (2006).  
Partial autocorrelation parameterization for subset autoregression.  
Journal of Time Series Analysis, 27, 599-612.

McLeod, A.I. and Zhang, Y. (2008a). Faster ARMA Maximum Likelihood Estimation, 
Computational Statistics and Data Analysis
52-4, 2166-2176.  
DOI link: \url{http://dx.doi.org/10.1016/j.csda.2007.07.020}. 

McLeod, A.I. and Zhang, Y. (2008b, Submitted).
Improved Subset Autoregression: With R Package.
Journal of Statistical Software.
}

\keyword{ ts }

\keyword{ package }

\examples{

}
