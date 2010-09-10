\name{SelectModel}
\alias{SelectModel}
\title{Select Best AR, ARz or ARp Model  }
\description{
The AIC/BIC/UBIC/EBIC/BICq criterion is used to select the best fitting AR or subset AR model.
When Best > 1, the result may be plotted using \code{plot}.
}
\usage{
SelectModel(z, lag.max = 15, ARModel = c("AR", "ARz", "ARp"), Criterion = "default", Best = 3, Candidates = 5, G = 0.5, Q = 0.25)
}
\arguments{
  \item{z}{ time series data }
  \item{lag.max}{ maximum order of autoregression }
  \item{ARModel}{ "AR" for full AR(p) or "ARp"/"ARz" corresponding to subset models}
  \item{Criterion}{ default is "BIC" for order selection and "BICq" for subset selection. 
  Options: "AIC", "BIC", "UBIC", "EBIC" and "BICq". }
  \item{Best}{ final number of models to be selected }
  \item{Candidates}{ number of models initially selected using the approximate criterion}
  \item{G}{ hyperparameter for EBIC }
  \item{Q}{ byperparameter for BICq }

}
\details{
  McLeod and Zhang (2006) outline an approximate AIC/BIC selection algorithm.
  This algorithm is a refinement of that method.  The refinement consists of automatically
  look for the best k candidates, where k = \code{Candidates}.  Then the exact likelihood
  is evaluated for all k candidates.  Out of these k candidates, the best q = \code{Best}
  are then selected.  This two-step procedure is needed because if k is too low,
  the approximate AIC/BIC rankings may not agree with the exact rankings. 
  This strategy is used for model selection for AR, ARz and ARp models.
  A plot method is available to graph the output.
  The UBIC developed by Chen and Chen (2007) is an important extension of the 
  BIC criterion for subset selection.  In the non-subset case UBIC is equivalent to BIC.
}
\value{
When \eqn{\code{Best} = 1}, a vector is returned indicated the lag or lags included
in the model.  The null model is indicated by returning 0 for the lag.
An object with class "Selectmodel" is returned when \eqn{\code{Best} >1}.
If ARModel = "AR", a matrix is return whose first column shows p and 
second AIC or BIC.
Otherwise for subset selection, the result is a list with q components, where q=BEST.
When Criterion = "UBIC", the components in this list are:
  \item{p}{lags present, a 0 indicates the null model}
  \item{UBIC}{exact UBIC}
similarly for the AIC/BIC case.

The components are arranged in order of the criterion used.

When ARModel = "ARp" or "ARz", an attribute "model" indicating "ARp" or "ARz"
is included.

}
\references{ 
McLeod, A.I. and Zhang, Y. (2006).  
Partial Autocorrelation Parameterization for Subset Autoregression.  
Journal of Time Series Analysis, 27, 599-612.

Chen, J. and Chen, Z. (2007).
Extended Bayesian Information Criteria for Model Selection with Large Model Space.
Preprint.

}
\author{ A.I. McLeod and Y. Zhang}

\note{For white noise, the best model is the null model, containing
no lags. This is indicating by setting the model order, \eqn{p=0}. }

\section{Warning }{
Setting \code{Candidates} too low can result in anomalous results.
For example if \code{Candidates} = 1, we find that the top ranking model
may depend on how large \code{Best} is set.
This phenomenon is due to the fact that among the best AIC/BIC models
there is sometimes very little difference in their AIC/BIC scores.
Since the initial ranking of the Candidates is done using the approximate
likelihood, the final ranking using the exact likelihood may change.
} 

\seealso{ 
\code{\link{plot.Selectmodel}},
\code{\link{PacfPlot}}, 
\code{\link{PacfPlot}}, 
\code{\link{FitAR}} 
}
\examples{
#Example 1: Find an ARp subset model for lynx data using BIC
z<-log(lynx)
out<-SelectModel(z, ARModel="ARp", Criterion="BIC", Best=5)
plot(out)
#
#Example 2: Find an ARz subset model for lynx data using BIC
out<-SelectModel(z, ARModel="ARz", Criterion="BIC", Best=5)
plot(out)
#
#Example 3: Select an AR(p) model
out<-SelectModel(z, ARModel="AR", Criterion="BIC", Best=5)
out
plot(out)
out<-SelectModel(z, ARModel="AR", Criterion="BIC", Best=1)
#
#Example 4: Fit subset models to lynx series
\dontrun{ #not run to save time in building package
z<-log(lynx)
#requires library leaps. Should be automatically when FitAR package is loaded.
#first fit ARp
pvec <- SelectModel(z, lag.max=11, ARModel="ARp", Criterion="AIC", Best=1)
ans1 <- FitAR(z, pvec, ARModel="ARp", MLEQ=FALSE)
#now fit ARz
pvec <- SelectModel(z, lag.max=11, ARModel="ARz", Criterion="AIC", Best=1)
ans2<-FitAR(z, pvec, ARModel="ARz")
#compare
summary(ans1)
summary(ans2)
#Use UBIC
pvec <- SelectModel(z, ARModel="ARp",lag.max=11,Best=1)
ans3<-FitAR(z, pvec, ARModel="ARp")
pvec <- SelectModel(z, ARModel="ARz",lag.max=11,Best=1)
ans4<-FitAR(z, pvec, ARModel="ARz")
#compare
summary(ans3)
summary(ans4)
}
}
\keyword{ ts }
