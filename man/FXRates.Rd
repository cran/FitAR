\name{FXRates}
\alias{FXRates}
\docType{data}
\title{ Foreign exchange rates}
\description{
Daily foreign exchange rates were obtained for: YenUS, DmUS, USGB, CanUS. 
From 1983-12-13 to 2008-11-12, 6330 values in each series.
}

\usage{data(FXRates)}

\format{
  A data frame with 6330 observations on the following 5 variables.
  \describe{
    \item{\code{Date}}{a character fector, dates}
    \item{\code{YenUS}}{Yen/US exchange rate}
    \item{\code{DmUS}}{Deutsche Mark/US Dollar, exchange rate}
    \item{\code{USGB}}{US/Great Britain exchange rate}
    \item{\code{CanUS}}{Canada/US exchange rate}
  }
}

\details{
The dates run from "1983-12-13" to "2008-11-24"  and were included in the downloaded file.
There were 48 missing values out of a total of 4*6330=25320 values. 
Missing values were replaced with the previous value. 
}

\source{
\url{http://www.econstats.com/fx/fx__d1.htm}
}



\examples{
data(FXRates)
head(FXRates)
}

\keyword{datasets}
