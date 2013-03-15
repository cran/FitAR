\name{toBinary}
\alias{toBinary}

\title{Binary representation of non-negative integer}
\description{
A non-negative integer is represented as a binary number.
The digits, 0 or 1, of this number are returned in a vector.
}
\usage{
toBinary(n, k = ceiling(logb(n+1,base=2)))
}
\arguments{
  \item{n}{
a non-negative integers
}
  \item{k}{
number of digits to be returned. 
}
}

\value{
A vector of length k. The first element is the least significant digit.
}
\author{
A.I. McLeod
}

\examples{
toBinary(63)
toBinary(64)
#sometimes we want to pad result with 'leading' 0's
toBinary(63, k=20)
toBinary(64, k=20)
}

\keyword{arith}
