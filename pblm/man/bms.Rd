\name{bms}
\alias{bms}
\docType{data}
\title{
A British male sample on occupational status.}
\description{
These data were analyzed by Goodman (1979) and concern the cross-classification of a sample of fathers and their sons according to the occupational status.
}
\usage{data(bms)}
\format{
  A data frame with 49 observations and 3 variables. 
  \describe{
    \item{\code{fathers}}{fathers'occupational status. A factor with levels from 1 to 7.}
    \item{\code{sons}}{sons'occupational status. A factor with levels from 1 to 7.}
    \item{\code{freq}}{a vector of integers representing the number of people cross-classified according to the occupational status of \code{fathers} and \code{sons}}
  }
}

\source{
   Goodman, L. A. (1979). Simple models for the analysis of 
   cross-classications having ordered categories.  \emph{Journal of the 
   American Statistical Association},  \bold{74}(367), 537-552.\cr
}

\examples{
data(bms)
xtabs(freq ~ fathers + sons, data=bms)
}
\keyword{datasets}
