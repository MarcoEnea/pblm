\name{ulcer}
\alias{ulcer}
\docType{data}
\title{
The ulcer data}
\description{
Data analyzed by Dale (1986)
}
\usage{data(ulcer)}
\format{
  A data frame with 48 observations and 4 variables. 
  \describe{
    \item{\code{medication}}{medication requirements. A factor with levels \code{never}; \code{seldom}; \code{occasionally}; and \code{regularly}.}
    \item{\code{pain}}{patients' post operative pain level. A factor with levels \code{none}, \code{slight} and \code{moderate}}. 
    \item{\code{operation}}{a factor representing the type of operation with levels: vagatomy drainage procedure (\code{vp}); vagatomy and distal antrectomy (\code{va}); vagatomy and hemigastrectomy (\code{vh});  and resection alone (\code{ra})}.
    \item{\code{freq}}{a numeric vector representing the number of patients classified for the corresponding levels of \code{pain}, \code{medication} and \code{operation}}
  }
}

\source{
            Dale, J. R. (1986) Global Cross-Ratio Models for Bivariate, Discrete, 
            Ordered Responses. \emph{Biometrics}, \bold{42}(4), 909-917.\cr
}

\examples{
data(ulcer)
}
\keyword{datasets}
