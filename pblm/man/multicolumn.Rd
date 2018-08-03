\name{multicolumn}
\alias{multicolumn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
transforming bivariate data in a multi-column format
}
\description{
This function transforms a grouped two-column response data frame into a multi-column one, another data format accepted by \code{pblm}
}
\usage{
multicolumn(formula,data)
}
\arguments{
  \item{formula}{a two-side formula with counts in the left side. }
  \item{data}{a data frame with two categorical responses, covariates (if any)
       and a count variable.}
}
\value{
A data frame with as many responses as the number cells from the underlying 
response table and covariates (if any).
}
\author{
Marco Enea
}

\examples{
#NOT RUN 
data(ulcer)
multicolumn(freq~medication+pain+operation,data=ulcer)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{ regression }
\keyword{ models }
\keyword{ nonparametric }
\keyword{ smooth }
