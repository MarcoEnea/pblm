\name{pblm.prop}
\alias{pblm.prop}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Auxiliary for specyfing category-dependent covariates in a \code{pblm} model }

\description{This is an auxiliary function which allows to specify partially proportional odds for one (or both) the marginals and with the association parameters which can depend (or not) on the categories of the responses. It simply returns a list with its arguments.  \cr 
}
\usage{
pblm.prop(prop1=NULL, prop2=NULL, prop12=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{prop1}{a \code{TRUE}/\code{FALSE} logical vector specifying which explanatory variables are category-dependent, including the intercepts. Each element of this vector must exactly match the order of the covariates appearing in \code{fo1}, leaving out any additive term. If a factor is present, the corresponding elements in this vector will have as many elements as the number of levels of the factor minus 1. }
  \item{prop2}{a logical vector like \code{prop1}, but now it refers to the covariates present in \code{fo2}.}
  \item{prop12}{a logical vector like \code{prop1}, but now it refers to the covariates present in \code{fo12}}   
}
\details{

The default specification will result in a model with category-dependent intercepts for both the marginal and the association, while all the covariates are assumed independent of the categories.\cr  
Note that, for ordered responses, setting category-independent intercepts for the marginals is not a good idea.
}

\value{
A list with the same arguments of the function, unless unlikely specified by the user.
}

\author{Marco Enea \email{marco.enea@unipa.it}}

\examples{
# an artificial data frame with two five-category responses and two factors
set.seed(10)
da <- expand.grid("Y1"=1:5,"Y2"=1:5,"fat1"=letters[1:3],"fat2"=letters[1:3])
da$Freq <- sample(1:20,5*5*3*3,replace=TRUE)

#A partial proportional-odds model with uniform association
m2 <- pblm(fo1=cbind(Y1,Y2) ~ fat1 + fat2, 
           fo2=~fat1,
           fo12=~1, 
           data=da, weights=da$Freq, 
           proportional=pblm.prop(prop1=c(FALSE,TRUE,TRUE,FALSE,FALSE),
           prop2=c(FALSE,TRUE,TRUE),
           prop12=c(TRUE)))
summary(m2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{ regression }
\keyword{ models }
\keyword{ nonparametric }
\keyword{ smooth }
