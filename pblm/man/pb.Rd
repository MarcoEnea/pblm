\name{pb}
\alias{pb}
\alias{pb.control}
\alias{pbs}
\title{
Specify a Penalised B-Spline Fit in a pblm Formula
}
\description{
Both \code{pb} and\code{pbs}are adaptations of function \code{\link[gamlss]{pb}} and \code{\link[gamlss]{ps}} from the \code{\link[gamlss]{gamlss}} package, respectively, to specify penalized B-spline.}
\usage{
pb(x, df = NULL, lambda = NULL, control = pb.control(...), ...)
pb.control(inter = 20, degree = 3, order = 2, quantiles = FALSE, ...) 
pbs(x, df = 3, lambda = NULL, ps.intervals = 20, degree = 3, order = 3)
}
\arguments{
  \item{x}{the univariate predictor. }
  \item{df}{the desidered equivalent number of degrees of freedom (trace of the smoother matrix minus two for the constant and linear fit.}  
  \item{lambda}{the smoothing parameter.}
  \item{control}{setting the control parameters}
  \item{ps.intervals}{the number of break points in the x-axis. }
  \item{inter}{the number of break points (knots) in the x-axis.}
  \item{degree}{the degree of the piecewise polynomials.}
  \item{order}{the required difference in the vector of coefficients.}
  \item{quantiles}{if \code{TRUE} the quantile values of x are used to determine the 
         knots.}
 \item{...}{for extra arguments.} 
}

\details{
Basically, \code{pb} is a reduced-functionality version of the original one specified in \code{gamlss} with no performance iteration methods (i.e. by \code{method} specification) implemented. The only method implemented minimizes the GAIC by internal \code{\link[stats]{optim}} calls.
}

\author{Marco Enea, based on the original versions of the corresponding functions contained in the \code{gamlss} package by Mikis Stasinopoulos and Bob Rigby.}

\examples{
#NOT RUN 

# an artificial data set: 
set.seed(1234)
da <- expand.grid("Y1"=1:3,"Y2"=1:3,"fat1"=0:9,"fat2"=0:1)
da$x1 <- seq(-5,5,l=180)
da$x2 <- rnorm(180)

da$Freq <- sample(5:30,180,replace=TRUE)

m1 <- pblm(fo1=cbind(Y1,Y2) ~ pbs(x1) + fat2, 
           fo2=~pb(x1) + x2, 
           fo12=~pb(x1) + x2, data=da, weights=Freq)

par(mfrow=c(2,3))
plot(m1)
}

\keyword{ multivariate }
\keyword{ regression }
\keyword{ models }
\keyword{ nonparametric }
\keyword{ smooth }
