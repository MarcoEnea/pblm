\name{pblm.control}
\alias{pblm.control}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Auxiliary for controlling the algorithm in a \code{pblm} model
%%  ~~function to do ... ~~
}
\description{
This is an auxiliary function for controlling the algorithm in a \code{pblm} model. 
}
\usage{
pblm.control(maxit = 30, maxit2 = 200, acc = 1e-07, acc2 = 1e-06, 
            zero.adj = 1e-06, l = NULL, restore.l = FALSE,    
            min.step.l = 1e-04, auto.select = FALSE, gaic.m = 2, 
            rss.tol = 1e-06, max.backfitting = 10, pgtol.df = 0.01, 
            factr.df = 1e+07, lmm.df = 5, parscale.df = 1, 
            max.gaic.iter = 500, pgtol.gaic = 1e-05, grad.tol = 1e-07, 
            factr.gaic = 1e+07, lmm.gaic = 5, parscale = 1, 
            conv.crit = c("dev", "pdev"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{maxit}{maximum number of Fisher-scoring iterations.}
  \item{maxit2}{maximum number of Newton-Raphson iterations for the inversion 
                \eqn{\eta}->\eqn{\pi}. }               
  \item{acc}{tolerance to be used for the estimation.}
  \item{acc2}{tolerance to be used for the inversion \eqn{\eta}->\eqn{\pi}. }
  \item{zero.adj}{adjustment factor for zeros in the probability vector \eqn{\pi}{pi}.}
  \item{l}{numerical, ranged in (0,1], representing the initial value of step lenght. By default \code{l=1}.}  
  \item{restore.l}{logical, should the step length be restored to its initial value after each iteration? This is an experimental option and may be changed in the future.}  
  \item{min.step.l}{numerical, minimum value fixed for the step length.}
  \item{auto.select}{logical, should the smoothing parameters be estimated by GAIC minimization? If \code{TRUE} The optimization will be performed numerically by using \code{optim}.
  }
\item{gaic.m}{the "penalty" per parameter of the generalized AIC. By default it is 2, corresponding to the classical AIC.
}
\item{rss.tol}{tolerance for the residual sum of squares used in the backfitting algorithm.
}
\item{max.backfitting}{maximum number of backfitting iterations.} 
\item{pgtol.df}{tolerance to be used  in order to get an amount of smoothing corresponding to the fixed degrees of freedom for the additive part. See argument \code{pgtol} from \code{\link{optim}}.
} 
\item{factr.df}{numerical. For degrees-of-freedom optimization in the additive part. See argument \code{factr} from \code{\link{optim}}.
} 

\item{lmm.df}{integer. For degrees-of-freedom optimization in the additive part. See argument \code{lmm} from \code{\link{optim}}.
}  
 
\item{parscale.df}{A vector of scaling parameters for vector lambda when optimizing lambda for fixed degrees of freedom. See argument \code{parscale} from \code{\link{optim}}.
}
\item{max.gaic.iter}{integer. Maximum number of iterations for automatic model optimization. See argument \code{maxit} from \code{\link{optim}}.
}
\item{pgtol.gaic}{numerical. Tolerance to be used for automatic selection of smoothing parameters. See argument \code{pgtol} from \code{\link{optim}}.
}   
\item{grad.tol}{numerical. Tolerance to be used when inverting the gradient matrix.
}
\item{factr.gaic}{numerical. For automatic selection of smoothing parameters. See argument \code{factr} from \code{\link{optim}}.
}
\item{lmm.gaic}{integer. For automatic selection of smoothing parameters. See argument \code{lmm} from \code{\link{optim}}.
}
\item{parscale}{A vector of scaling parameters for vector lambda for automatic model optimization. See argument \code{parscale} from \code{\link{optim}}.
}
\item{conv.crit}{Convergence criterion for model estimation. The default is "dev", corresponding to log-likelihood maximization. Alternatively, "pdev" is concerned with maximum penalized log-likelihood.
}                   

}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A list with the same arguments of the function, unless unlikely specified by the user.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Marco Enea %%  ~~who you are~~
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{pblm}}
}
%\examples{
%
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{ regression }
\keyword{ models }
\keyword{ nonparametric }
\keyword{ smooth }