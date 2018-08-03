\name{summary.pblm}
\alias{summary.pblm}
\alias{print.pblm}
\alias{print.summary.pblm}
\alias{AIC.pblm}
\alias{logLik.pblm}
\alias{vcov.pblm}
\alias{coef.pblm}
\alias{coefficients.pblm}
\alias{edf.pblm}
\alias{residuals.pblm}
\alias{resid.pblm}
\alias{fitted.pblm}
\alias{predict.pblm}
\alias{deviance.pblm}
\alias{se.smooth.pblm}
\alias{chisq.test.pblm}
\alias{Rsq.pblm}
\title{Summarizing methods for bivariate additive logistic regression}
\description{
  Summarizing methods anf functions for objects of class \code{pblm}.
}
\usage{
\method{summary}{pblm}(object,...)
 
\method{print}{pblm}(x,digits = max(3, getOption("digits") - 3),...)

\method{print}{summary.pblm}(x,digits = max(3, getOption("digits") - 3),...) 
 
\method{AIC}{pblm}(object,...,k=2) 

\method{logLik}{pblm}(object, penalized=FALSE,...)

\method{vcov}{pblm}(object,...)

\method{coef}{pblm}(object, digits = max(3, getOption("digits") - 3), ...)

\method{coefficients}{pblm}(object, digits = max(3, getOption("digits") - 3), ...)

\method{residuals}{pblm}(object, type = c("working", "pearson"),...)

\method{resid}{pblm}(object, type = c("working", "pearson"),...)

\method{fitted}{pblm}(object,...)

\method{predict}{pblm}(object, newdata, type=c("link","response","terms","count"), 
                         se.fit=FALSE, digits= max(6, getOption("digits") - 3),...)

\method{deviance}{pblm}(object, penalized=FALSE,...)

edf.pblm(object, which.var=1, which.eq=1)

se.smooth.pblm(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ an object of class \code{pblm}. }
  \item{digits}{controls the number of digits printed in the output.}
  \item{x}{an object to be printed produced by  \code{pblm} or \code{summary.pblm}.}
  \item{k}{numeric, the penalty per parameter to be used; by default k = 2, that is 
  the classical AIC.}
  \item{penalized}{logical indicating whether the value of the penalized log-likehihood 
  is required}
  \item{which.var}{variable index indicating the position of the smoother 
  as this appears in the model formula.}
  \item{which.eq}{equation index for the smoothers.}
  \item{newdata}{if provided this represents a data.frame with new values of model covariates.}
  \item{type}{the type of desidered residuals or prediction.}
  \item{se.fit}{logical. Do you want prediction standard errors to be returned? 
               Currently allowed only for \code{type="link"} when \code{newdata} is not specified.}
  \item{\dots}{further arguments. }
}
%\details{
  
%}
%\value{
%A list the following components:
%  \item{results}{a p x 4 data frame with columns for the estimated coefficient, 
%                 its standard error, z-statistic and corresponding (two-sided) 
%                 p-value.}
%  \item{convergence}{a logical flag indicating whether the algorithm converged.}
%  \item{logLikp}{the penalized log-likelihood of the model fitted.}  
%  \item{AIC}{the AIC of the model fitted.}  
%  \item{BIC}{the BIC of the model fitted.}  
%  \item{names}{two-vector containing the names of the response variables.}
%  \item{Df}{degree of freedom of the model fitted.}  
%  \item{tr}{number of parameter effectively estimated.}
%  \item{call}{the mathed call.}  
%}
%\references{ ~put references to the literature/web site here ~ }
\author{Marco Enea}
\seealso{ \code{\link{pblm}}, \code{\link{pb}}}
\examples{

#NOT RUN 
# an artificial data set: 
set.seed(10)
da <- expand.grid("Y1"=1:3,"Y2"=1:3,"fat1"=0:4,"fat2"=0:1)
da$Freq <- sample(1:20,3*3*5*2,replace=TRUE)
da$x1 <- rnorm(90)
#the bivariate additive proportional-odds model
m7 <- pblm(fo1=cbind(Y1,Y2) ~ fat1 + pb(x1), data=da, weights=da$Freq)
summary(m7)
par(mfrow=c(1,3))
plot(m7,which.eq=1,xlab="x1",main="marginal equation 1")
plot(m7,which.eq=2,xlab="x1",main="marginal equation 2")
plot(m7,which.eq=3,xlab="x1",main="association")

}
\keyword{ multivariate }
\keyword{ regression }
\keyword{ models }
\keyword{ nonparametric }
\keyword{ smooth }
