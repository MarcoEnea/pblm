\name{pblm}
\alias{pblm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bivariate Additive Regression for Categorical Responses}

\description{This function allows to fit bivariate additive marginal logistic regression models for nominal, 
ordinal or mixed nominal/ordinal responses. 
}
\usage{
pblm(fo1=NULL, fo2=NULL, fo12=NULL, RC.fo=NULL, data, weights=NULL, 
     type="gg", verbose=FALSE, contrasts=NULL, start=NULL, x=FALSE, 
     center=FALSE, scale=FALSE, plackett=NULL, ncat1=NULL, ncat2=NULL, 
     fit=TRUE, proportional=pblm.prop(...), penalty=pblm.penalty(...), 
     control=pblm.control(...), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fo1}{a two hand-side formula for the logit(s) of the first response.  
             Depending on the data structure, the left-hand side can be a 
             two-column or a multi-column response vector. In the latter case, argument \code{ncat1} must be specified. See Examples.}
  \item{fo2}{a one hand-side formula for the logit(s) of the second response. 
             If omitted, it will be assumed to be equal to \code{fo1}. }
  \item{fo12}{a one hand-side formula for the log-odds ratio(s) of the 
              association between the two responses. If omitted, it will 
              be assumed to be equal to \code{fo1}.}                        
  \item{RC.fo}{a Row/Column type formula specifying the association structure. See Details.}   
  \item{data}{a data frame.}
  \item{weights}{An optional vector containing the observed frequencies.}
  \item{type}{a two-length string specifying the type of logits to use in the model fit. See Details.}
   \item{verbose}{logical. Should information about convergence be printed during model estimation? 
}
  \item{contrasts}{
   the Row/Column contrasts to be used in \code{RC.fo}. See argument \code{contrasts.arg} in \code{\link{model.matrix.default}}.  
}
  \item{x}{logical. Should the model matrix used in the fitting process be returned as component of the fitted model?
}
  \item{center}{logical. Should the covariates be centered with respect their mean before the fit?    
}
  \item{scale}{logical. Should the covariates be scaled with respect their standard deviation before the fit?
}
  \item{start}{an optional vector of starting values for the coefficients of the non-additive part.
}
 \item{plackett}{
     logical. Should a Plackett-based formula be used for the inversion \eqn{\eta}->\eqn{\pi}. Actually, this is allowed (and it is the default) for binary or ordered responses only. The more general method uses the Newton-Raphson inversion algorithm described in Glonek and McCullagh (1995).
} 
  \item{ncat1}{an integer indicating the number of levels of the first response. Mandatory the data is in multicolumn format. See example below.
}
  \item{ncat2}{an integer indicating the number of levels of the second response. 
}
  \item{fit}{logical. If TRUE (default), the model will be estimated, otherwise only some objects created prior to estimation will be returned. 
}
  \item{control}{this sets the control parameters for the Fisher-scoring and the inner Newton-Raphson and backfitting algorithms. The default setting is specified by the 
  \code{\link{pblm.control}} function. 
  }
\item{proportional}{this sets a list of logical vectors specifying which explanatory variables depend on the categories of the responses. The default setting is specified by the \code{\link{pblm.prop}} function 
}
\item{penalty}{this sets the penalty terms and smoothing parameters for a non-parametric "vertical" smoothing across the categories of the responses. The default setting is specified by the \code{\link{pblm.penalty}} function. 
}
\item{\dots}{ further arguments. }
}
\details{
It is possible to fit partial proportional odds models and specify several association
structures like the Goodman's (1979) model (interactions are allowed though these are of
linear type), the Dale's (1986) model and their additive version (Bustami \emph{et al.}, 2001), Enea \emph{et al.}(2014).
Furthermore, the association structure can also be smoothed by using penalty terms of
polynomial type as those considered in Enea and Attanasio (2015). That allows to enlarge
the range of possible parametrizations of the association structure as an alternative
to the Dale's Row-Column parametrization. Further details on the penalty terms specified through \code{\link{pblm.penalty}} can also found in Enea and Lovison (2014). \cr 

The algorithm is based on the bivariate version of the model by Glonek and McCullagh (1995),
that is by using  the \emph{multivariate logit transform} \deqn{\bm{C}'\log(\bm{M}\bm{\pi})=
\bm{\eta}=\bm{X}\bm{\beta}.}{C'log(M\pi)=\eta=X\beta.}

Once \code{fo1} has been specified, if \code{fo2} and \code{fo12} are left unspecified,
these are assumed to be equal to \code{fo1}. By default, the function fits a proportional
odds model for ordered responses, in which only the marginal and the association intercepts
are assumed to be category-dependent (Glonek and McCullagh, 1995).\cr


Model formulae using a Row-Column type parametrization, like in the Goodman or the
Dale models, need to be specified by using \code{RC.fo}. The right-hand side of such
formula only recognizes  \code{Row} and \code{Col} to specify row and/or column effects.
Covariates must be specified separately by using \code{fo12}. See Examples.\cr
The logits implemented are \emph{local}; \emph{global}; \emph{continuation};
\emph{reverse continuation}; (Colombi and Forcina, 2001) and \emph{basic}.
By using argument \code{type}, several log-odds ratios can be specified, among the
permutations of the local-local type (\code{type="ll"}), local-global (\code{"lg"}),...,
basic-basic (\code{type="bb"}). Furthermore, if the responses are binary, setting
\code{type="ss"} will correspond to classical logit
\eqn{\pi = \log P(Y=1)/P(Y=0)}{\pi = log P(Y=1)/P(Y=0)}
for both responses, while other specifications will produce
\eqn{\log P(Y=0)/P(Y=1)}{log P(Y=0)/P(Y=1)}. \cr
         
The vector of the starting values must be set in the following order: itercepts of
the first logit; covariates of the first logit; intercepts of the second logit; 
covariates of the second logit; association intercepts; association covariates.\cr

For what concerning the additive part, p-slines can be fitted by using \code{\link{pb}}
and/or \code{\link{pbs}}, which are adaptations (reduced versions) of \code{\link[gamlss]{pb}}
and \code{\link[gamlss]{ps}} from the \code{gamlss} package, respectively.  
}

\value{
 A list of components from the model fit. Some of these could be redundant and removed in a
 future version of the package. When \code{fit=TRUE} the returned components are:\cr
  \item{coef}{a named column vectors of coefficients.}
  \item{n}{the total number of observations.}
  \item{m}{the number of observed configurations of the responses given covariates (if any).
           It corresponds to the number of rows of the dataset.}
  \item{p}{fitted probability matrix given the observations.}
  \item{Y}{the weighted matrix of the responses in multinomial format.}
  \item{x}{if requested, the model matrix.}
  \item{xx1, xx2, xx12}{vectors, for internal use.}
  \item{ynames}{vector with the names of the responses.}
  \item{tol}{the accuracy reached at the convergence.}
  \item{llp}{the penalized log-likelihood value at the convergence.}
  \item{ll}{the log-likelihood value at the convergence.}
  \item{devp}{the penalized deviance value at the convergence.}
  \item{dev}{the deviance value at the convergence.}
  \item{IM}{the estimated Information Matrix at the convergence.}
  \item{IMp}{the estimated penalized Information Matrix at the convergence.
             Its inverse is used to calculate the standard error of estimates.}
  \item{convergence}{logical indicating whether convergence criteria were satisfied.}
  \item{iter}{the number of iterations performed in the model fitting.}
  \item{maxB}{the number of smoothers present in the equation with the maximum number
              of smoothers.This also represents the number of outer iterations of the
              backfitting algorithm.}
  \item{ncat1, ncat2}{the number of levels for the two responses.}
  \item{ncat}{this is simply \code{ncat1*ncat2}.}
  \item{weights}{the prior weights, that is the weights initially supplied, a vector
                 of 1s if none were.}
  \item{P}{a \eqn{p \times p}{p x p} penalty matrix, where \emph{p} is the number of
           parameters estimated (including all the intercepts), concerning the penalty
           term specified (if any), but excluding the additive part (if any). If no
           penalty terms are specified, this will be a matrix of zeros.}
  \item{gaic.m}{the \emph{penalty} per parameter of the generalized AIC initially
                provided. It is 2 by default.}
  \item{lam1, lam2, \cr lam3, lam4}{ vectors of smoothing parameters initially supplied
        with the specification of a penalty term. If not, these are \code{NULL}. }
  \item{opt}{the object returned by \code{\link{optim}} if automatic selection of smoothing
             parameters has been set. \code{NULL} otherwise.}
   \item{etasmooth}{ the matrix of predicted values for the additive part, \code{NULL} otherwise.}
  \item{eta}{ the \code{n}\eqn{\times}{x}(\code{ncat}-1) matrix of predicted values
             on the observed data.}
   \item{fsmooth}{ a list of objects initially created by the smoothers (if any),
                  \code{NULL} otherwise.}
   \item{one.smooth}{for internal use.}
   \item{df.fix}{ the degrees of freedom of the parametric part of the model.
                  Note that if a penalty term is used by specifying \code{penalty},
                  the corresponding degrees of freedoms will be counted in the
                  parametric part of the model.}
   \item{df.fix.vect}{for internal use.}
   \item{df.smooth}{ the degrees of freedom of the additive part of the model (if any),
                    otherwise zero.}
   \item{w.res}{ a \code{n}\eqn{\times}{x}(\code{ncat}-1) matrix of working residuals.}
   \item{W2}{a \code{m}-length list of matrices of working weights.}
   \item{z}{a (\code{m}*\code{ncat}-1)-length working vector of responses.}
   \item{any.smoother}{logical indicating whether smoothers have been used in the model fit.}
   \item{Bmat}{list of bases of smoothers used in the model fit (if any), \code{NULL} otherwise.}
   \item{wh.eq}{list of logical vectors indicating which terms in the linear predictor are smoothers.}
   \item{wh.eq2}{list of logical vectors indicating which terms in the linear predictor are smoothers. For internal use.}
   \item{PBWB}{list of numerical matrices if smoothers are used, \code{NULL} otherwise. For internal use.} 
   \item{BWB}{list of numerical matrices if smoothers are used, \code{NULL} otherwise. For internal use.}
   \item{etalist}{list of numerical vectors if smoothers are used. For internal use. }
   \item{spec.name}{list of numerical vectors if smoothers are used. For internal use. }
   \item{beta.smooth}{list of estimated coefficients for the smoothers (if any), \code{NULL} otherwise.}
   \item{n.smoothers}{the number of smoothers used.}
   \item{PPmat}{list of penalty matrices for the smoothers multiplied by the smoothing parameters. \code{NULL} if smoothers are not used.}
   \item{pnl.type}{the type of penalty used.}
   \item{xnames}{list of vectors with the names of the covariates.}
   \item{prop.smooth}{for internal use.}
   \item{GAIC}{the value of the generalized AIC from the fitted model for the specified value of \code{gaic.m} in \code{\link{pblm.control}}.}
\item{ta}{the underlying observed contingency table of the responses, marginally to the covariates.}
 \item{set0}{the parameter setting from \code{\link{pblm.control}} as specified by the user.}
\item{fo.list}{list of formulas used.}
\item{center}{logical indicating whether argument if the covariates were centered}
\item{scale}{logical indicating whether argument if the covariates were scaled}
\item{type}{the type of logit used.}
\item{plackett}{logical indicating whether the plackett inversion formula was used.}
\item{RC.fo}{the Row-Column formula as specified by the user.}
\item{contrasts}{the Row-Column formula contrasts as specified by the user.}
\item{acc2}{tolerance to be used for the estimation as specified by the user.}
\item{maxit}{maximum number of Fisher-scoring iterations as specified by the user.}
} 
  

\seealso{ \code{\link{summary.pblm}}}

\references{Bustami, R., Lesaffre, E., Molenberghs, G., Loos, R., Danckaerts, M. and 
            Vlietinck, R. 2001. Modelling bivariate ordinal responses smoothly with  
            examples from ophtalmology and genetics. \emph{Statistics in Medicine}, 
            \bold{20}, 1825-1842.\cr 
                      
            Colombi, R. and Forcina, A. 2001. Marginal regression models for the
            analysis of positive association of ordinal response variables. 
             \emph{Biometrika},  \bold{88}(4), 1007-1019.\cr
            
            Dale, J. R. (1986) Global Cross-Ratio Models for Bivariate, Discrete, 
            Ordered Responses.  \emph{Biometrics},  \bold{42}(4), 909-917.\cr

	          Enea, M. and Attanasio, M. 2015. A model for bivariate data with application 
            to the analysis of university students success. \emph{Journal of Applied 
            Statistics}, http://dx.doi.org/10.1080/02664763.2014.998407 \cr
            
            Enea, M. and Lovison, G. 2014. A Penalized penalized approach to the 
            bivariate logistic regression model for the association between ordinal 
            responses. arXiv:1407.1751 [math.ST]\cr 
            
	          Enea, M.,  Stasinopoulos, M., Rigby, R., and Plaia, A. 2014. The pblm package: 
            semiparametric regression for bivariate categorical responses in R. In 
            Thomas Kneib, Fabian Sobotka, Jan Fahrenholz, Henriette Irmer (Eds.)
            \emph{Proceeding of the 29th  International Workshop of Statistical Modelling},
             Volume 2, 47-50.\cr

            Glonek, G. F. V. and McCullagh, P. (1995) Multivariate logistic models.
             \emph{Journal of the Royal Statistical Society, Series B},  
             \bold{57}, 533-546.\cr 
                                  
            Goodman, L. A. (1979). Simple models for the analysis of 
            cross-classications having ordered categories.  \emph{Journal of the 
            American Statistical Association},  \bold{74}(367), 537-552.\cr
            
            Lapp, K., Molenberghs, G., and Lesaffre, E. (1998) Models for the 
            association between ordinal variables. \emph{Computational Statistics 
            and Data Analysis}, \bold{27}, 387-411.\cr 
            
}

\author{Marco Enea \email{marco.enea@unipa.it} with contribution by Mikis Stasinopoulos and Bob Rigby}
\note{
Please note that specifying a formula with interaction terms in \code{RC.fo}
corresponds to a model with association structure of the type
\eqn{\alpha + \beta_{r} + \gamma_{c} + \delta_{rc}}{\alpha + \beta[r] + \gamma[c] + \delta[rc]}.
In the current version of the package, non linear interaction terms of the type
\eqn{\alpha + \beta_r + \gamma_c + \delta_{1r}\delta_{2c}}{\alpha + \beta[r] + \gamma[c] + \delta[1r]\delta[2c]},
as considered for example in Lapp et al. (1998), are not implemented here.\cr
Furthermore, unlikely from the Dale's paramaterization, \code{pblm} does not put a minus
sign to covariates.\cr  
}
\section{Warning}{
Please be sure that the results you get are really what you are "expecting" when using
uncommon specifications of argument \code{type}, mainly those involving an
\code{s} logit type and its combinations with other logits. A part from the binary
case, many of those have been not checked yet in the current version of the package.\cr

The estimation of category-dependent p-splines, as outlined in Enea et al. (2014), is not 
allowed (work in progress).  

}
\examples{
## Example 1
#The Dale's model
data(ulcer)
m1 <- pblm(fo1=cbind(pain,medication)~1, fo12=~I(operation=="vh"), RC.fo=~Col,
           data=ulcer, weights=ulcer$freq,  contrasts=list(Col="contr.SAS"))

#compare with Dale (1986), Table 3           
summary(m1)
deviance(m1)

#the same data in another format
dat <- multicolumn(freq~medication+pain+operation,data=ulcer)
dat

fo <- as.formula(paste(attributes(dat)$"resp","~1",sep=""))

m1bis <- pblm(fo1=fo, fo12=~I(operation=="vh"), RC.fo=~Col, verbose=TRUE,
              data=dat, ncat1=3, contrasts=list(Col="contr.SAS"))

deviance(m1bis)

## Example 2
#NOT RUN 
# an artificial example: 
set.seed(10)
da <- expand.grid("Y1"=1:3,"Y2"=1:3,"fat1"=0:4,"fat2"=0:1)
da$Freq <- sample(0:20,3*3*5*2,replace=TRUE)
da$x1 <- rnorm(90)

#the bivariate additive proportional-odds model
m2 <- pblm(fo1=cbind(Y1,Y2) ~ fat1 + pb(x1), data=da, weights=da$Freq)
par(mfrow=c(1,3))
plot(m2)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{ regression }
\keyword{ models }
\keyword{ nonparametric }
\keyword{ smooth }
