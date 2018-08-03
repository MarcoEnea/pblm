\name{pblm.penalty}
\alias{pblm.penalty}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Auxiliary for specifying penalty terms in a \code{pblm} model 
}
\description{
This is an auxiliary function for specifying penalty terms in a \code{pblm} model. 
}
\usage{
pblm.penalty(pnl.type=c("none","ARC1","ARC2","ridge","lasso","lassoV", "equal"),
             lam1=NULL, lam2=NULL, lam3=NULL, lam4=NULL, 
             s1=NULL, s2=NULL, s3=NULL, s4=NULL, min.lam.fix=0.1, 
             constraints=FALSE, lamv1=1e10, lamv2=1e10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
                
  \item{pnl.type}{The type of penalty term to be used. By default \code{"none"} 
                  of them is used. See Details.}
  \item{lam1,lam2}{vectors of smoothing parameters for the marginals. By default 
                   they are zero. It is common to all the penalty terms 
                   implemented.} 
  \item{lam3}{vector of smoothing parameters for the association. It is common 
              to all the penalty terms implemented. If \code{"ARC2"} is selected, 
              it will act on differences of by-row adjacent parameters.}
  \item{lam4}{vector smoothing parameters for the association. Specific 
                   for the "ARC2" penalty, for which the differences of by-column 
                   adjacent parameters are penalized}
  \item{s1,s2,s3,s4}{the orders of the difference operator. Specific for 
                     the "ARC2" penalty.}    
  \item{min.lam.fix}{the minimum value for any penalty parameters in order to consider any 
                     lamdas smaller than this value as fixed. See details.  }
  \item{constraints}{Should inequality constraints be applied to the marginal probabilities? 
                     This is done through a penalty term and intended for ordered resposes.}
  \item{lamv1, lamv2}{penalty parameters to be applied to both the marginal probabilities in
                      order to mimicking inequality constraints.}                                   

}
\details{
Some penalty terms implemented in \code{pblm} are described in Enea and Lovison (2014) and Enea and Attanasio (2015). \cr

Just one penalty per model is allowed.\cr

Penalty "ARC1" acts on first order differences of category-adjacent parameters. By 
increasing the smoothing parameters, the resulting marginal and/or association parameters,
will tend to be equal among the categories. When the underlying contingency table 
cross-classifying the responses contains zero cells, This penalty may be useful to stabilize 
the estimates, for example to get a more "regular" association structure. \cr

Penalty "ARC2" generalizes in a certain sense 
"ARC1", since it acts on high order differences of Ajacent Row and/or Column parameters, 
but it is maily used for ordered responses. For high smoothing values it constraints the 
marginal parameters to lie onto a polynomial curve, and/or the association structure to lie 
onto a polynomial surface. The degrees of the marginal polynomials are determined by \code{s1}-1,
for the first marginal, and \code{s2}-1 for the second. The degree of the association polynomial 
surface is determined by \code{s3}+\code{s4}-2, in which \code{s3}-1 is the by-row polynomial 
degree and \code{s4}-1 the by-column one.\cr

Penalty "ridge" constraints the regression parameters towards zero (horizontal penalty), so providing equal estimates for high penalty values.\cr 

The current implementation of the "lasso" and "lassoV" penalty terms is to be considerer provisional and needs to be better checked. \cr

Penalty "lasso" is acts similarly to "ridge" and it is based on absolute values of the regression parameters. \cr 

Penalty "lassoV" vertically penalizes the absolute value of differences of adjacent row and column parameters. It is similar to ARC1. By increasing the smoothing parameter, the resulting marginal and/or association regression parameters, will tend to be equal among the response categories. This\cr



Penalty "equal" constraints the marginal equations to be equal, so providing equal estimates. 
This could be useful, for example, in eyes or twins studies. The tuning parameter to be specified for this penalty is \code{lam1}.\cr  

Furthermore, if global logits are specified, inequality constraints on marginal predictors are mimicked by using penalty terms. Argument \code{lamv1} is the penalty parameter for the marginal predictor of the first response, \code{lamv2} is
that for the second one.\cr       

Argument \code{min.lam.fix} is useful when an automatic selection of penalty parameters is desidered for certain lambdas only. The remaining lambdas can be either 0 or less than \code{min.lam.fix}, and excluded from the automatic selection. Fixing some lambdas to assume values in \eqn{[}{\[}0,\code{min.lam.fix}\eqn{)}{\)} may be useful, for example, for parameter space regularization.  
}
\value{
A list with the same arguments of the function, unless unlikely specified by the user.
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Enea, M. and Attanasio, M. 2015. A model for bivariate data with application to the analysis of university students success. Journal of Applied Statistics, 
http://dx.doi.org/10.1080/02664763.2014.998407 \cr

Enea, M. and Lovison, G. 2014. A Penalized penalized approach to the bivariate logistic regression model for the association between ordinal responses. arXiv:1407.1751 [math.ST]\cr 
}
\author{Marco Enea \email{marco.enea@unipa.it}}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{pblm}}
}
\examples{

#Example 1
# A British male sample on occupational status. 
data(bms)

# A third degree polynomial surface with equally-spaced integer scores
m1 <- pblm(fo1=cbind(fathers,sons)~1, data=bms, weights=bms$freq,
             penalty=pblm.penalty(pnl.type="ARC2",lam3=c(1e7), lam4=c(1e7), 
                                  s3=c(4), s4=c(4)))

require(lattice)
g <- expand.grid("sons"=1:6,"fathers"=1:6)
g$logGOR <- m1$coef[13:48]

wireframe(logGOR ~ sons*fathers, data = g, zlim=c(min(g$logGOR-1),max(g$logGOR+1)), 
          scales = list(arrows = FALSE), screen = list(z = -130, x = -70),
          col.regions="magenta")



#Example 2

# an artificial data frame with two binary responses and two factors
set.seed(12)
da <- expand.grid("Y1"=0:1,"Y2"=0:1,"fat1"=0:1,"fat2"=0:1)
da$Freq <- sample(0:20,2*2*2*2,replace=TRUE)

# A quasi-independence model obtained by strongly penalizing the association intercept 
# through a ridge-type penalty term
m3 <- pblm(fo1=cbind(Y1,Y2) ~ fat1 + fat2, 
           fo12=~ 1, 
           data=da, weights=da$Freq, type="ss",
           proportional=pblm.prop(prop12=c(TRUE)),
           penalty=pblm.penalty(pnl.type="ridge",lam3=1e12))
summary(m3)

# notice that the last coefficient is not exactly zero
coef(m3)

m3.1 <- glm(Y1 ~ fat1 + fat2, data=da, weights=Freq, family=binomial)
m3.2 <- glm(Y2 ~ fat1 + fat2, data=da, weights=Freq, family=binomial)

all.equal(logLik(m3), logLik(m3.1)[1]+logLik(m3.2)[1])



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{ regression }
\keyword{ models }
\keyword{ nonparametric }
\keyword{ smooth }