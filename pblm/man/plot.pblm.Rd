\name{plot.pblm}
\alias{plot.pblm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plotting smoothers for a \code{pblm} object
}
\description{
Plotting smoothers for a \code{pblm} object
}
\usage{
\method{plot}{pblm}(x, which.eq = 1:3, which.var = 1:x$maxNpred, add.bands = TRUE, 
                    col.bands = "lightblue", col.line = "black", type = "l", 
                    dashed.bands = FALSE, pause = FALSE, ylim, xlim, ylab, xlab, 
                    main, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Object of class \code{pblm}. }
  \item{which.var}{variable index for the smoothers indicating the rank of the smoother as this appear in the model formula.}
  \item{which.eq}{equation index for the smoothers.}
  \item{add.bands}{logical, should the confidence bands for the smoother be added to the graph?}
  \item{col.bands}{sets the color for the smoother confidence bands.}
  \item{col.line}{sets the color for the smoother line.}
  \item{type}{graphical parameter.}
  \item{dashed.bands}{logical. If TRUE (as well as add.bands) dashed (instead of colour filled) confidence bands will be drawn.}
  \item{pause}{logical. If TRUE the user will be asked to press a buttom to plot each graph.}
  \item{ylim}{graphical parameter.}
  \item{xlim}{graphical parameter.}
  \item{ylab}{graphical parameter.}
  \item{xlab}{graphical parameter.}
  \item{main}{graphical parameter.}
  \item{\dots}{ further arguments. } 
}
%%\details{
%%bels for different  
%%}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Marco Enea
}
\seealso{
\code{\link{pb}}, \code{\link{pbs}}
}
\examples{
#NOT RUN 
# an artificial data set: 
set.seed(10)
da <- expand.grid("Y1"=1:3,"Y2"=1:3,"fat1"=0:4,"fat2"=0:1)
da$Freq <- sample(0:20,3*3*5*2,replace=TRUE)
da$x1 <- rnorm(90)
#the bivariate additive proportional-odds model
m7 <- pblm(fo1=cbind(Y1,Y2) ~ factor(fat1) + pb(x1), data=da, weights=da$Freq)
par(mfrow=c(1,3))
plot(m7,which.eq=1,main="marginal equation 1")
plot(m7,which.eq=2,main="marginal equation 2")
plot(m7,which.eq=3,main="association")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{ regression }
\keyword{ models }
\keyword{ nonparametric }
\keyword{ smooth }
