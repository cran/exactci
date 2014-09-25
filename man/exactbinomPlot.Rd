\name{exactbinomPlot}
\alias{exactbinomPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot p-value function for one binomial response.
}
\description{
Plots p-values as a function of different point null hypothesis values for p. 
For two-sided p-values, can plot three types of p-values: the minimum likelihood method (default for binom.test), the central 
method (i.e., twice the one-sided exact p-values), and Blaker's exact.
}
\usage{
exactbinomPlot(x, n, p = NULL, ndiv = 1000, 
    tsmethod = "central", pRange = c(0, 1), 
    dolines = FALSE, dopoints = TRUE, doci=TRUE, 
    alternative=c("two.sided","less","greater"),
    relErr=1+10^(-7), conf.level=.95, alphaline=TRUE, 
    newplot = TRUE, midp=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{number of successes, or a vector of length 2 giving the numbers of successes and failures, respectively}
  \item{n}{number of trials, ignored if x has length 2}
  \item{p}{null values of p for plot, if NULL divides pRange into ndiv pieces}
  \item{ndiv}{number of pieces to divide up range of x-axis}
  \item{tsmethod}{two-sided method for p-value calculation, either "minlike","blaker" or "central"}
  \item{pRange}{range for plotting null hypothesis values of p}
  \item{dolines}{logical, add lines to a plot?}
  \item{dopoints}{logical, add points to a plot?}
  \item{doci}{logical, add lines for confidence interval?}
  \item{alternative}{type of alternative for p-values}
  \item{relErr}{number close to 1, avoids problems with ties, see \code{\link{binomControl}} }
  \item{conf.level}{confidence level for use when doci=TRUE}
\item{alphaline}{logical, if doci=TRUE should a line be drawn at significance level}
  \item{newplot}{logical,start a new plot?}
  \item{midp}{logical, use mid-pvalues? See \code{link{binom.exact}} }
  \item{\dots}{values passed to plot, points, or lines statement}
}
%\details{}
%\value{}
%\references{%% ~put references to the literature/web site here ~}
%\author{%%  ~~who you are~~}
%\note{%%  ~~further notes~~}
\seealso{\code{\link{binom.exact}}}
\examples{
## plot central two-sided p-values (double one-sided p-values)
## for 5 positive responses out of 17 tries
exactbinomPlot(5,17)
## add blakers exact p-values
## pch option acts on points, lty acts on ci lines
exactbinomPlot(5,17,tsmethod="blaker",col="blue",pch=".",lty=2,newplot=FALSE)

## can plot one-sided p-values, tsmethod is ignored
exactbinomPlot(5,17,alternative="less")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
