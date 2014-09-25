\name{exactpoissonPlot}
\alias{exactpoissonPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot p-value function for single or pair of poisson responses.
}
\description{
Plots p-values as a function of different point null hypothesis values for rate. 
For two-sided p-values, can plot three types of p-values: the minimum likelihood method (default for poisson.test), the central 
method (i.e., twice the one-sided exact p-values), and Blaker's exact.
}
\usage{
exactpoissonPlot(x,
    T=1,
    r=NULL,
    ndiv=1000,
    tsmethod="central",
    rRange=NULL,
    dolog=TRUE,
    dolines=FALSE,
    dopoints=TRUE,
    doci=TRUE, 
    alternative = c("two.sided", "less", "greater"), 
    relErr=1 + 10^(-7), 
    conf.level=.95,
    alphaline=TRUE,
    newplot=TRUE,
    midp=FALSE,...)
}
\arguments{
  \item{x}{number of events. A vector of length one or two}
  \item{T}{time base for event count. A vector of length one or two}
  \item{r}{null values of rate for plot, if NULL divides rRange into ndiv pieces}
  \item{ndiv}{number of pieces to divide up range of x-axis}
  \item{tsmethod}{two-sided method for p-value calculation, either "minlike","blaker" or "central"}
  \item{rRange}{range for plotting null hypothesis values of rate, if null then uses confidence interval to determine range}
  \item{dolog}{logical, plot horizontal axis in log scale?}
  \item{dolines}{logical, add lines to a plot?}
  \item{dopoints}{logical, add points to a plot?}
  \item{doci}{logical, add lines for confidence interval?}
  \item{alternative}{type of alternative for p-values}
  \item{relErr}{number close to 1, avoids problems with ties, see \code{\link{binomControl}} }
  \item{conf.level}{confidence level for use when doci=TRUE}
\item{alphaline}{logical, if doci=TRUE should line be drawn at significance level} 
 \item{newplot}{logical,start a new plot?}
\item{midp}{logical, use mid-p for p-values? Not available for tsmethod='blaker' or 'minlike'}
  \item{\dots}{values passed to plot, points, or lines statement}
}
%\details{}
\value{
Does graph or adds lines or points. Returns  (invisibly, see \code{\link{invisible}}) a list with elements r (null hypothesis values) and p.value (associated p-values).
}
%\references{%% ~put references to the literature/web site here ~}
%\author{%%  ~~who you are~~}
%\note{%%  ~~further notes~~}
\seealso{\code{\link{binom.exact}}}
\examples{
## single Poisson response
exactpoissonPlot(2,17877)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
