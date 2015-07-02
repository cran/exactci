\name{exactci-internal}
\alias{exactci-internal}
\alias{exactpoissonPval}
\alias{exactpoissonPvals}
\alias{exactbinomCI}
\alias{exactpoissonCI}
\alias{exactbinomPvals}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Internal functions, not to be called directly
}
\description{These functions give p-values for exact binomial or poisson test. Not to be called directly, use \code{\link{binom.exact}} or 
\code{\link{poisson.exact}}}
\usage{
exactbinomPvals(x, n, p, relErr = 1 + 10^(-7), tsmethod = "minlike",midp=FALSE)
exactpoissonPval(x,T=1,r=1,relErr=1+1e-07,tsmethod="central")
exactpoissonPvals(x,T=1,r=1,relErr=1+1e-07,tsmethod="central")
exactbinomCI(x, n, tsmethod="minlike",conf.level=.95,tol=.00001,
    pRange=c(1e-10,1-1e-10))
exactpoissonCI(x, tsmethod="minlike",conf.level=.95,tol=.00001,
    pRange=c(1e-10,1-1e-10))
}
\arguments{
  \item{x}{number of successes or counts, vectors not allowed}
  \item{n}{number at risk, vectors not allowed}
  \item{p}{binomial parameter for null hypothesis, may be vector}
  \item{relErr}{used in calculation to avoid ties, slightly bigger than 1}
  \item{tsmethod}{two-sided method, one of "minlike", "blaker" ("central" only allowed for  \code{exactpoissonPval} or \code{exactpoissonPvals}) }
  \item{T}{number at risk or person-years at risk, vectors not allowed}
 \item{r}{rate parameter for null hypothesis, null is E(x*T)=r, vectors only allowed for \code{exactpoissonPvals}}
 \item{conf.level}{number between 0 and 1 for level of confidence interval}
 \item{pRange}{range to search for confidence intervals, between 0 and 1 (even for poisson where it is transformed to a 0 to Inf-like range)}
 \item{tol}{tolerance for precision of confidence interval, very small number}
\item{midp}{logical, use mid-p for p-values?}
}
\details{
The function \code{exactbinomPvals} tests point null hypotheses for a single binomial observation.
The function \code{exactpoissonPvals} tests point null hypotheses for a single Poisson observation. To get p-values for the two-sample 
Poisson test save results from \code{exactpoissonPlot}. 
The functions \code{exactbinomCI} and \code{exactpoissonCI} calculate the "minlike" and "blaker" confidence intervals.

}
\value{ Returns either a confidence interval with attributes giving precision, or a pvalue (exactpoissonPval), 
or a list with pvals and parameters (r,T for poisson and p for binomial).}
%\references{}
%\author{}
%\note{}
\seealso{ \code{\link{poisson.exact}} and \code{\link{binom.exact}}}
\examples{
exactbinomPvals(3,10,c(.3,.4,.5),tsmethod="minlike")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}