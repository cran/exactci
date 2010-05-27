\name{poisson.exact}
\alias{poisson.exact}
\title{Exact Poisson tests with Matching Confidence Intervals}
\description{
     Performs an exact test of a simple null hypothesis about the
     rate parameter in Poisson distribution, or for the
     ratio between two rate parameters. This is different from \code{\link{poisson.test}}
     in that 3 different types of exact two-sided tests (and the matching confidence intervals)
     are offered. The one-sided tests are the same as in \code{\link{poisson.test}}.
}
\usage{
poisson.exact(x, T = 1, r = 1,
    alternative = c("two.sided", "less", "greater"), tsmethod=c("central","minlike","blaker"),
    conf.level = 0.95, control=binomControl(), plot=FALSE)
}
\arguments{
  \item{x}{number of events. A vector of length one or two.}
  \item{T}{time base for event count. A vector of length one or two. }
  \item{r}{hypothesized rate or rate ratio}
  \item{alternative}{indicates the alternative hypothesis and must be
    one of \code{"two.sided"}, \code{"greater"} or \code{"less"}.
    You can specify just the initial letter.}
  \item{tsmethod}{character giving two-sided method, one of "central", "minlike" or "blaker", ignored if alternative not equal "two.sided"}
  \item{conf.level}{confidence level for the returned confidence
    interval.}
  \item{control}{list with settings to avoid problems with ties, etc, should not need to change this for normal use, 
   see \code{\link{binomControl}} }
  \item{plot}{logical, plot p-value function? For finer control on plot see \code{\link{exactpoissonPlot}}}
}
\details{
  Confidence intervals are computed similarly to those of
  \code{\link{binom.exact}} in the one-sample case, in that there are three two-sided options depending on the 
\code{tsmethod}. For the one-sample case the 
default intervals use \code{tsmethod="central"} giving the Garwood (1936) exact central confidence intervals.  
For the two-sample case we condition on the total counts and then use binomial methods,
see Lehmann and Romano (2005) for that motivation and \code{vignette("exactci")} for description 
of the three different two-sided methods for calculating p-values and confidence intervals.
}
\value{
  A list with class \code{"htest"} containing the following components:
  \item{statistic}{the number of events (in the first sample if there
    are two.)}
  \item{parameter}{the corresponding expected count}
  \item{p.value}{the p-value of the test.}
  \item{conf.int}{a confidence interval for the rate or rate ratio.}
  \item{estimate}{the estimated rate or rate ratio.}
  \item{null.value}{the rate or rate ratio under the null,
    \code{r}.}
  \item{alternative}{a character string describing the alternative
    hypothesis.}
  \item{method}{the character string \code{"Exact Poisson test"} or
    \code{"Comparison of Poisson rates"} as appropriate.}
  \item{data.name}{a character string giving the names of the data.}
}
\references{
Garwood, F (1936). Fiducial limits for the Poisson distribution. Biometrika, 437-442.

Lehmann, EL, and Romano, JP (2005). Testing Statistical Hypotheses, third edition. Springer:New York.
}
\note{
  The rate parameter in Poisson data is often given based on a
  \dQuote{time on test} or similar quantity (person-years, population
  size). This is the role of the \code{T} argument.
}
\seealso{
  \code{\link{poisson.test}},  \code{\link{exactpoissonPlot}}, 
}
\examples{
### Suppose you have observed rates of  2 out of 17877 in group A 
### and 10 out of 20000 in group B
### poisson.test gives non-matching confidence intervals
### i.e., p-value using  'minlike' criteria but confidence interval using 'central' criteria
poisson.test(c(2,10),c(17877,20000))
### poisson.exact gives matching CI to the p-values
### defaults to 'central' two-sided method
poisson.exact(c(2,10),c(17877,20000))
### other options
poisson.exact(c(2,10),c(17877,20000),tsmethod="minlike")
poisson.exact(c(2,10),c(17877,20000),tsmethod="blaker")
}
\keyword{htest}

