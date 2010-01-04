\name{binom.exact}
\alias{binom.exact}
\title{Exact tests with matching confidence intervals for single binomial parameter}
\description{
Calculates exact p-values and confidence intervals for a single binomial parmeter. 
This is different from \code{\link{binom.test}} only when alternative='two.sided', in which 
case \code{binom.exact} gives three choices for tests based on the 'tsmethod' option. 
The resulting p-values and confidence intervals will match.
}
\usage{
binom.exact(x, n, p = 0.5, 
   alternative = c("two.sided", "less", "greater"), 
   tsmethod = c("central", "minlike", "blaker"), 
   conf.level = 0.95, 
   control=binomControl())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{number of successes, or a vector of length 2 giving the
          numbers of successes and failures, respectively.
}
  \item{n}{
number of trials, ignored if x has length 2.
}
  \item{p}{
hypothesized probability of success.
}
  \item{alternative}{
indicates the alternative hypothesis and must be one of
          "two.sided", "greater" or "less".  You can
          specify just the initial letter.
}
  \item{tsmethod}{
indicates the method for a two-sided alternative hypothesis
 and must be one of
          "minlike", "central" or "blaker".  You can
          specify just the initial letter.
}
  \item{conf.level}{
confidence level for the returned confidence interval.
}
  \item{control}{list with settings to avoid problems with ties, etc, should not need to change this for normal use, 
   see \code{\link{binomControl}} 
}
}
\details{
Traditionally, hypothesis tests and confidence intervals are treated separately.
A more unified approach suggested by Hirji (2006) is to use the same p-value 
function to create confidence intervals. There is essentially only one way to 
calculate one-sided p-values and confidence intervals so these methods are the same 
in \code{\link{binom.test}} and \code{binom.exact}. However,  
there are three main ways that 
\code{binom.exact} allows for defining two-sided p-values.
\tabular{ll}{
minlike: \tab sum probabilities of all likelihoods equal or less than observed \cr
central: \tab double minimum one-sided p-value \cr
blaker: \tab combine smaller observed tail probability with opposite tail not greater than observed tail \cr
} 
The 'minlike' method is the p-value that has been used in \code{\link{binom.test}},
and 'blaker' is described in Blaker (2000) or Hirji (2006), where it is called the 
'combined tails' method. 
Once the p-value function is defined we can invert the test to create 'matching' 
confidence intervals defined as the smallest interval that contains all parameter values for which the two-sided
hypothesis test does not reject. There are some calculation issues for the 
'minlike' and 'blaker' methods which are the same as for exact tests 
for 2x2 tables (see Fay, 2009).
}
\value{
An object of class 'htest': a list with items
\item{p.value}{p-value}
\item{conf.int}{ confidence interval, see attributes 'conf.level' and perhaps 'conf.limit.prec' }
\item{statistic}{ number of successes}
\item{parameter}{number of trials}
\item{estimate}{observed proportion of success}
\item{null.value}{null hypothesis probability of success, 'p'} 
\item{alternative}{a character string describing alternative hypothesis}
\item{method}{a character string  describing method}
\item{data.name}{a character string giving the names of the data}
}
\references{
Blaker, H. (2000) Confidence curves and improved exact confidence intervals for discrete distributions.
Canadian Journal of Statistics 28: 783-798.

Fay, M. P. (2009). Confidence intervals that Match Fisher's exact and Blaker's exact tests.
  (to appear Biostatistics).

Hirjim K. F. (2006). Exact analysis of discrete data. Chapman and Hall/CRC. New York.

Stern, T (1954). Some remarks on confidence and fiducial limits.
Biometrika, 275-278.

Vos, P.W. and Hudson, S. (2008). Problems with binomial two-sided tests and the associated
confidence interals. Aust. N.Z. J. Stat. 50: 81-89.
}
\author{
M.P. Fay
}
\note{
The 'central' method gives the Clopper-Pearson intervals, and the 'minlike'
method gives confidence intervals proposed by Stern (1954) (see Blaker, 2000).
The 'blaker' method is guaranteed to be more powerful than the 
'central' method (see Blaker, 2000, Corollary 1), 
but both the 'blaker' method and 'minlike' method may have 
some undesireable properties. For example, there are cases where adding an 
additional  Bernoulli observation REGARDLESS OF THE RESPONSE will increase the p-value,
see Vos and Hudson (2008). The 'central' method does not have 
those undesireable properties.

The Blyth-Still-Casella intervals  given in StatXact (and not by \code{binom.exact}) are the shortest possible 
intervals, but those intervals are not nested. This means that the Blyth-Still-Casella 
intervals are not guaranteed 
to have the 95 percent interval contain the 90 percent interval. See Blaker (2000) Theorem 2.   
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{binom.test}}, for two-sample exact binomial tests see \code{\link[exact2x2]{exact2x2}}
}
\examples{
## Notice how binom.test p-value is given by tsmethod='minlike'
## but the confidence interval is given by tsmethod='central'
## in binom.exact p-values and confidence intervals match
binom.test(10,12,p=20000/37877)
binom.exact(10,12,p=20000/37877,tsmethod="minlike")
binom.exact(10,12,p=20000/37877,tsmethod="central")
binom.exact(10,12,p=20000/37877,tsmethod="blaker")
## two-sided methods are also available
## as in binom.test

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
