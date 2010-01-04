\name{exactci-package}
\alias{exactci-package}
\alias{exactci}
\docType{package}
\title{
Exact binomial and Poisson tests with Matching Confidence Intervals
}
\description{
Calculates exact binomial and Poisson tests giving matching confidence 
intervals. There are 3 different methods for defining the two-sided 
p-values.  
}
\details{
\tabular{ll}{
Package: \tab exactci\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2010-01-04\cr
License: \tab GPL2\cr
LazyLoad: \tab yes\cr
}
Although \code{\link{binom.test}} and \code{\link{poisson.test}} give exact 
tests and confidence intervals, for two-sided tests the confidence intervals (CI) are not formed by inverting 
the tests. Thus, there may be test-CI inconsistencies whereby the test rejects 
but the confidence interval contains the null parameter. The \code{exactci} package eliminates many of these 
test-CI inconsistencies for  
two-sided tests, by outputing the matching confidence interval with each test.
The package uses one of three different 
methods for defining the two-sided p-value. The main functions of the package 
are \code{\link{binom.exact}} and \code{\link{poisson.exact}} which follow the same format 
as  \code{binom.test} and \code{poisson.test} except have the option 'tsmethod'
to define the two-sided method for calculating the p-values, and give matching 
confidence intervals (i.e., ones that come from the inversion of the p-values). 
}
\author{
M.P. Fay

Maintainer: Mike Fay <mfay@niaid.nih.gov>
}
\references{
Blaker, H. (2000) Confidence curves and improved exact confidence intervals for discrete distributions.
Canadian Journal of Statistics 28: 783-798.

Fay, M. P. (2009). Confidence intervals that Match Fisher's exact and Blaker's exact tests.
  (to appear Biostatistics).

Hirjim K. F. (2006). Exact analysis of discrete data. Chapman and Hall/CRC. New York.
}
\keyword{ package }
\keyword{htest}
\seealso{
For comparisons of two binomial groups see \code{\link[exact2x2:exact2x2-package]{exact2x2}}
}
\examples{
## Note binom.test calculates p-values using priciple of minimum likelihood
## while it calculates the central confidence intervals. That is why the 
## inferences do not match in this example.
binom.test(10,12,p=20000/37877)
binom.exact(10,12,p=20000/37877,tsmethod="minlike")
binom.exact(10,12,p=20000/37877,tsmethod="central")
## We also allow the method studied in Blaker (2000)
binom.exact(10,12,p=20000/37877,tsmethod="blaker")
}