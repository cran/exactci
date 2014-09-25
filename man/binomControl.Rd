\name{binomControl}
\alias{binomControl}
\title{Tuning parameters for binom.exact function}
\description{
This function produces a list of tuning parameters used in the calculations done by \code{\link{binom.exact}}
and \code{\link{poisson.exact}}. These will not need to be changed by most ordinary users.
}
\usage{
binomControl(relErr=1+1e-07,tol=.00001,
    pRange=c(1e-10,1-1e-10),
    minn=1, maxn=10^4, 
    PRINT.STEPS=FALSE)
}
\arguments{
  \item{relErr}{
value very close to 1, used in calculation of two-sided p-values
}
  \item{tol}{
value very close to 0, used in calculation of two-sided confidence intervals
}
  \item{pRange}{
range close to [0,1], but excluding the endpoints, used in calculation of two-sided confidence intervals
}
\item{minn}{
minimum n used by powerBinom. The search for the lowest n that gives power is crude and goes from minn to maxn by ones.
}
\item{maxn}{
maximum n used by powerBinom
}
\item{PRINT.STEPS}{
logical, print steps of uniroot.integer used to find the n when type='cilength' in powerBinom?
}
}
\details{
See the code for \code{fisher.test}, where the term \code{relErr} is hard-coded into the function. The purpose 
is to avoid problems with ties. It serves the same purpose in this package and probably need not be changed.
The value \code{tol} indicates the tolerance for the precision of the confidence limits.
The value \code{pRange} is input into \code{uniroot} to give bounds when searching for confidence limits. For poisson 
limits \code{pRange} is transformed using the \code{qgamma} function (see code in \code{exactpoissonCI}).
}
\value{
  A list with containing the following components:
  \item{relErr}{a number larger than 1}
  \item{tol}{a number greater than 0}
  \item{pRange}{a vector with 2 elements between 0 and 1, exclusive}
}
%\note{}
%\seealso{}
%\examples{}
\keyword{misc}
