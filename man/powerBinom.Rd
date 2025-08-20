\name{powerBinom}
\alias{powerBinom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Exact Power and Sample Size methods for Bernoulli responses 
}
\description{
Calculates sample sizes or power for a study of n independent Bernoulli responses (i.e., a binomial response with parameter n). There are three types of calculations. For type='standard' calculate the usual sample size or power under a given alternative. The option type='cilength' bases calculations on the expected length of the confidence intervals. For type='obs1ormore' calculate the sample size or power to observe 1 or more positive responses. All calculations use exact methods. 
}
\usage{
powerBinom(n = NULL, p0 = 0.5, p1 = NULL, 
    sig.level = 0.05, power = NULL, 
    alternative = c("two.sided", "one.sided"), strict = FALSE,
    type = c("standard","cilength","obs1ormore"), 
    cilength = NULL, conf.level = 0.95, 
    direction = c("greater","less"), 
    control = binomControl(), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
sample size
}
  \item{p0}{
probability of success under null
}
  \item{p1}{
probability of success under alternative
}
  \item{sig.level}{
significance level of test
}
  \item{power}{
target power of test (type='standard') or probability of observing at least 1 success (type='obs1ormore')
}
  \item{alternative}{
one- or two-sided test
}
  \item{strict}{
use strict interpretation in two-sided case (count rejections in opposite direction)
}
  \item{type}{
either 'standard', 'cilength', or 'obs1ormore' (see details)
}
  \item{cilength}{
average length of confidence interval (used when type='cilength')
}
  \item{conf.level}{
confidence level (used when type='cilength')
}
  \item{direction}{
direction of alternative, either 'greater' or 'less' (used when type='standard' and p1 is NULL) 
}
  \item{control}{
a list with arguments that control algorithms, see 
\code{\link{binomControl}}
}
  \item{\dots}{
further arguments passed to \code{\link{binom.exact} (see details)}
}
}
\details{
Type='standard' calculates the power to reject the null hypothesis with the parameter p0 on the boundary between the null and alternative. In other words, the null could be either of the three, H0: p=p0 (for alternative='two.sided') or H0: p<=p0 or H0:p>=p0 (for alternative='one.sided'). For one-sided alternatives, p1 determines the alternative (and hence determines the null). For example, p0=.5 and p1=.7 will test H0: p<=0.5 vs H1: p> 0.5. 

Type='cilength' calculates the expected length of the confidence intervals when p=p1. If p1=NULL, then it uses p1=0.5 because this will give the largest CI lengths. For sample size calculations, the n value is found using the \code{\link[ssanv]{uniroot.integer}} function searching the range from control$minn to control$maxn, finding the smallest n that has expected CI length at least as large as cilength.

Type='obs1ormore' calculates sample sizes related to the probability (i.e., the power) to observe at least one success. 

Here are some details of the calculation method for type='standard' and type='cilength'. We use control$pRange (which should be from a vector of length two giving an interval very close to (0,1) but excluding the ends) to save computation time. We do not need to calculate the powers or expected CI lengths for all the possible values of X from 0:n. The algorithm only uses X from the likely values, using control$pRange to determine the quantiles that dominate the power or expectation calculation. For example, if n=1000 and pRange[1]=10^-10 and p1=0.5, then there is no need to calculate the CI lengths for x=0:399, because pbinom(399,1000,.5)< 10^-10.
    
For the \dots argument, you cannot pass the 'alternative' argument to \code{\link{binom.exact}}, since it has a different form  in \code{powerBinom} since p1 determines which one-sided alternative will be used. 

}
\value{
a object of class \code{power.htest}. List with input arguments plus calculated values, together with 'note' and 'method' character vectors describing the method. 
}
%\references{}
\author{
Michael P. Fay
}
\note{Often you will get the same power with midp=FALSE (default) and midp=TRUE. This is because the rejection region may be the same for both cases. For example, for type='standard' with n=20 using the default two.sided 0.05 exact central test, regardless of whether midp=FALSE or TRUE, we reject p0=0.5 for x in (0:5 and 15:20).}

%% \section{Warning }{....} 

%\seealso{%% ~~objects to See Also as \code{\link{help}}, ~~~}
\examples{
# find power to reject H0: p = p0
# when p=p1 (only counting rejections in the correct direction)     
# defaults to exact central test of binom.exact      
powerBinom(n=40,p0=.4,p1=.6)
# The following calculates the sample size
# to have the average length of confidence intervals 
# be at least 0.4 (regardless of the true proportion)
powerBinom(type="cilength",cilength=.4)
# The following answers the question:
# if the true proportion is .001, how many 
# observations do you need to sample to 
# have the probability of seeing at least 1 success be
# at least .9?
powerBinom(power=.9,p1=.001,type="obs1ormore")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}