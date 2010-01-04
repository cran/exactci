binom.exact<-function (x, n, p = 0.5, alternative = c("two.sided", "less", 
    "greater"), tsmethod=c("central","minlike","blaker"), conf.level = 0.95, 
    control=binomControl()) 
{
    relErr<-control$relErr
    tol<-control$tol
    pRange<-control$pRange
    ## copied setup code from binom.test, Version 2.10.1
    DNAME <- deparse(substitute(x))
    xr <- round(x)
    if (any(is.na(x) | (x < 0)) || max(abs(x - xr)) > 1e-07) 
        stop("'x' must be nonnegative and integer")
    x <- xr
    if (length(x) == 2) {
        n <- sum(x)
        x <- x[1L]
    }
    else if (length(x) == 1) {
        nr <- round(n)
        if ((length(n) > 1) || is.na(n) || (n < 1) || abs(n - 
            nr) > 1e-07 || (x > nr)) 
            stop("'n' must be a positive integer >= 'x'")
        DNAME <- paste(DNAME, "and", deparse(substitute(n)))
        n <- nr
    }
    else stop("incorrect length of 'x'")
    if (!missing(p) && (length(p) > 1 || is.na(p) || p < 0 || 
        p > 1)) 
        stop("'p' must be a single number between 0 and 1")
    alternative <- match.arg(alternative)
    if (!((length(conf.level) == 1) && is.finite(conf.level) && 
        (conf.level > 0) && (conf.level < 1))) 
        stop("'conf.level' must be a single number between 0 and 1")
    ## end of copied setup code

    tsmethod<-match.arg(tsmethod)
    if (tsmethod!="central" & tsmethod!="minlike" & tsmethod!="blaker") stop("tsmethod must be one of 'central', 'minlike', or 'blaker' ")

    PVAL <- switch(alternative, less = pbinom(x, n, p), greater = pbinom(x-1,n,p,lower.tail=FALSE), 
        two.sided = exactbinomPvals(x,n,p,relErr=relErr,method=tsmethod)$pvals)

    p.L <- function(x, alpha) {
        if (x == 0) 
            0
        else qbeta(alpha, x, n - x + 1)
    }
    p.U <- function(x, alpha) {
        if (x == n) 
            1
        else qbeta(1 - alpha, x + 1, n - x)
    }
    if (alternative=="less"){
        CINT<-c(0, p.U(x, 1 - conf.level))
    } else if (alternative=="greater"){
        CINT<-c(p.L(x, 1 - conf.level), 1)
    } else {
        if (tsmethod=="central"){
            alpha <- (1 - conf.level)/2
            CINT<-c(p.L(x, alpha), p.U(x, alpha))
        } else {
            CINT<-exactbinomCI(x,n,method=tsmethod,conf.level=conf.level,
                tol=tol,pRange=pRange)
        }
    }
 
    attr(CINT,"conf.level")<-conf.level

    ESTIMATE <- x/n
    names(x) <- "number of successes"
    names(n) <- "number of trials"
    names(ESTIMATE) <- names(p) <- "probability of success"
    methodphrase<-"Exact one-sided binomial test"
    if (alternative=="two.sided")
        methodphrase<-switch(tsmethod,
           minlike="Exact two-sided binomial test (sum of minimum likelihood method)",
           central="Exact two-sided binomial test (central method)",
           blaker="Exact two-sided binomial test (Blaker's method)")
           
    structure(list(statistic = x, parameter = n, p.value = PVAL, 
        conf.int = CINT, estimate = ESTIMATE, null.value = p, 
        alternative = alternative, method = methodphrase, 
        data.name = DNAME), class = "htest")
}
