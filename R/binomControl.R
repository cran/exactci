binomControl<-function(relErr=1+1e-07,tol=.00001,pRange=c(1e-10,1-1e-10)){
    if (relErr<=1) stop("relErr should be greater than 1")
    if (tol<=0) stop("tol should be greater than 0")
    if (pRange[1]>pRange[2] | pRange[1]<0 | pRange[1]>1 | pRange[2]<0 | pRange[2]>1) stop("pRange should be a range of possible values for the binomial parameter p")
    list(relErr=relErr,tol=tol,pRange=pRange)
}