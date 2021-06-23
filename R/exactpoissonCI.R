exactpoissonCI<-function(x, tsmethod="minlike",conf.level=.95,tol=.00001,
    pRange=c(1e-10,1-1e-10)){
    if (conf.level<=0 | conf.level>=1) stop("conf.level must be 0<conf.level<1")
    # function to calculate theta tilde a,b
    # jumps for minlike p-value function
    ttheta<-function(a,b){ exp( (lfactorial(b) - lfactorial(a))/
                                  (b-a) ) }
    
    # function to calculate theta hat xlo, xhi
    # jumps for blaker p-value function
    htheta<-function(xlo,xhi, Tol=tol, PvalRange=pRange){
      ## root is the parameter where the 
      ## tails are equal, i.e.,
      ## where ppois(xlo,root)=
      ##      ppois(xhi-1,root,lower.tail=FALSE)
      ## we solve using uniroot
      rootfunc<-function(beta){
        ppois(xlo,beta) - 
          ppois(xhi-1,beta,lower.tail=FALSE)
      }
      ## note that for intergers,i, 
      ## ppois(i-1,x,lower.tail=FALSE)=pgamma(x,a)
      ## so we use qgamma to get range from pvalRange
      root<-uniroot(rootfunc,
                    c(qgamma(PvalRange[1],xlo),
                      qgamma(PvalRange[2],xhi)),tol=Tol)$root
      
      root
    }
    exactpoissonCI.minlike<-function(x, alpha=1-conf.level, Tol=tol){
      
      ##################
      # upper limit
      #################
      # no need to calculate p-value for ttheta(x,x+1), it is always 1
      j<-1
      pval.hi<-1 
      while (pval.hi>alpha){
        j<- j+1
        pval.hi<- exactpoissonPval(x,r=ttheta(x,x+j), tsmethod="minlike")
      }
      # exit the loop at the first pval.hi<=alpha
      # so go back one
      j<- j-1
      upperCL<- ttheta(x,x+j)
      ###################
      # lower limit
      ###################
      if (x==0){
        lowerCL<-0
      } else {
        # for each j, there is an interval
        # ttheta(x-j,x) <= theta0 < ttheta(x-j+1,x)
        # we need to calculate the p-value at each end 
        # of the interval
        # let those two p-value be 
        # pval.lowLeft and pval.lowRight
        #
        # for j=1
        # ttheta(x-1,x)=x
        # and pval.lowLeft=pm(x,x)=1 = pval.lowRight
        tthetaLeft<- x
        pval.lowLeft<- 1
        for (j in 2:(x+1)){
          # to get pval.lowRight, take the previous pval.lowLeft and subtract 
          #  f(x-j+1;theta=ttheta(x-j+1,x))
          # previous tthetaLeft becomes tthetaRight
          tthetaRight<- tthetaLeft
          pval.lowRight<- pval.lowLeft - dpois(x-j+1,tthetaRight)
          if (x-j<0){
            # if x-j<0, then F(x-j;theta0)=0 and we just need to 
            # find theta0 so that Fbar(x, theta0)=alpha
            # Fbar(x,B) =ppois(x-1, B, lower.tail=FALSE) = alpha
            # when  qgamma(alpha,x,1) = B
            # so we get the lower confidence limit 
            # by
            lowerCL<- qgamma(alpha, x, 1)
            break()
          } 
          tthetaLeft<- ttheta(x-j,x)
          pval.lowLeft<- exactpoissonPval(x,r=tthetaLeft, tsmethod="minlike")
          if (pval.lowLeft==alpha){
            lowerCL<- tthetaLeft
            break()
          } else if (pval.lowLeft<alpha & pval.lowRight>alpha){
            # assume p-value monotonically increasing and use uniroot
            # to find when pval=alpha
            rootfunc<-function(R, Alpha=alpha){
              exactpoissonPval(x,r=R, tsmethod="minlike") - Alpha
            }
            lowerCL<-uniroot(rootfunc, lower=tthetaLeft, upper=tthetaRight-Tol, tol=Tol)$root
            break()
          } else if (pval.lowLeft<alpha & pval.lowRight<=alpha){
            # assume p-value increases from pval.lowLeft to pval.lowRight
            # if pval.lowRight<=alpha, then tthetaRight+epsilon > alpha
            # and tthetaRight+epsilon= tthetaRight=lowerCL
            lowerCL<- tthetaRight
            break()
          } # else continue to next j
        }
        
      }
      
      
      c(lowerCL,upperCL)
    }
    

    exactpoissonCI.blaker<-function(x, alpha=1-conf.level, Tol=tol,pvalRange=c(1e-10,1-1e-10)){
      
      ##################
      # lower limit
      #################
      # 
      if (x==0){
        lowerCL<-0
      } else {
        # p-value function is piecewise continuous
        # for each interval let the range of the interval be 
        # hthetaLeft = htheta(x-j,x)   to 
        # hthetaRight = htheta(x-j+1,x) - epsilon
        # with p-values at those ends 
        # pval.loLeft ----pval.loRight
        # for the lower p-values (for theta0<x)
        # the p-values are strictly increasing in the 
        # piecewise continuous region 
        # (conjecture for now, but all plots show this)
        # at j=1 the left and right p-value is 1
        pval.loLeft<-pval.loRight<-1
        j<-1
        hthetaRight<- htheta(x-j+1,x)
        hthetaLeft<- htheta(x-j,x)
        
        
        for (j in 2:(x+1)){
          # calculate hthetaRight = previous j value's hthetaLeft
          hthetaRight<- hthetaLeft
          # take previous j values's pvalue and subtract off pmf
          # to get the pval.loRight
          pval.loRight<- pval.loLeft - dpois(x-j+1,hthetaRight)
          if (pval.loRight<=alpha){
            lowerCL<- hthetaRight
            break()
          }
          if (x-j<0){
            # when x-j<0 and pval.loRight>alpha, 
            # that means ppois(x-j,theta0)=0, so the solution
            # of the p-value equal to alpha is the value theta0 such that
            # Fbar(x,theta0) = ppois(x-1, theta0, lower.tail=FALSE) = alpha
            # So solve for theta0: 
            # ppois(x-1, theta0, lower.tail=FALSE) = alpha
            # Do that using relationship between Poisson
            # and gamma distributions:
            # ppois(x-1, theta0, lower.tail=FALSE) = alpha
            # when  qgamma(alpha,x,1) = theta0
            # so we get the lower confidence limit 
            # by
            lowerCL<- qgamma(alpha, x, 1)
            break()
          } else {
            hthetaLeft<- htheta(x-j,x)
            pval.loLeft<- ppois(x-j, hthetaLeft) + ppois(x-1,hthetaLeft, lower.tail=FALSE)
            if (pval.loLeft<=alpha){
              # pval.loLeft<= alpha < pval.loRight
              # conjectur that piecewise p-value function is monotonic
              # use uniroot
              rootfunc<-function(theta0,Alpha=alpha){
                Alpha - ppois(x-j, theta0) - ppois(x-1,theta0, lower.tail=FALSE)
              }
              lowerCL<- uniroot(rootfunc, lower=hthetaLeft, upper=hthetaRight, tol=Tol)$root
              break()
            } 
          } # end x-j>=0 
        } # end for loop
      } # end: x>0
      ###################
      # upper limit
      ###################
      # for the piecewise continuous p-value functions on the upper 
      # part of the curve (when theta0>x)
      # let the two bounds on the interval be
      # hthetaLeft and hthetaRight
      # and the p-values approaching those ends from within the interval, be 
      # pval.hiLeft and pval.hiRight
      # an issue with the upper limits for the Blaker interval is that
      # the p-value function is not monotonic within these intervals
      # the p-value function is decreasing until some value, 
      # say pval.hiMid, then it increases until the end of the interval
      # We can solve for pval.hiMid by taking the derivative and 
      # setting that to zero and solving for theta0, and 
      # it turns out that is equal to one of the ends of the minlike interval
      # thus we have
      # hthetaLeft < hthetaMid < hthetaRight
      #    and  (when the p-values<1)
      # pval.hiLeft > pval.hiMid < pval.hiRight
      pval.hiMid<-pval.hiLeft<- 1
      pval.hiRight<-  1
      j<-1
      hthetaRight<- htheta(x,x+j)
      while (pval.hiMid>alpha){
        hthetaLeft<- hthetaRight
        pval.hiLeft<- pval.hiRight - dpois(x+j, hthetaLeft)
        hthetaMid<- ttheta(x,x+j)
        j<- j+1
        hthetaRight<- htheta(x,x+j)
        pval.hiMid<- ppois(x, hthetaMid) + ppois(x+j-1, hthetaMid, lower.tail=FALSE)
        pval.hiRight<- ppois(x, hthetaRight) + ppois(x+j-1, hthetaRight, lower.tail=FALSE)   
      }
      # after while loop, pval.hiMid<=alpha
      if (pval.hiMid==alpha){
        upperCL<- hthetaMid
      } else if (pval.hiLeft< alpha){
        # since 
        # pval.hiMid < pval.hiRight < pval.hiLeft
        # if pval.hiLeft< alpha
        upperCL<- hthetaLeft
      } else if (pval.hiRight==alpha){
        upperCL<- hthetaRight
      } else if (pval.hiLeft>alpha){
        if (pval.hiRight>alpha){
          # pval.hiMid < alpha < pval.hiRight
          upperCL<-hthetaRight
        } else {
          # pval.hiRight < alpha < pval.hiLeft
          rootfunc<-function(theta0, Alpha=alpha){
            Alpha -  ppois(x, theta0) - ppois(x+j-1, theta0, lower.tail=FALSE) 
          }
          upperCL<-uniroot(rootfunc, lower=hthetaLeft, upper=hthetaMid, tol=Tol)$root
        }
      } 
      c(lowerCL,upperCL)
    }
    
    ##########################################
    # end of function definitions
    ##########################################
    if (tsmethod=="minlike"){
      CINT<-exactpoissonCI.minlike(x)
    } else if (tsmethod=="blaker"){
      CINT<-exactpoissonCI.blaker(x)
    }
    CINT
}

#exactpoissonCI(5,tsmethod="minlike")
#exactpoissonCI(1,tsmethod="minlike",conf.level=.366)
#exactpoissonCI(1,tsmethod="blaker",conf.level=.366)

### Check answers for x=1
###
#theta<-40:1250/1000
#pvalb<-cib<-pvalm<-cim<-rep(NA,length(theta))
#for (i in 1:length(theta)){
# pvalm[i]<-exactpoissonPvals(1,theta[i],tsmethod="minlike")$pvals
# if (pvalm[i]<1) cim[i]<-exactpoissonCI(1,tsmethod="minlike",conf.level=1-pvalm[i])[1]
# pvalb[i]<-exactpoissonPvals(1,theta[i],tsmethod="blaker")$pvals
# if (pvalb[i]<1) cib[i]<-exactpoissonCI(1,tsmethod="blaker",conf.level=1-pvalb[i])[1]
#}
# plot lower Conf Limit (at conf.level=1-pval(theta)), 
# should get theta=lowerCL, when pval(theta)<1 
#par(mfrow=c(2,2))
#plot(theta,pvalm,main="minlike")
#plot(theta,pvalb,main="blaker")
#plot(theta,cim)
#plot(theta,cib)