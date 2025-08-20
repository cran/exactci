## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo=FALSE,
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(exactci)

## -----------------------------------------------------------------------------
# We compare our results to those from the function by Reiczigel
# referenced in Kaschka and Reiczigel (2020)
# here is that function copied without changes....
###################################################################
# beginning of  function by Reiczigel
####################################################################
poisson.Sterne.CI = function(x,conflevel=.95,epsilon=1e-8){
	# Sterne's CI for Poisson data
	# epsilon is used in interval halving and computing left-right limits
	# reiczigel.jeno@univet.hu (07.07.2019)

	alpha=1-conflevel
	lam = function(k,m) exp((lfactorial(m)-lfactorial(k))/(m-k))

	# first going downwards (except if x==0)
	if(x==0) LCL=0 else {
		i=x-1
		repeat{
			if(i<0) L=0 else L=lam(i,x)
			#pvr=sterne.poi.test.orig(x,L+epsilon)
			############################################
			lambda0=L+epsilon
			threshold=dpois(x,lambda0); sum.gt=0
			ii=x-1
			repeat{
				probab=dpois(ii,lambda0)
				if(probab > threshold) {
					sum.gt=sum.gt+probab
					ii=ii-1
				} else break
			}
			ii=x+1
			repeat{
				probab=dpois(ii,lambda0)
				if(probab > threshold) {
					sum.gt=sum.gt+probab
					ii=ii+1
				} else break
			}
			pvr=1-sum.gt
			############################################

			#cat(i,L,pvr,"\n")
			if(pvr>alpha){
				i=i-1 
				L.previous=L
				pvl.previous=pvr-dpois(x,L+epsilon)
			} else break
		}
		#cat("   ",L.previous,pvl.previous,"\n")
		if(pvl.previous<alpha) LCL=L.previous else {
			#interval halving
			left.end=L; right.end=L.previous
			repeat{
				midpoint=(left.end+right.end)/2
				 if(ppois(i,midpoint)+1-ppois(x-1,midpoint)>alpha){
					right.end=midpoint
				} else {
					left.end=midpoint
				}
				if(right.end-left.end < epsilon) break
			}
			LCL=left.end
		}
	}
	#cat("------\n")
	
	# now going upwards
	i=x+1
	repeat{
		L=lam(x,i)
		#pvl=sterne.poi.test.orig(x,L-epsilon)
		############################################
		lambda0=L-epsilon
		threshold=dpois(x,lambda0); sum.gt=0
		ii=x-1
		repeat{
			probab=dpois(ii,lambda0)
			if(probab > threshold) {
				sum.gt=sum.gt+probab
				ii=ii-1
			} else break
		}
		ii=x+1
		repeat{
			probab=dpois(ii,lambda0)
			if(probab > threshold) {
				sum.gt=sum.gt+probab
				ii=ii+1
			} else break
		}
		pvl=1-sum.gt
		############################################

		#cat(i,L,pvl,"\n")
		if(pvl>alpha){
			i=i+1
			L.previous=L
			pvr.previous=pvl-dpois(x,L-epsilon)
		} else break
	}
	
	#cat("   ",L.previous,pvr.previous,"\n")
	# my conjecture is that this part can be replaced simply by
	# UCL = lam(x,i-1)
	if(pvr.previous<alpha) UCL=L.previous else {
		#interval halving - I think never needed
		cat("!!! BIG SURPRISE !!! in 'sterne.poi.ci' \n")
		left.end=L.previous; right.end=L
		repeat{
			midpoint=(left.end+right.end)/2
			############################################
			lambda0=midpoint
			threshold=dpois(x,lambda0); sum.gt=0
			ii=x-1
			repeat{
				probab=dpois(ii,lambda0)
				if(probab > threshold) {
					sum.gt=sum.gt+probab
					ii=ii-1
				} else break
			}
			ii=x+1
			repeat{
				probab=dpois(ii,lambda0)
				if(probab > threshold) {
					sum.gt=sum.gt+probab
					ii=ii+1
				} else break
			}
			sterne.poi.orig=1-sum.gt
			############################################
			#if(sterne.poi.test.orig(x,midpoint)>alpha){
			
			if(sterne.poi.orig>alpha){
				left.end=midpoint
			} else {
				right.end=midpoint
			}
			if(right.end-left.end < epsilon) break
		}
		UCL=right.end
	}
	return(c(LCL,UCL))
}
###################################################################
# end of  function by Reiczigel
####################################################################

X<-c(1,0,8,8,8)
CL<-c(.70,.95,.99,.9,.01)



checkData<- matrix(NA,2*length(X),4,dimnames=list(rep(c("Reiczigel","Fay"),length(X)),
                                                c("x","conf level","lower CL","upper CL")))

checkData[,"x"]<- rep(X,each=2)
checkData[,"conf level"]<-rep(CL,each=2)
for (i in 1:length(X)){
  cir<-poisson.Sterne.CI(X[i], conflevel=CL[i])
  checkData[2*i-1,3:4]<- cir
  cif<-exactpoissonCI(X[i],tsmethod="minlike",conf.level=CL[i])
  checkData[2*i,3:4]<-cif  
}

library(knitr)
kable(checkData,caption="Compare Results")

## ----fig.width=8, fig.height=6------------------------------------------------
library(exactci)
exactpoissonPlot(8,tsmethod="central", dolog=FALSE, dolines = FALSE, dopoints=TRUE,pch=16,cex=1, col="blue")
exactpoissonPlot(8,tsmethod="blaker", dolines = FALSE, dopoints=TRUE, newplot=FALSE, col=gray(0.2), pch=16,cex=1.5)
exactpoissonPlot(8,tsmethod="minlike", dolines = FALSE, dopoints=TRUE, newplot=FALSE, col=gray(0.8), pch=16, cex=.75)
legend("topright",legend=c("central","minlike","blaker"),pch=c(16,16,16),col=c("blue",gray(.8),gray(.2)))

## ----fig.width=8, fig.height=6------------------------------------------------
par(mfrow=c(1,2))

exactpoissonPlot(8,tsmethod="central", dolines = FALSE, dopoints=TRUE,pch=16,cex=1, col="blue",xlim=c(3.5,5.5),ylim=c(.05,.11))
exactpoissonPlot(8,tsmethod="blaker", dolines = FALSE, dopoints=TRUE, newplot=FALSE, col=gray(0.2), pch=16,cex=1.5)
exactpoissonPlot(8,tsmethod="minlike", dolines = FALSE, dopoints=TRUE, newplot=FALSE, col=gray(0.8), pch=16, cex=.75)
legend("topright",legend=c("central","minlike","blaker"),pch=c(16,16,16),col=c("blue",gray(.8),gray(.2)))

exactpoissonPlot(8,tsmethod="central", dolines = FALSE, dopoints=TRUE,pch=16,cex=1, conf.level=1-0.0558, col="blue",xlim=c(15.4,15.6),ylim=c(.055,.058))
exactpoissonPlot(8,tsmethod="blaker", dolines = FALSE, dopoints=TRUE, newplot=FALSE, col=gray(0.2), pch=16,cex=1.5)
exactpoissonPlot(8,tsmethod="minlike", dolines = FALSE, dopoints=TRUE, newplot=FALSE, col=gray(0.8), pch=16, cex=.75)

legend("topleft",legend=c("central","minlike","blaker"),pch=c(16,16,16),col=c("blue",gray(.8),gray(.2)))


par(mfrow=c(1,1))


## -----------------------------------------------------------------------------

library(BlakerCI)

X<-c(0,1,8,8,8)
CL<-c(.70,.95,1-0.0558,1-0.056,1-0.057)



checkData<- matrix(NA,2*length(X),4,dimnames=list(rep(c("Klaschka","Fay"),length(X)),
                                                c("x","conf level","lower CL","upper CL")))

checkData[,"x"]<- rep(X,each=2)
checkData[,"conf level"]<-rep(CL,each=2)
for (i in 1:length(X)){
  cik<-poisson.blaker.limits(X[i], level=CL[i])
  checkData[2*i-1,3:4]<- cik
  cif<-exactpoissonCI(X[i],tsmethod="blaker", conf.level=CL[i])
  checkData[2*i,3:4]<-cif  
}

library(knitr)
kable(checkData,caption="Compare Results")

