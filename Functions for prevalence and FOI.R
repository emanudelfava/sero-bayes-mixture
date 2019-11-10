#Various functions to be used in all files

#################################################################
##### FUNCTION TO FIND THE BEST CONSTRAINED FP(M=1) #####
#################################################################

sel_foi_1<-function(data){
#Function to compute all fractional polynomials
p<-seq(-2,3,0.05)
p1<-numeric(1)
y<-data$x
n<-data$n
age<-data$age
#lmax<-sum(dbinom(y,n,y/n,log = TRUE)) # log-likelihood of saturared model
tab_dev<-matrix(0,length(p),2)
tab_dev[,1]<-p
colnames(tab_dev)<-c("p1","SDev")

for(i in 1:dim(tab_dev)[1]){
		X<-bfp(age,powers=tab_dev[i,1])
		fp1<-glm(cbind(y,n-y)~X,data=data,family=binomial(link="logit"))
		#lmod<-logLik(fp1) # log-likelihood of fitted model
		tab_dev[i,2]<-deviance(fp1) # scaled deviance (diff. saturated vs. fitted model)
}

#Function to select the best monotone fracpoly
tab_dev_ord<-tab_dev[order(tab_dev[,2]),]

check<-1
i<-1
while(check>0){
  p1<-tab_dev_ord[i,1]

#Estimation of the force of infection for a FP(m=1)-logit
  X1<-bfp(age,p1)
  final1<-glm(cbind(y,n-y)~X1,data=data,family=binomial(link="logit"))

  j.1<- function(t,a,b,p1) {
    if(p1!=0)
      expr<-expression(a+b*t^(p1))
    if(p1==0)
      expr<-expression(a+b*log(t))	
    c(eval(D(expr, "t")))
    }

  b.1<-list(a=coef(final1)[1],b=coef(final1)[2])
  t<-age
  prev<-fitted(final1)
  foi<-prev*do.call("j.1",c(list(t=t),b.1,p1))
    
# Check for FOI strictly positive
  if(all(foi>5e-3)) check<-0
    else{
      check<-1
      i<-i+1
    }
}

z<-list(p=p1,foi=foi,prev=prev,lmod=logLik(final1),
        dev=deviance(final1),b=c(coef(final1)))
z
}

#################################################################
##### FUNCTION TO FIND THE BEST CONSTRAINED FP(M=2) #####
#################################################################

sel_foi_2<-function(data){
#Function to compute all fractional polynomials
p<-seq(-2,3,0.05)
p1<-numeric(1)
p2<-numeric(1)
y<-data$x
n<-data$n
age<-data$age
#lmax<-sum(dbinom(y,n,y/n,log = TRUE)) # log-likelihood of saturared model
tab_dev<-matrix(0,length(p)*length(p),3)
tab_dev[,1]<-rep(p,each=length(p))
tab_dev[,2]<-rep(p,length(p))
colnames(tab_dev)<-c("p1","p2","SDev")

for(i in 1:dim(tab_dev)[1]){
  X<-bfp(age,powers=c(tab_dev[i,1],tab_dev[i,2]))
  fp2<-glm(cbind(y,n-y)~X,data=data,family=binomial(link="logit"))
  #lmod<-logLik(fp2) # log-likelihood of fitted model
  tab_dev[i,3]<-deviance(fp2) # scaled deviance
}

#Function to select the best monotone fracpoly
tab_dev_ord<-tab_dev[order(tab_dev[,3]),]

check<-1
i<-1
while(check>0){
  p1<-tab_dev_ord[i,1]
  p2<-tab_dev_ord[i,2]
  
  #Estimation of the force of infection for a FP(m=2)-logit
  X2<-bfp(age,c(p1,p2))
  final2<-glm(cbind(y,n-y)~X2,data=data,family=binomial(link="logit"))
  
  j.2<- function(t,a,b,c,p1,p2) {
    if(p1!=0&p2==p1)
      expr<-expression(a+b*t^(p1)+c*t^(p2)*log(t))
    if(p1!=0&p2!=0&p2!=p1)
      expr<-expression(a+b*t^(p1)+c*t^(p2))
    if(p1==0&p2==0)
      expr<-expression(a+b*log(t)+c*log(t)*log(t))
    if(p1==0&p2!=0)
      expr<-expression(a+b*log(t)+c*t^(p2))
    if(p1!=0&p2==0)
      expr<-expression(a+b*t^(p1)+c*log(t))
    c(eval(D(expr, "t")))
  }
  
  b.2<-list(a=coef(final2)[1],b=coef(final2)[2],c=coef(final2)[3])
  t<-age
  prev<-fitted(final2)
  foi<-prev*do.call("j.2",c(list(t=t),b.2,p1,p2))
  
  # Check for FOI strictly positive
  if(all(foi>5e-3)) check<-0
  else{
    check<-1
    i<-i+1
  }
}

z<-list(p=c(p1,p2),foi=foi,prev=prev,lmod=logLik(final2),
        dev=deviance(final2),b=c(coef(final2)))
z
}

#############################################
##### "PAVIT" FUNCTION FOR MONOTONICITY #####
#############################################

# PAVA algorithm (Hens et al. 2012)
pavit<- function(datai){
pai1<-pai2<-datai
N<-length(pai1)
for(i in 1:(N-1)){
	if (pai2[i] > pai2[i+1]){
     pool<-(pai1[i]+pai1[i+1])/2
     pai2[i:(i+1)]<-pool
     k<-i+1	
       for(j in (k-1):1){
       	if (pai2[j]>pai2[k]){
	       pool.2<-sum(pai1[j:k])/length(pai1[j:k])
	       pai2[j:k]<-pool.2}   
}
}
}	
return(list(pai1=pai1,pai2=pai2))
}

#################################################################
##### FUNCTION TO CONSTRAIN NEGATIVE FOI TO 0 #####
#################################################################

foi0<-function(obj){
	x<-fitted(obj)
	for(i in 1:length(x)){
		if(x[i]<0) x[i]<-0
		}
	return(x)
	}

#################################################################
##### FUNCTION FOR THE PIECEWISE CONSTANT FOI MODEL #####
#################################################################
	
# Function to obtain piecewise force of infection. 
# The breaks should reflect the school system in the specific country.
const_foi<-function(lambda,pos,tot,age){
	S<-numeric(class)
	foi<-numeric(max(a))
	for(i in 1:class){
		for(j in a[i]:(a[i+1]-1)){
			foi[j]<-lambda[i]
		}
	}
	S[1]<-exp(-foi[1])
	for(i in 2:max(age)){
		S[i]<-exp(-foi[i])*S[i-1]
	}
	R<-1-S[age]
	lmod<-sum(dbinom(pos,tot,R,log=TRUE))
	lmax<-sum(dbinom(pos,tot,pos/tot,log=TRUE))
	deviance<-2*(lmax-lmod)
	return(deviance)
}
		
#################################################################
#### FUNCTION FOR THE PREV FROM PIECEWISE CONSTANT FOI MODEL ####
#################################################################
	
#Function to obtain the prevalence, given the estimate for the piecewise force of infection
F_const_foi<-function(lambda,data){
	pos<-data$x
	tot<-data$n
	age<-data$age
	S<-numeric(class)
	foi<-numeric(max(a))
	S[1]<-1
	for(i in 1:class){
		for(j in a[i]:(a[i+1]-1)){
			foi[j]<-lambda[i]
		}
	}
	S[1]<-exp(-foi[1])
	for(i in 2:max(age)){
		S[i]<-exp(-foi[i])*S[i-1]
	}
	R<-1-S[age]
	foi<-foi[age]
	list(prev=R,lambda=foi)
	}

#### BOOTSTRAP ####

#Parametric bootstrap for q and R0
fp1.boot.fun<-function(dat){
	fit<-sel_foi_1(dat)
	m<-1
	res.fp<-waifw.fp(fit,m)
	R0<-res.fp$R0
	R0
}

vzv.rg<-function(data,mle){
#Function to generate random binomial numbers of infections.
#data<-data.frame(y,n)
	out<-data
	for(i in 1:nrow(out)){
		out$x[i]<-rbinom(1,out$n[i],mle[i])
		}
	return(out)
}

##### INFORMATIVE PRIOR DISTRIBUTIONS #####

## Compute the value of parameters (alpha, beta) for a Beta distribution to have mean and sd (m,s)
BetaPar <- function(m,s) {
  ss <- s^2
  alpha <- ((1 - m) / ss - 1 / m) * m ^ 2
  beta <- alpha * (1 / m - 1)
  list(alpha = alpha, beta = beta)
}

##### Variance of the ratio estimator #####
VarRatio<-function(x,y){
  sx<-sd(x)
  sy<-sd(y)
  s2x<-var(x)
  s2y<-var(y)
  mx<-mean(x)
  my<-mean(y)
  m2x<-mx^2
  m2y<-my^2
  x2<-x^2
  y2<-y^2
  
  var<-(m2y/m2x)*((s2y/m2y)-(2*cov(y,x)/(my*mx))+(s2x/m2x))
  c(var)
}

##
##  PURPOSE:  Function to calculate (posterior) summary
##            statistics for a difference of two quantities supplied as (MCMC) samples
##            * primarily used to calculate posterior summary for the difference of the
##              deviances of two competing models
##
##  PROGRAMMER:  Arnost Komarek
##
##  CREATED:     19/05/2010
##
## ===================================================

## \item{x}{a numeric vector with the sample of the first quantity}
## \item{y}{a numeric vector with the sample of the second quantity to be subtracted from \code{x}}
## \item{prob}{a numeric vector of probabilities for quantiles to be calculated from the sample of differences}
## \item{cut}{numeric value(s) which specify the cutoff(s) we are interested in estimating
##   \eqn{\mbox{P}(x - y < \mbox{cut})}{P(x - y < cut)} from the sample. The default values are motivated by
##   the arguments given in Section 4 of Aitkin, Liu and Chadwick (2009).}
## \item{na.rm}{logical indicating on how to handle \code{NA}'s}

summaryDiff <- function(x, y, prob=c(0.025, 0.5, 0.975), cut=c(-2*log(9), 0), na.rm=TRUE)
{
  x <- as.numeric(x)
  y <- as.numeric(y)
  if (length(x) != length(y)) stop("x and y must be of the same length")
  
  if (any(prob < 0)) stop("all prob values must be nonnegative")
  if (any(prob > 1)) stop("all prob values must not exceed 1")  
  
  diff <- x - y
  
  Ediff <- mean(diff[!(diff == Inf | diff == -Inf)], na.rm=na.rm)
  Qdiff <- quantile(diff, prob=prob, na.rm=na.rm)
  Pcut <- numeric(length(cut))
  for (i in 1:length(cut)) Pcut[i] <- prop.table(table(diff < cut[i]))["TRUE"]
  names(Pcut) <- paste("P(diff < ", round(cut, 2), ")", sep="")
  summ <- c(Ediff, Qdiff)
  names(summ)[1] <- "Mean"
  
  RET <- list(summary=summ, Pcut=Pcut)
  return(RET)    
}  
