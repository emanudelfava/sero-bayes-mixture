#####################################################################################
###################### MIXTURE MODELLING FOR VZV IN ENGLAND AND WALES ###############
#####################################################################################
rm(list=ls())
library(R2jags)
library(ggplot2)
library(sn)
library(MCMCpack)
library(R.utils)
library(mgcv)
library(snow)
library(dclone)
library(survey)
library(gplots)
source("codes/Functions for prevalence and FOI.R")
options(max.print=100000)

EnglandWales<-read.table("data/UK_polymod.txt",header=T,sep="\t",
                    na.strings = c("",NA))
EnglandWales<-EnglandWales[is.na(EnglandWales$vzvmiuml)==F,]
head(EnglandWales) # N=2091 # Those with no missing values for VZV antibodies
# These are the same data analysed by Vyse et al. (2004) with mixture models
# Collected in 1996.
summary(EnglandWales$age)

### Building design weights by age group
EnglandWales$agestrata<-cut(EnglandWales$age,
                      c(1:10,15,20),
                      right=F,include.lowest=T)
table(EnglandWales$agestrata)

EnglandWales$agesample<-cut(EnglandWales$age,
                      c(1:21),
                      right=F,include.lowest=T)
table(EnglandWales$agesample)

### Creating age groups for analysis (aggregating data from 60 years onward)
EnglandWales<-EnglandWales[order(EnglandWales$agesample),]
EnglandWales$agep<-rep(c(1:20),table(EnglandWales$agesample))
prop.table(table(EnglandWales$vzvres)) # overall prevalence from cut-off point
prop.table(table(EnglandWales$vzvres,EnglandWales$agesample),2) # age-specific prevalence from cut-off point
EnglandWales$Y<-log10(EnglandWales$vzvmiuml)

head(EnglandWales)

# 2.8% of cases are equivocal
# 67.8% prevalence of VZV

EnglandWales1<-subset(EnglandWales,age!=0) #I removed the equivocal cases and age 0
N<-dim(EnglandWales1)[1] # N=2091
#EnglandWales1.dw<-svydesign(ids=~id,data=EnglandWales1,weights=~dweight) # weighted sample
EnglandWales1.uw<-svydesign(ids=~id,data=EnglandWales1,weights=~1) # unweighted sample

##### FIT WITH SIR MODEL #####

### Other data
logod<-EnglandWales1$Y
Nsub<-length(logod)
age.g<-EnglandWales1$agep
age1<-unique(age.g)
Nage<-length(age1)
a<-c(1,5,10,15,max(EnglandWales1$age)+1) # age breaks
#a1<-c(1,5,10,15,21) # age breaks
J<-length(a)-1
h1<-diff(a)
h2<-rep(1,20)
eps<-0 # duration of maternal antibody protection (in years)
nu<-1-0 # 1- (PB prevalence estimate of VZV among mothers at age 28.5) (US Census 2014)
D<-7 # duration of infectiousness
lag<-0
vzv.data<-list("logod"=logod,"Nsub"=Nsub,"age.g"=age.g,"Nage"=Nage,"lag"=lag,
               "D"=D,"a"=a,"J"=J,"h1"=h1,"h2"=h2,"nu"=nu,"eps"=eps)

hist(logod,breaks=100,xlab=expression("IgG antibodies to VZV ("*log[10]*"[OD+1])"),
     freq=F,cex.lab=1.5, cex.main=2, cex.sub=2, mgp=c(2.5,1,0),
     cex.axis=1.75,col=2,main="")

#Parameters
vzv.params<-c("f","xi","omega","gamma","T","prev","foi","ppp")

#Initial values
vzv.inits<-function(){
  list("tau"=c(11,6),"xi0"=c(0.3,2))}

#Fit JAGS
date()
#cl<-makeCluster(3,type="SOCK")
vzv.pcfoi.sn.orig.new<-jags(vzv.data,vzv.inits,vzv.params,n.chains=3,
                                 model.file="codes/pcfoi_sir_sn_final.txt",n.iter=20000,n.thin = 10)
#vzv.pcfoi.sn.orig.new.up<-update(vzv.pcfoi.sn.orig.new,n.iter=10000)
#stopCluster(cl)
date()
print(vzv.pcfoi.sn.orig.new,digits=2)
save(vzv.pcfoi.sn.orig.new,file="results/vzv_pcfoi_sir_sn_new.RData")
load("results/vzv_pcfoi_sir_sn_new.RData")
recompile(vzv.pcfoi.sn.orig.new)
dic.samples(vzv.pcfoi.sn.orig.new$model,type="popt",n.iter=1000) 
# PED = 10,617 MSIR
# PED = 9,538 SIR

##### ANALYSIS OF RESULTS OF  JAGS MODEL #####
fit<-vzv.pcfoi.sn.orig.new

#Checking convergence
traceplot(fit,varname="xi")
xi1<-mcmc(fit$BUGSoutput$sims.matrix[,"xi[1]"])
geweke.diag(xi1)
heidel.diag(xi1)
xi1<-mcmc.list(mcmc(fit$BUGSoutput$sims.matrix[1:1000,"xi[1]"]),
               mcmc(fit$BUGSoutput$sims.matrix[1001:2000,"xi[1]"]),
               mcmc(fit$BUGSoutput$sims.matrix[2001:3000,"xi[1]"]))
gelman.diag(xi1)
plot(xi1)

xi2<-mcmc(fit$BUGSoutput$sims.matrix[,"xi[2]"])
geweke.diag(xi2)
heidel.diag(xi2)
xi2<-mcmc.list(mcmc(fit$BUGSoutput$sims.matrix[1:1000,"xi[2]"]),
               mcmc(fit$BUGSoutput$sims.matrix[1001:2000,"xi[2]"]),
               mcmc(fit$BUGSoutput$sims.matrix[2001:3000,"xi[2]"]))
gelman.diag(xi2)
plot(xi2)

traceplot(fit,varname="omega")
omega1<-mcmc(fit$BUGSoutput$sims.matrix[,"omega[1]"])
geweke.diag(omega1)
heidel.diag(omega1)
omega1<-mcmc.list(mcmc(fit$BUGSoutput$sims.matrix[1:1000,"omega[1]"]),
                  mcmc(fit$BUGSoutput$sims.matrix[1001:2000,"omega[1]"]),
                  mcmc(fit$BUGSoutput$sims.matrix[2001:3000,"omega[1]"]))
gelman.diag(omega1)
plot(omega1)

omega2<-mcmc(fit$BUGSoutput$sims.matrix[,"omega[2]"])
geweke.diag(omega2)
heidel.diag(omega2)
omega2<-mcmc.list(mcmc(fit$BUGSoutput$sims.matrix[1:1000,"omega[2]"]),
                  mcmc(fit$BUGSoutput$sims.matrix[1001:2000,"omega[2]"]),
                  mcmc(fit$BUGSoutput$sims.matrix[2001:3000,"omega[2]"]))
gelman.diag(omega2)
plot(omega2)

traceplot(fit,varname="gamma")
gamma1<-mcmc(fit$BUGSoutput$sims.matrix[,"gamma[1]"])
geweke.diag(gamma1)
heidel.diag(gamma1)
gamma1<-mcmc.list(mcmc(fit$BUGSoutput$sims.matrix[1:1000,"gamma[1]"]),
                  mcmc(fit$BUGSoutput$sims.matrix[1001:2000,"gamma[1]"]),
                  mcmc(fit$BUGSoutput$sims.matrix[2001:3000,"gamma[1]"]))
gelman.diag(gamma1)
plot(gamma1)

gamma2<-mcmc(fit$BUGSoutput$sims.matrix[,"gamma[2]"])
geweke.diag(gamma2)
heidel.diag(gamma2)
gamma2<-mcmc.list(mcmc(fit$BUGSoutput$sims.matrix[1:1000,"gamma[2]"]),
                  mcmc(fit$BUGSoutput$sims.matrix[1001:2000,"gamma[2]"]),
                  mcmc(fit$BUGSoutput$sims.matrix[2001:3000,"gamma[2]"]))
gelman.diag(gamma2)
plot(gamma2)

f<-mcmc(fit$BUGSoutput$sims.matrix[,"f"])
geweke.diag(f)
heidel.diag(f)
f<-mcmc.list(mcmc(fit$BUGSoutput$sims.matrix[1:1000,"f"]),
             mcmc(fit$BUGSoutput$sims.matrix[1001:2000,"f"]),
             mcmc(fit$BUGSoutput$sims.matrix[2001:3000,"f"]))
gelman.diag(f)
plot(f)

#Saving main parameters from piecewise-constant FOI model
Z<-fit$BUGSoutput$median$T
lb.Z<-fit$BUGSoutput$summary[1:2091,"2.5%"]
ub.Z<-fit$BUGSoutput$summary[1:2091,"97.5%"]
f<-fit$BUGSoutput$mean$f
lb.f<-fit$BUGSoutput$summary[2093,"2.5%"]
ub.f<-fit$BUGSoutput$summary[2093,"97.5%"]
foi<-fit$BUGSoutput$mean$foi
lb.foi<-fit$BUGSoutput$summary[2094:2113,"2.5%"]
ub.foi<-fit$BUGSoutput$summary[2094:2113,"97.5%"]
alpha<-fit$BUGSoutput$mean$gamma
lb.alpha<-fit$BUGSoutput$summary[2114:2115,"2.5%"]
ub.alpha<-fit$BUGSoutput$summary[2114:2115,"97.5%"]
sigma<-fit$BUGSoutput$mean$omega
lb.sigma<-fit$BUGSoutput$summary[2116:2117,"2.5%"]
ub.sigma<-fit$BUGSoutput$summary[2116:2117,"97.5%"]
prev<-fit$BUGSoutput$mean$prev
lb.prev<-fit$BUGSoutput$summary[2119:2138,"2.5%"]
ub.prev<-fit$BUGSoutput$summary[2119:2138,"97.5%"]
mu<-fit$BUGSoutput$mean$xi
lb.mu<-fit$BUGSoutput$summary[2139:2140,"2.5%"]
ub.mu<-fit$BUGSoutput$summary[2139:2140,"97.5%"]

EnglandWales1$mix.res<-Z
head(EnglandWales1)

#EnglandWales2.dw<-svydesign(ids=~id,data=EnglandWales1,weights=~dweight) # weighted sample
EnglandWales2.uw<-svydesign(ids=~id,data=EnglandWales1,weights=~1) # unweighted sample
#prop.table(svytable(~mix.res+sex,subset(EnglandWales2.dw,age>=20 & age<46)),1)
#prop.table(table(EnglandWales1[EnglandWales1$age>=20 & EnglandWales1$age<46,]$mix.res,
#                 EnglandWales1[EnglandWales1$age>=20 & EnglandWales1$age<46,]$sex),1)
#prop.table(table(EnglandWales1[EnglandWales1$age>=60,]$mix.res))
# 63.2% of susceptible people in the childbearing period are women.

x<-seq(min(EnglandWales1$Y,na.rm=T),max(EnglandWales1$Y,na.rm=T),length.out = 100)
mix.sn<-(1-mean(prev))*dsn(x,mu[1],sigma[1],alpha[1])+
  mean(prev)*dsn(x,mu[2],sigma[2],alpha[2])
  
hist(EnglandWales1$Y,breaks=100,xlab=expression("IgG antibodies to VZV ("*log[10]*"[OD+1])"),
     freq=F,cex.lab=1.5, cex.main=2, cex.sub=2, mgp=c(2.5,1,0),
     cex.axis=1.75,col=2,main="",ylim=c(0,2))
abline(v=log10(15+1),lty=2,lwd=2)
abline(v=log10(20+1),lty=2,lwd=2)
lines(x,dsn(x,mu[1],sigma[1],alpha[1]),col=1,lty=3,lwd=2)
lines(x,dsn(x,mu[2],sigma[2],alpha[2]),col=1,lty=3,lwd=2)
lines(x,mix.sn,col=4,lwd=3)

#Create a table with n° trials and n° successes
prop.table(table(Z)) # POS=72.5% (CI: 65.9% - 81.3%)
prop.table(table(round(lb.Z)))
prop.table(table(round(ub.Z)))

#Create a table with n° trials and n° successes, by age
n.mix<-table(EnglandWales1$agep) 
prop.mix<-prop.table(table(Z,EnglandWales1$agep),2)[2,] # Seropositive proportions
prop.mix.lb<-prop.table(table(round(lb.Z),EnglandWales1$agep),2)[2,]
prop.mix.ub<-prop.table(table(round(ub.Z),EnglandWales1$agep),2)[2,]

#DataGigi<-data.frame(age1,prop=prop.mix,lb=prop.mix.lb,ub=prop.mix.ub)
#write.csv(DataGigi,file="tmp/DataSeropositivePropFromMixture.csv",row.names = F)

plot(age1,prop.mix,ty="b",ylim=c(0,1),xlim=c(0,20))
lines(age1,prop.mix.lb,lty=2)
lines(age1,prop.mix.ub,lty=2)
polygon(c(age1,rev(age1)),c(prop.mix.ub,rev(prop.mix.lb)),
        col=rgb(0, 0, 0, 0.5), border=NA)

EnglandWales1$agegrp2<-cut(EnglandWales1$age,
                     c(1:9,10,15,20),
                     right=F,include.lowest=T)
#EnglandWales1<-EnglandWales1[order(EnglandWales1$agegrp2),]
#EnglandWales2.dw<-svydesign(ids=~id,data=EnglandWales1,weights=~dweight) # weighted sample
EnglandWales2.uw<-svydesign(ids=~id,data=EnglandWales1,weights=~1) # weighted sample

#pos.cut<-svyby(~mres==1,~agegrp2,EnglandWales2.dw,svytotal)[,3]
#tot.cut<-svyby(~mres!=2,~agegrp2,EnglandWales2.dw,svytotal)[,3]
pos.cut.uw<-svyby(~vzvres=="pos",~agegrp2,EnglandWales2.uw,svytotal,na.rm = T)[,3]
tot.cut.uw<-svyby(~vzvres!="equi",~agegrp2,EnglandWales2.uw,svytotal,na.rm = T)[,3]
pos.mix<-tapply(Z==1,EnglandWales1$agegrp2,sum)
pos.mix.lb<-tapply(lb.Z==1,EnglandWales1$agegrp2,sum)
pos.mix.ub<-tapply(ub.Z==1,EnglandWales1$agegrp2,sum)
tot.mix<-tapply(Z!=2,EnglandWales1$agegrp2,sum)
F.mix<-pos.mix/tot.mix
F.mix.lb<-pos.mix.lb/tot.mix
F.mix.ub<-pos.mix.ub/tot.mix
#F.cut<-pos.cut/tot.cut
F.cut.uw<-pos.cut.uw/tot.cut.uw
#F.GR<-c(0.112,0.163,0.402,0.485,0.653,0.698,.714,0.828,.781,.814,.895,.864,.90,.919,.957,.918,.948,.959,.958,.937)
age_point<-1:11

agelabel<-c("1y","2y","3y","4y","5y","6y","7y","8y","9y","10-14y","15-20y")

svg("figs/EnglandWales/Fig2_ew-pcfoi.svg",width=8.,height=5.)
par(mgp=c(3.5,1,0), mar=c(6,6,0,1),cex.lab=1.5,cex.axis=1.5)
fig1.bar<-barplot(F.cut.uw,col="gray",ylim=c(0,1.1),
                  names.arg="",las=3,axes=F,
                  ylab="Percentage seropositive",
                  cex.names = 1.5)
axis(1,at=fig1.bar,labels=agelabel,las=2)
axis(2,at=seq(0,1,0.1),labels=seq(0,100,10),las=2)
lines(fig1.bar,F.mix)
points(fig1.bar,F.mix,pch=19)
polygon(c(fig1.bar,rev(fig1.bar)),c(F.mix.ub,rev(F.mix.lb)),
        col=rgb(0, 0, 0, 0.5), border=NA)
dev.off()

agelabel1<-seq(0,20,5)
svg("figs/EnglandWales/Fig3_ew-pcfoi.svg",width=5.,height=5.)
#tiff("figs/Fig3.tif",width=480,height=480,compression = "lzw")
# Plot prevalence
par(mgp=c(4,1,0), mar=c(5.5,5.5,0,1),cex.lab=1.5,cex.axis=1.5)
plot(age1,prev,main="",lwd=3,
     xlab="Age (years)",ylab="Seroprevalence (%)",
     ylim=c(0,1.),xlim=c(0,20),ty="l",col="black",
     axes=F)
lines(age1,lb.prev,lty=2)
lines(age1,ub.prev,lty=2)
polygon(c(age1,rev(age1)),c(ub.prev,rev(lb.prev)),
        col=rgb(0, 0, 0, 0.5), border=NA)
points(age1,prop.mix,cex=0.01*n.mix)
axis(side=1,at=seq(0,20,5),agelabel1,las=3)
axis(side=2,at=seq(0,1.,0.1),seq(0,100,10),las=2)  
dev.off()

svg("figs/EnglandWales/Fig4_ew-pcfoi.svg",width=5.,height=5.)
# Plot force of infection
par(mgp=c(4,1,0), mar=c(5.5,5.5,0,1),cex.lab=1.5,cex.axis=1.5)
plot(age1[-(Nage-1):-Nage],foi[-(Nage-1):-Nage],main="",lwd=3,
     xlab="Age (years)",ylab="FOI",
     ylim=c(0,1.),xlim=c(0,20),ty="l",col="black",
     axes=F)
lines(age1[(Nage-2):Nage],foi[(Nage-2):Nage],lty=3,lwd=2)
lines(age1,lb.foi,lty=2)
lines(age1,ub.foi,lty=2)
polygon(c(age1,rev(age1)),c(ub.foi,rev(lb.foi)),
        col=rgb(0, 0, 0, 0.5), border=NA)
axis(side=1,at=seq(0,20,5),agelabel1,las=3)
axis(side=2,at=seq(0,1.,0.1),las=2)  
dev.off()

##### Comparison between models #####
load("results/vzv_nonparam_sn_new.RData")
dev.np<-vzv.prodbeta.sn.orig.new$BUGSoutput$sims.array[,,"deviance"]
load("results/vzv_pcfoi_sir_sn_new.RData")
dev.pc<-vzv.pcfoi.sn.orig.new$BUGSoutput$sims.array[,,"deviance"]

summaryDiff(dev.np,dev.pc) # the two models fit similarly

svg("figs/diffdev_vzv.svg",width=5.,height=5.)
par(mgp=c(3.5,1,0), mar=c(5,5,1.5,1.5),cex.lab=1.5,cex.axis=1.5,las=1)
plot(ecdf(dev.np),lwd=3,xlab="Posterior deviance",main="",axes=F)
lines(ecdf(dev.pc),lty=2,col=2,lwd=3)
legend("right",c("PB","SIR"),lwd=3,col=1:2,lty=1:2,bty="n",cex=1.5)
axis(side=1,at=seq(-3000,-1500,500))
axis(side=2,at=seq(0,1.,0.1)) 
dev.off()

##### Saving results for further analyses #####

### Mixture parameters
DataMix<-data.frame(status=c("susceptible","immune"),
                    mu,lb.mu,ub.mu,
                    sigma,lb.sigma,ub.sigma,
                    alpha,lb.alpha,ub.alpha,
                    f=c(f,NA),lb.f=c(ub.f,NA),ub.f=c(lb.f,NA),
                    country=rep("ew",2))
print(DataMix,digits=2)
write.csv(DataMix,file="results/DataMixtureParams_sir_vzv_ew.csv",row.names = F)

##### Plotting non-parametric and PC-FOI models together #####

n.mix<-as.numeric(n.mix)
PC<-data.frame(age1,prop.mix,n.mix,prev,lb.prev,ub.prev,foi,lb.foi,ub.foi)
write.csv(PC,file="results/DataPrevFOI_pcfoi_sir_vzv.csv",row.names=F)

# Reading data
BP<-read.csv("results/DataPrevFOI_nonparam_vzv.csv",header=T)
PC<-read.csv("results/DataPrevFOI_pcfoi_sir_vzv.csv",header=T)

agelabel1<-seq(0,20,5)
svg("figs/Fig2d.svg",width=5.,height=5.)
#tiff("figs/Fig3.tif",width=480,height=480,compression = "lzw")
# Plot prevalence
par(cex.lab=1.5, cex.main=1.5, cex.sub=1.5,cex.axis=1.5,mgp=c(3.5,1,0),mar=c(4.75,4.75,3,5.))
plot(BP$age1,BP$prev,main="England and Wales",lwd=3,
     xlab="Age (years)",ylab="Seroprevalence VZV (%)",
     ylim=c(0,1),xlim=c(0,20),ty="l",col="black",
     axes=F)
lines(BP$age1,BP$lb.prev,lty=2)
lines(BP$age1,BP$ub.prev,lty=2)
polygon(c(BP$age1,rev(BP$age1)),c(BP$ub.prev,rev(BP$lb.prev)),
        col=rgb(0, 0, 0, 0.25), border=NA)
points(BP$age1,BP$prop.mix,cex=0.02*BP$n.mix)

lines(PC$age1,PC$prev,lwd=3,lty=2,col=2)
lines(PC$age1,PC$lb.prev,lty=2,col=2)
lines(PC$age1,PC$ub.prev,lty=2,col=2)
polygon(c(PC$age1,rev(PC$age1)),c(PC$ub.prev,rev(PC$lb.prev)),
        col=rgb(1, 0, 0, 0.25), border=NA)
points(PC$age1,PC$prop.mix,cex=0.02*PC$n.mix,pch=20,col=2)

# Plot force of infection
Nage<-dim(BP)[1]
lines(BP$age1[-(Nage-1):-Nage],BP$foi[-(Nage-1):-Nage],lwd=3,
      ty="l",col="black")
lines(BP$age1[(Nage-2):Nage],BP$foi[(Nage-2):Nage],lty=3,lwd=2)
lines(BP$age1,BP$lb.foi,lty=2)
lines(BP$age1,BP$ub.foi,lty=2)
polygon(c(BP$age1,rev(BP$age1)),c(BP$ub.foi,rev(BP$lb.foi)),
        col=rgb(0, 0, 0, 0.25), border=NA)

Nage<-dim(PC)[1]
lines(PC$age1[-(Nage-1):-Nage],PC$foi[-(Nage-1):-Nage],lwd=3,lty=2,col=2)
lines(PC$age1[(Nage-2):Nage],PC$foi[(Nage-2):Nage],lty=2,lwd=2,col=2)
lines(PC$age1,PC$lb.foi,lty=2,col=2)
lines(PC$age1,PC$ub.foi,lty=2,col=2)
polygon(c(PC$age1,rev(PC$age1)),c(PC$ub.foi,rev(PC$lb.foi)),
        col=rgb(1, 0, 0, 0.25), border=NA)

axis(side=1,at=seq(0,20,5),agelabel1,las=3)
axis(side=2,at=seq(0,1,0.1),seq(0,100,10),las=2)
axis(side=4,at=seq(0,1,0.1),las=2)
mtext(text="FOI",side=4,cex=1.5,line=3.5,las=3)
legend("right",lty=1:2,col=1:2,lwd=3,c("PB","SIR"),bty="n",cex = 1.5)
dev.off()
