#####################################################################################
###################### MIXTURE MODELLING FOR B19 IN BELGIUM ##########################
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
source("Functions for prevalence and FOI.R")
options(max.print=100000)

Dat<-read.csv("serodata_vzv_b19_be.csv",header=T,sep=",",na.strings = "")
Belgium<-Dat[is.na(Dat$parvouml)==F,c(1,2,3,4,5,8,9)]
head(Belgium) # N=3098 # Those with no missing values for b19 antibodies
range(Belgium$age)
# Data collected in 2001-2003.

### ESEN2 data for Belgium are a subset of data for POLYMOD
#Belgium<-read.table("original/Belgium_esen.txt",header=T,sep="\t",
#                    na.strings = c("",NA))
#dim(Belgium)
#Belgium<-Belgium[is.na(Belgium$stduni)==F,]
#head(Belgium) # N=2762 # Those with no missing values for b19 antibodies
#range(Belgium$age)
# Data collected in 2002.

### Building design weights by age group
Belgium$agestrata<-cut(Belgium$age,
                       c(0:10,15,20,25,30,35,40,45,50,60,70,83),
                       right=F,include.lowest=T)
table(Belgium$agestrata)

Belgium$agesample<-cut(Belgium$age,
                       c(0:60,83),
                       right=F,include.lowest=T)
table(Belgium$agesample)

### Creating age groups for analysis (aggregating data from 60 years onward)
Belgium<-Belgium[order(Belgium$agesample),]
Belgium$agep<-rep(c(0:60),table(Belgium$agesample))
prop.table(table(Belgium$parvores)) # overall prevalence from cut-off point
prop.serop<-prop.table(table(Belgium$parvores,Belgium$age),2)[2,] # age-specific prevalence from cut-off point
Belgium$Y<-log10(Belgium$parvouml+1)

plot(prop.serop)

head(Belgium)

# 3.8% of cases are equivocal
# 85.2% prevalence of b19

Belgium1<-subset(Belgium,age>0 & age<=65) #I removed the equivocal cases and age 0
N<-dim(Belgium1)[1] # N=2760
range(Belgium1$age)
#Belgium1.dw<-svydesign(ids=~id,data=Belgium1,weights=~dweight) # weighted sample
Belgium1.uw<-svydesign(ids=~yid,data=Belgium1,weights=~1) # unweighted sample

##### FIT WITH PIECEWISE-CONSTANT FOI MODEL #####

### Other data
logod<-Belgium1$Y
Nsub<-length(logod)
age.g<-Belgium1$agep
age1<-unique(age.g)
Nage<-length(age1)
a<-c(1,5,10,15,20,25,30,35,40,45,50,55,60,max(Belgium1$age)+1) # age breaks
J<-length(a)-1
h1<-diff(a)
h2<-c(rep(1,59),6) # size of one-year FOI age groups (plus final age group after age 60)
eps<-0 # duration of maternal antibody protection (in years)
nu<-1-0 # 1- (PB prevalence estimate of VZV among mothers at age 28.6) (US Census 2013)
D<-6 # duration of infectiousness
lag<-0
b19.data<-list("logod"=logod,"Nsub"=Nsub,"age.g"=age.g,"Nage"=Nage,"lag"=lag,
               "D"=D,"a"=a,"J"=J,"h1"=h1,"h2"=h2,"nu"=nu,"eps"=eps)

hist(logod,breaks=100,xlab=expression("IgG antibodies to B19 ("*log[10]*"[OD+1])"),
     freq=F,cex.lab=1.5, cex.main=2, cex.sub=2, mgp=c(2.5,1,0),
     cex.axis=1.75,col=2,main="")

#Parameters
b19.params<-c("f","xi","omega","gamma","T","prev","foi","ppp")
 
#Initial values
b19.inits<-function(){
  list("tau"=c(11.11,4.00),"xi0"=c(0.8,2.3))}

#Fit JAGS
date()
#cl<-makeCluster(3,type="SOCK")
b19.pcfoi.sn.orig.new<-jags(b19.data,b19.inits,b19.params,n.chains=2,
                            model.file="pcfoi_sir_sn_model.txt",n.iter=30000,n.thin = 10)
#b19.pcfoi.sn.orig.new.up<-update(b19.pcfoi.sn.orig.new,n.iter=10000,n.thin = 10)
#stopCluster(cl)
date()
print(b19.pcfoi.sn.orig.new,digits=2)
save(b19.pcfoi.sn.orig.new,file="b19_pcfoi_sir_sn_res.RData")
load("b19_pcfoi_sir_sn_res.RData")
recompile(b19.pcfoi.sn.orig.new)
dic.samples(b19.pcfoi.sn.orig.new$model,type="popt",n.iter=1000) 
# PED = 11,732

##### ANALYSIS OF RESULTS OF JAGS MODEL #####
fit<-b19.pcfoi.sn.orig.new

#Checking convergence
traceplot(fit,varname="xi")
xi1<-mcmc(fit$BUGSoutput$sims.matrix[,"xi[1]"])
geweke.diag(xi1)
heidel.diag(xi1)
xi1<-mcmc.list(mcmc(fit$BUGSoutput$sims.matrix[1:1000,"xi[1]"]),
               mcmc(fit$BUGSoutput$sims.matrix[1001:2000,"xi[1]"]))
gelman.diag(xi1)
plot(xi1)

xi2<-mcmc(fit$BUGSoutput$sims.matrix[,"xi[2]"])
geweke.diag(xi2)
heidel.diag(xi2)
xi2<-mcmc.list(mcmc(fit$BUGSoutput$sims.matrix[1:1000,"xi[2]"]),
               mcmc(fit$BUGSoutput$sims.matrix[1001:2000,"xi[2]"]))
gelman.diag(xi2)
plot(xi2)

traceplot(fit,varname="omega")
omega1<-mcmc(fit$BUGSoutput$sims.matrix[,"omega[1]"])
geweke.diag(omega1)
heidel.diag(omega1)
omega1<-mcmc.list(mcmc(fit$BUGSoutput$sims.matrix[1:1000,"omega[1]"]),
                  mcmc(fit$BUGSoutput$sims.matrix[1001:2000,"omega[1]"]))
gelman.diag(omega1)
plot(omega1)

omega2<-mcmc(fit$BUGSoutput$sims.matrix[,"omega[2]"])
geweke.diag(omega2)
heidel.diag(omega2)
omega2<-mcmc.list(mcmc(fit$BUGSoutput$sims.matrix[1:1000,"omega[2]"]),
                  mcmc(fit$BUGSoutput$sims.matrix[1001:2000,"omega[2]"]))
gelman.diag(omega2)
plot(omega2)

traceplot(fit,varname="gamma")
gamma1<-mcmc(fit$BUGSoutput$sims.matrix[,"gamma[1]"])
geweke.diag(gamma1)
heidel.diag(gamma1)
gamma1<-mcmc.list(mcmc(fit$BUGSoutput$sims.matrix[1:1000,"gamma[1]"]),
                  mcmc(fit$BUGSoutput$sims.matrix[1001:2000,"gamma[1]"]))
gelman.diag(gamma1)
plot(gamma1)

gamma2<-mcmc(fit$BUGSoutput$sims.matrix[,"gamma[2]"])
geweke.diag(gamma2)
heidel.diag(gamma2)
gamma2<-mcmc.list(mcmc(fit$BUGSoutput$sims.matrix[1:1000,"gamma[2]"]),
                  mcmc(fit$BUGSoutput$sims.matrix[1001:2000,"gamma[2]"]))
gelman.diag(gamma2)
plot(gamma2)

f<-mcmc(fit$BUGSoutput$sims.matrix[,"f"])
geweke.diag(f)
heidel.diag(f)
f<-mcmc.list(mcmc(fit$BUGSoutput$sims.matrix[1:1000,"f"]),
             mcmc(fit$BUGSoutput$sims.matrix[1001:2000,"f"]))
gelman.diag(f)
plot(f)

#Saving main parameters from piecewise-constant FOI model
Z<-fit$BUGSoutput$median$T
lb.Z<-fit$BUGSoutput$summary[1:3085,"2.5%"]
ub.Z<-fit$BUGSoutput$summary[1:3085,"97.5%"]
f<-fit$BUGSoutput$mean$f
lb.f<-fit$BUGSoutput$summary[3087,"2.5%"]
ub.f<-fit$BUGSoutput$summary[3087,"97.5%"]
foi<-fit$BUGSoutput$mean$foi[1:60]
lb.foi<-fit$BUGSoutput$summary[3088:3147,"2.5%"]
ub.foi<-fit$BUGSoutput$summary[3088:3147,"97.5%"]
alpha<-fit$BUGSoutput$mean$gamma
lb.alpha<-fit$BUGSoutput$summary[3153:3154,"2.5%"]
ub.alpha<-fit$BUGSoutput$summary[3153:3154,"97.5%"]
sigma<-fit$BUGSoutput$mean$omega
lb.sigma<-fit$BUGSoutput$summary[3155:3156,"2.5%"]
ub.sigma<-fit$BUGSoutput$summary[3155:3156,"97.5%"]
prev<-fit$BUGSoutput$mean$prev
lb.prev<-fit$BUGSoutput$summary[3158:3217,"2.5%"]
ub.prev<-fit$BUGSoutput$summary[3158:3217,"97.5%"]
mu<-fit$BUGSoutput$mean$xi
lb.mu<-fit$BUGSoutput$summary[3218:3219,"2.5%"]
ub.mu<-fit$BUGSoutput$summary[3218:3219,"97.5%"]

Belgium1$mix.res<-Z
head(Belgium1)

#Belgium2.dw<-svydesign(ids=~id,data=Belgium1,weights=~dweight) # weighted sample
Belgium2.uw<-svydesign(ids=~yid,data=Belgium1,weights=~1) # unweighted sample
#prop.table(svytable(~mix.res+sex,subset(Belgium2.dw,age>=20 & age<46)),1)
prop.table(table(Belgium1[Belgium1$age>=20 & Belgium1$age<46,]$mix.res,
                 Belgium1[Belgium1$age>=20 & Belgium1$age<46,]$gender),1)
prop.table(table(Belgium1[Belgium1$age>=60,]$mix.res))
# 48.5% of susceptible people in the childbearing period are women.

x<-seq(min(Belgium1$Y,na.rm=T),max(Belgium1$Y,na.rm=T),length.out = 100)
mix.sn<-(1-mean(prev))*dsn(x,mu[1],sigma[1],alpha[1])+
  mean(prev)*dsn(x,mu[2],sigma[2],alpha[2])
  
hist(Belgium1$Y,breaks=100,xlab=expression("IgG antibodies to b19 ("*log[10]*"[OD+1])"),
     freq=F,cex.lab=1.5, cex.main=2, cex.sub=2, mgp=c(2.5,1,0),
     cex.axis=1.75,col=2,main="",ylim=c(0,2))
abline(v=log10(20+1),lty=2,lwd=2)
abline(v=log10(24+1),lty=2,lwd=2)
lines(x,(1-mean(prev))*dsn(x,mu[1],sigma[1],alpha[1]),col=1,lty=3,lwd=2)
lines(x,mean(prev)*dsn(x,mu[2],sigma[2],alpha[2]),col=1,lty=3,lwd=2)
lines(x,mix.sn,col=4,lwd=3)

#Create a table with n° trials and n° successes
prop.table(table(round(Z))) # POS=67.4% (CI: 66.1% - 69.5%)
prop.table(table(round(lb.Z)))
prop.table(table(round(ub.Z)))

#Create a table with n° trials and n° successes, by age
n.mix<-table(Belgium1$agep) 
prop.mix<-prop.table(table(round(Z),Belgium1$agep),2)[2,] # Seropositive proportions
prop.mix.lb<-prop.table(table(round(lb.Z),Belgium1$agep),2)[2,]
prop.mix.ub<-prop.table(table(round(ub.Z),Belgium1$agep),2)[2,]

#DataGigi<-data.frame(age1,prop=prop.mix,lb=prop.mix.lb,ub=prop.mix.ub)
#write.csv(DataGigi,file="tmp/DataSeropositivePropFromMixture.csv",row.names = F)

plot(age1,prop.mix,ty="b",ylim=c(0,1))
lines(age1,prop.mix.lb,lty=2)
lines(age1,prop.mix.ub,lty=2)
polygon(c(age1,rev(age1)),c(prop.mix.ub,rev(prop.mix.lb)),
        col=rgb(0, 0, 0, 0.5), border=NA)

Belgium1$agegrp2<-cut(Belgium1$age,
                      c(1:10,15,20,25,30,35,40,45,50,60,70),
                     right=F,include.lowest=T)
#Belgium1<-Belgium1[order(Belgium1$agegrp2),]
#Belgium2.dw<-svydesign(ids=~id,data=Belgium1,weights=~dweight) # weighted sample
Belgium2.uw<-svydesign(ids=~yid,data=Belgium1,weights=~1) # weighted sample

#pos.cut<-svyby(~mres==1,~agegrp2,Belgium2.dw,svytotal)[,3]
#tot.cut<-svyby(~mres!=2,~agegrp2,Belgium2.dw,svytotal)[,3]
pos.cut.uw<-svyby(~parvores==1,~agegrp2,Belgium2.uw,svytotal,na.rm = T)[,3]
tot.cut.uw<-svyby(~is.na(parvores)==F,~agegrp2,Belgium2.uw,svytotal,na.rm = T)[,3]
pos.mix<-tapply(Z==1,Belgium1$agegrp2,sum)
pos.mix.lb<-tapply(lb.Z==1,Belgium1$agegrp2,sum)
pos.mix.ub<-tapply(ub.Z==1,Belgium1$agegrp2,sum)
tot.mix<-tapply(Z!=2,Belgium1$agegrp2,sum)
F.mix<-pos.mix/tot.mix
F.mix.lb<-pos.mix.lb/tot.mix
F.mix.ub<-pos.mix.ub/tot.mix
#F.cut<-pos.cut/tot.cut
F.cut.uw<-pos.cut.uw/tot.cut.uw
#F.GR<-c(0.112,0.163,0.402,0.485,0.653,0.698,.714,0.828,.781,.814,.895,.864,.90,.919,.957,.918,.948,.959,.958,.937)
age_point<-1:19

agelabel<-c("1y","2y","3y","4y","5y","6y","7y","8y","9y","10-14y","15-19y","20-24y",
            "25-29y","30-34y","35-39y","40-44y","45-49y","50-59y","60-69y")

#svg("figs/Belgium/Fig2_be-poly-farr.svg",width=8.,height=5.)
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
#dev.off()

agelabel1<-c(paste(seq(0,55,5),"y",sep=""),">60y")
#svg("figs/Belgium/Fig3_be-poly-farr.svg",width=5.,height=5.)
#tiff("figs/Fig3.tif",width=480,height=480,compression = "lzw")
# Plot prevalence
par(mgp=c(4,1,0), mar=c(5.5,5.5,0,1),cex.lab=1.5,cex.axis=1.5)
plot(age1,prev,main="Belgium",lwd=3,
     xlab="Age (years)",ylab="Seroprevalence (%)",
     ylim=c(0,1),ty="l",col="black",
     axes=F)
lines(age1,lb.prev,lty=2)
lines(age1,ub.prev,lty=2)
polygon(c(age1,rev(age1)),c(ub.prev,rev(lb.prev)),
        col=rgb(0, 0, 0, 0.5), border=NA)
points(age1,prop.mix,cex=0.01*n.mix)
axis(side=1,at=seq(0,60,5),agelabel1,las=3)
axis(side=2,at=seq(0,1,0.1),seq(0,100,10),las=2)  
#dev.off()

#svg("figs/Belgium/Fig4_be-poly-farr.svg",width=5.,height=5.)
# Plot force of infection
par(mgp=c(4,1,0), mar=c(5.5,5.5,0,1),cex.lab=1.5,cex.axis=1.5)
plot(age1[-(Nage-1):-Nage],foi[-(Nage-1):-Nage],main="",lwd=3,
     xlab="Age (years)",ylab="FOI",
     ylim=c(0,1.),xlim=c(0,60),ty="l",col="black",
     axes=F)
lines(age1[(Nage-2):Nage],foi[(Nage-2):Nage],lty=3,lwd=2)
lines(age1,lb.foi,lty=2)
lines(age1,ub.foi,lty=2)
polygon(c(age1,rev(age1)),c(ub.foi,rev(lb.foi)),
        col=rgb(0, 0, 0, 0.5), border=NA)
axis(side=1,at=seq(0,60,5),agelabel1,las=3)
axis(side=2,at=seq(0,1.,0.1),las=2)  
#dev.off()

##### Comparison between models #####
load("b19_prodbeta_sn_res.RData")
dev.np<-b19.prodbeta.sn.orig.new$BUGSoutput$sims.array[,,"deviance"]
load("b19_pcfoi_sir_sn_res.RData")
dev.pc<-b19.pcfoi.sn.orig.new$BUGSoutput$sims.array[,,"deviance"]

summaryDiff(dev.np,dev.pc) # the two models fit similarly

#svg("figs/diffdev_b19.svg",width=5.,height=5.)
par(mgp=c(4.5,1,0), mar=c(5.75,5.75,1.5,1.5),cex.lab=1.5,cex.axis=1.5,las=2)
plot(ecdf(dev.np),lwd=3,xlab="Posterior deviance",main="",axes=F,xlim=c(-9000,-6000))
lines(ecdf(dev.pc),lty=2,col=2,lwd=3)
legend("right",c("PB","SIR"),lwd=3,col=1:2,lty=1:2,bty="n",cex=1.5)
axis(side=1,at=seq(-9000,-6000,500))
axis(side=2,at=seq(0,1.,0.1)) 
#dev.off()

##### Saving results for further analyses #####

### Mixture parameters
DataMix<-data.frame(status=c("susceptible","immune"),
                    mu,lb.mu,ub.mu,
                    sigma,lb.sigma,ub.sigma,
                    alpha,lb.alpha,ub.alpha,
                    f=c(f,NA),lb.f=c(ub.f,NA),ub.f=c(lb.f,NA),
                    country=rep("be",2))
#write.csv(DataMix,file="DataMixtureParams_sir_b19_be.csv",row.names = F)
print(DataMix,digits=2)
#DataMix<-read.csv("tmp/DataMixtureParams_pb_b19_be.csv",header = T)

n.mix<-as.numeric(n.mix)
PC<-data.frame(age1,prop.mix,n.mix,prev,lb.prev,ub.prev,foi,lb.foi,ub.foi)
write.csv(PC,file="DataPrevFOI_sir_b19_be.csv",row.names=F)

##### Plotting non-parametric and PC-FOI models together #####

# Reading data
BP<-read.csv("results/DataPrevFOI_pb_b19_be.csv",header=T)
PC<-read.csv("results/DataPrevFOI_sir_b19_be.csv",header=T)

agelabel1<-c(seq(0,55,5),">60")
#svg("figs/Fig2b.svg",width=5.,height=5.)
#tiff("figs/Fig3.tif",width=480,height=480,compression = "lzw")
# Plot prevalence
par(cex.lab=1.5, cex.main=1.5, cex.sub=1.5,cex.axis=1.5,mgp=c(3.5,1,0),mar=c(5,5.,3,5.))
plot(BP$age1,BP$prev,main="Belgium",lwd=3,
     xlab="Age (years)",ylab="Seroprevalence B19 (%)",
     ylim=c(0,1),xlim=c(0,60),ty="l",col="black",
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

axis(side=1,at=seq(0,60,5),agelabel1,las=3)
axis(side=2,at=seq(0,1,0.1),seq(0,100,10),las=2)
axis(side=4,at=seq(0,1,0.1),las=2)
mtext(text="FOI",side=4,cex=1.5,line=3.5,las=3)
legend("right",lty=1:2,col=1:2,lwd=3,c("PB","SIR"),bty="n",cex = 1.5)
#dev.off()
