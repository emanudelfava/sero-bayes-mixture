#####################################################################################
###################### MIXTURE MODELLING FOR VZV IN TUSCANY 2005 ####################
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

Tuscany<-read.csv("serodata_mev_it_0506.csv",header=T,dec=".",sep=";")
Tuscany<-Tuscany[Tuscany$od>0.001,]
head(Tuscany) # N=927
summary(Tuscany$age)

### Building design weights by age group
Tuscany$agestrata<-cut(Tuscany$age,
                      c(0:10,15,20,25,30,35,40,45,50),
                      right=F,include.lowest=T)
table(Tuscany$agestrata)

Tuscany$agesample<-cut(Tuscany$age,
                      c(0:50),
                      right=F,include.lowest=T)
table(Tuscany$agesample)

Tuscany$id<-1:dim(Tuscany)[1]

### Creating age groups for analysis (aggregating data from 60 years onward)
Tuscany<-Tuscany[order(Tuscany$agesample),]
Tuscany$agep<-rep(c(0:49),table(Tuscany$agesample))
prop.table(table(Tuscany$cutoff)) # F=0.86
Tuscany$Y<-Tuscany$logod

head(Tuscany)

Tuscany1<-subset(Tuscany,age!=0) #I removed the equivocal cases and age 0
N<-dim(Tuscany1)[1] # N=927
#Tuscany1.dw<-svydesign(ids=~id,data=Tuscany1,weights=~dweight) # weighted sample
Tuscany1.uw<-svydesign(ids=~id,data=Tuscany1,weights=~1) # unweighted sample

##### AGE-INDEPENDENT MODEL #####

### Other data
Y<-Tuscany1$Y
Nsub<-length(Y)
K<-4
e0<-rep(1,K)
mea05.data<-list("Y"=Y,"Nsub"=Nsub,"e0"=e0,"K"=K)

#Parameters
mea05.params<-c("xi","sigma","p","Z","ppp")

#Initial values
mea05.inits<-function(){
  list("tau"=c(4.58,9.65,18.11,66.67),"xi0"=c(1.85,2.97,3.55,3.93),
       "S"=sample(1:K,Nsub,replace=TRUE),
       "p"=as.numeric(rdirichlet(1,e0)))}

#Fit JAGS
date()
#cl<-makeCluster(3,type="SOCK")
mv.const.n.4comp.dv<-jags(mea05.data,mea05.inits,mea05.params,n.chains=2,
                           model.file="mev_const_n_dv.txt",n.iter=30000,n.thin = 10)
print(mv.const.n.4comp.dv,digits=2)
#mv.nonparam.sn.4comp.cv.up<-update(mv.nonparam.sn.4comp.cv,n.iter=30000)
#stopCluster(cl)
date()
save(mv.const.n.4comp.dv,file="mev_const_n_4comp_dv.RData")

fit<-mv.const.n.4comp.dv
traceplot(fit,varname="xi")
traceplot(fit,varname="sigma")
traceplot(fit,varname="p")

mu<-fit$BUGSoutput$mean$xi
sigma<-fit$BUGSoutput$mean$sigma
Z<-fit$BUGSoutput$median$Z
p<-fit$BUGSoutput$mean$p

x<-seq(min(Y),max(Y),length.out = 100)
mix.n<-p[1]*dnorm(x,mu[1],sigma[1])+
  p[2]*dnorm(x,mu[2],sigma[2])+
  p[3]*dnorm(x,mu[3],sigma[3])+
  p[4]*dnorm(x,mu[4],sigma[4])

#svg("figs/Histmix_4comp_mev.svg",width=5.,height=5.)
par(cex.lab=1.5, cex.main=1.5, cex.sub=1.5,cex.axis=1.5,mgp=c(3.5,1,0),mar=c(5,5.5,3,2),las=1)
hist(Tuscany1$Y,breaks=100,xlab=expression("MeV ("*log[10]*"[OD+1])"),
     freq=F,main="Tuscany (Italy)",ylim=c(0,1))
abline(v=2.2,lty=2,lwd=2)
abline(v=2.5,lty=2,lwd=2)
lines(x,p[1]*dnorm(x,mu[1],sigma[1]),col=1,lty=3,lwd=2)
lines(x,p[2]*dnorm(x,mu[2],sigma[2]),col=1,lty=3,lwd=2)
lines(x,p[3]*dnorm(x,mu[3],sigma[3]),col=1,lty=3,lwd=2)
lines(x,p[4]*dnorm(x,mu[4],sigma[4]),col=1,lty=3,lwd=2)
lines(x,mix.n,col=1,lwd=3)
#dev.off()

##### NONPARAMETRIC MODEL #####

### Other data
Y<-Tuscany1$Y
Nsub<-length(Y)
age.g<-Tuscany1$agep
age1<-unique(age.g)
Nage<-length(age1)
K<-4
e0<-rep(1,K)
#mean.y<-mean(Y)
#var.y<-var(Y)
#G0<-1/(0.5*var.y)
mea05.data<-list("Y"=Y,"Nsub"=Nsub,"Nage"=Nage,"age.g"=age.g,"e0"=e0,"K"=K)

hist(Y,breaks=100,xlab=expression("IgG antibodies to MeV ("*log[10]*"[OD+1])"),
     freq=F,cex.lab=1.5, cex.main=2, cex.sub=2, mgp=c(2.5,1,0),
     cex.axis=1.75,col=2,main="")

#Parameters
mea05.params<-c("xi","sigma","prevs","Z","ppp")

#Initial values
mea05.inits<-function(){
  list("tau"=c(4.58,9.65,18.11,66.67),"xi0"=c(1.85,2.97,3.55,3.93),
       "S"=sample(1:K,Nsub,replace=TRUE),
       "p"=rdirichlet(Nage,e0))}
#mea05.inits<-function(){
#  list("tau"=runif(K,0,100),"xi0"=runif(K,min(Y),max(Y)),
#       "S"=sample(1:K,Nsub,replace=TRUE),
#       "p"=rdirichlet(Nage,e0))}

#Fit JAGS
date()
#cl<-makeCluster(3,type="SOCK")
mv.nonparam.n.4comp.dv<-jags(mea05.data,mea05.inits,mea05.params,n.chains=2,
                  model.file="mev_nonparam_n_dv.txt",n.iter=30000,n.thin = 10)
#mv.nonparam.n.3comp.dv.up<-update(mv.nonparam.n.3comp.dv,n.iter=10000)
#stopCluster(cl)
date()
print(mv.nonparam.n.4comp.dv,digits=2)
save(mv.nonparam.n.4comp.dv,file="mev_nonparam_n_4comp_dv_res.RData")
load("mev_nonparam_n_4comp_dv_res.RData")
recompile(mv.nonparam.n.4comp.dv)
dic.samples(mv.nonparam.n.4comp.dv$model,type="popt",n.iter=1000) 
# Inversa gamma for variance, PED=8552 (Mean dev = 272.9)

##### ANALYSIS OF RESULTS OF JAGS MODEL #####
fit<-mv.nonparam.n.4comp.dv

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

xi3<-mcmc(fit$BUGSoutput$sims.matrix[,"xi[3]"])
geweke.diag(xi3)
heidel.diag(xi3)
xi3<-mcmc.list(mcmc(fit$BUGSoutput$sims.matrix[1:1000,"xi[3]"]),
               mcmc(fit$BUGSoutput$sims.matrix[1001:2000,"xi[3]"]))
gelman.diag(xi3)
plot(xi3)

xi4<-mcmc(fit$BUGSoutput$sims.matrix[,"xi[4]"])
geweke.diag(xi4)
heidel.diag(xi4)
xi4<-mcmc.list(mcmc(fit$BUGSoutput$sims.matrix[1:1000,"xi[4]"]),
               mcmc(fit$BUGSoutput$sims.matrix[1001:2000,"xi[4]"]))
gelman.diag(xi4)
plot(xi4)

traceplot(fit,varname="sigma")
sigma1<-mcmc(fit$BUGSoutput$sims.matrix[,"sigma[1]"])
geweke.diag(sigma1)
heidel.diag(sigma1)
sigma1<-mcmc.list(mcmc(fit$BUGSoutput$sims.matrix[1:1000,"sigma[1]"]),
                  mcmc(fit$BUGSoutput$sims.matrix[1001:2000,"sigma[1]"]))
gelman.diag(sigma1)
plot(sigma1)

sigma2<-mcmc(fit$BUGSoutput$sims.matrix[,"sigma[2]"])
geweke.diag(sigma2)
heidel.diag(sigma2)
sigma2<-mcmc.list(mcmc(fit$BUGSoutput$sims.matrix[1:1000,"sigma[2]"]),
                  mcmc(fit$BUGSoutput$sims.matrix[1001:2000,"sigma[2]"]))
gelman.diag(sigma2)
plot(sigma2)

sigma3<-mcmc(fit$BUGSoutput$sims.matrix[,"sigma[3]"])
geweke.diag(sigma3)
heidel.diag(sigma3)
sigma3<-mcmc.list(mcmc(fit$BUGSoutput$sims.matrix[1:1000,"sigma[3]"]),
               mcmc(fit$BUGSoutput$sims.matrix[1001:2000,"sigma[3]"]))
gelman.diag(sigma3)
plot(sigma3)

sigma4<-mcmc(fit$BUGSoutput$sims.matrix[,"sigma[4]"])
geweke.diag(sigma4)
heidel.diag(sigma4)
sigma4<-mcmc.list(mcmc(fit$BUGSoutput$sims.matrix[1:1000,"sigma[4]"]),
                  mcmc(fit$BUGSoutput$sims.matrix[1001:2000,"sigma[4]"]))
gelman.diag(sigma4)
plot(sigma4)

#Saving main parameters from piecewise-constant FOI model
Z<-fit$BUGSoutput$median$Z
lb.Z<-round(fit$BUGSoutput$summary[1:927,"2.5%"])
ub.Z<-round(fit$BUGSoutput$summary[1:927,"97.5%"])
prev<-fit$BUGSoutput$mean$prevs
lb.prev<-fit$BUGSoutput$summary[930:978,"2.5%"]
ub.prev <-fit$BUGSoutput$summary[930:978,"97.5%"]
sigma<-fit$BUGSoutput$mean$sigma
lb.sigma<-fit$BUGSoutput$summary[979:982,"2.5%"]
ub.sigma<-fit$BUGSoutput$summary[979:982,"97.5%"]
mu<-fit$BUGSoutput$mean$xi
lb.mu<-fit$BUGSoutput$summary[983:986,"2.5%"]
ub.mu<-fit$BUGSoutput$summary[983:986,"97.5%"]

Tuscany1$mix.res<-Z
head(Tuscany1)

#Tuscany2.dw<-svydesign(ids=~id,data=Tuscany1,weights=~dweight) # weighted sample
Tuscany2.uw<-svydesign(ids=~id,data=Tuscany1,weights=~1) # unweighted sample
#prop.table(svytable(~mix.res+sex,subset(Tuscany2.dw,age>=20 & age<46)),1)
prop.table(table(Tuscany1[Tuscany1$age>=20 & Tuscany1$age<46,]$mix.res,
                 Tuscany1[Tuscany1$age>=20 & Tuscany1$age<46,]$gender),1)
prop.table(table(Tuscany1[Tuscany1$age>=60,]$mix.res))
# 51.7% of susceptible people in the childbearing period are women.

#Create a table with n° trials and n° successes
prop.table(table(Z)) # POS=76.2% (CI: 25.7% - 87.1%)
prop.table(table(lb.Z))
prop.table(table(ub.Z))

#Create a table with n° trials and n° successes, by age
n.mix<-table(Tuscany1$agep) 
prop.mix<-prop.table(table(Z,Tuscany1$agep),2)[2,] # Seropositive proportions
prop.mix.lb<-prop.table(table(lb.Z,Tuscany1$agep),2)[2,]
prop.mix.ub<-prop.table(table(ub.Z,Tuscany1$agep),2)[2,]

#DataGigi<-data.frame(age1,prop=prop.mix,lb=prop.mix.lb,ub=prop.mix.ub)
#write.csv(DataGigi,file="tmp/DataSeropositivePropFromMixture.csv",row.names = F)

plot(age1,prop.mix,ty="b",ylim=c(0,1))
lines(age1,prop.mix.lb,lty=2)
lines(age1,prop.mix.ub,lty=2)
polygon(c(age1,rev(age1)),c(prop.mix.ub,rev(prop.mix.lb)),
        col=rgb(0, 0, 0, 0.5), border=NA)

Tuscany1$agegrp2<-cut(Tuscany1$age,
                     c(1:9,10,15,20,25,30,35,40,45,50),
                     right=F,include.lowest=T)
#Tuscany1<-Tuscany1[order(Tuscany1$agegrp2),]
#Tuscany2.dw<-svydesign(ids=~id,data=Tuscany1,weights=~dweight) # weighted sample
Tuscany2.uw<-svydesign(ids=~id,data=Tuscany1,weights=~1) # weighted sample

#pos.cut<-svyby(~mres==1,~agegrp2,Tuscany2.dw,svytotal)[,3]
#tot.cut<-svyby(~mres!=2,~agegrp2,Tuscany2.dw,svytotal)[,3]
pos.cut.uw<-svyby(~cutoff==1,~agegrp2,Tuscany2.uw,svytotal,na.rm = T)[,3]
tot.cut.uw<-svyby(~cutoff!=2,~agegrp2,Tuscany2.uw,svytotal,na.rm = T)[,3]
pos.mix<-tapply(Z==1,Tuscany1$agegrp2,sum)
pos.mix.lb<-tapply(lb.Z>0,Tuscany1$agegrp2,sum)
pos.mix.ub<-tapply(ub.Z>0,Tuscany1$agegrp2,sum)
tot.mix<-table(Tuscany1$agegrp2)
F.mix<-pos.mix/tot.mix
F.mix.lb<-pos.mix.lb/tot.mix
F.mix.ub<-pos.mix.ub/tot.mix
#F.cut<-pos.cut/tot.cut
F.cut.uw<-pos.cut.uw/tot.cut.uw
#F.GR<-c(0.112,0.163,0.402,0.485,0.653,0.698,.714,0.828,.781,.814,.895,.864,.90,.919,.957,.918,.948,.959,.958,.937)
age_point<-1:17

agelabel<-c("1y","2y","3y","4y","5y","6y","7y","8y","9y","10-14y","15-19y","20-24y",
            "25-29y","30-34y","35-39y","40-44y","45-49y")

#svg("figs/Fig2_tu05-3comp_dv.svg",width=8.,height=5.)
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

agelabel1<-seq(0,50,5)
#svg("figs/prev_nonparam_m4.svg",width=5.,height=5.)
#tiff("figs/Fig3.tif",width=480,height=480,compression = "lzw")
# Plot prevalence
par(cex.lab=1.5, cex.main=1.5, cex.sub=1.5,cex.axis=1.5,mgp=c(3.5,1,0),mar=c(5,5.5,3,2))
plot(age1,prev,main="Tuscany (Italy)",lwd=3,
     xlab="Age (years)",ylab="Seroprevalence MeV (%)",
     ylim=c(0,1),ty="l",col="black",
     axes=F)
lines(age1,lb.prev,lty=2)
lines(age1,ub.prev,lty=2)
polygon(c(age1,rev(age1)),c(ub.prev,rev(lb.prev)),
        col=rgb(0, 0, 0, 0.5), border=NA)
points(age1,prop.mix,cex=0.075*n.mix)
axis(side=1,at=seq(0,50,5),agelabel1,las=3)
axis(side=2,at=seq(0,1,0.1),seq(0,100,10),las=2)  
#dev.off()

### Mixture parameters
DataMix<-data.frame(status=c("c1","c2","c3","c4"),
                    mu,lb.mu,ub.mu,
                    sigma,lb.sigma,ub.sigma,
                    country=rep("tuscany05",4))
print(DataMix,digits=3)
#write.csv(DataMix,file="results/DataMixtureParams_prob_m4_tu05.csv",row.names = F)

### Seroprevalence and FOI
DataRes<-data.frame(age1,prop.mix,n.mix=as.numeric(n.mix),
                    prev,lb.prev,ub.prev)
#write.csv(DataRes,file="results/DataPrevFOI_prob_m4_tu05.csv",row.names = F)
