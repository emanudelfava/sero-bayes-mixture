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

# Reading measles data from Tuscany (2005)
ItMeV<-read.csv("data/serodata_mev_it_0506.csv",header=T,dec=".",sep=",")
ItMeV<-ItMeV[ItMeV$od>0.001,]
head(ItMeV) # N=927
summary(ItMeV$age)

# Reading VZV data from Belgium (2001-2003)
BeVZV<-read.table("data/serodata_vzv_b19_be.csv",header=T,sep=",",na.strings = c("",NA))
BeVZV<-BeVZV[is.na(BeVZV$vzvmiuml)==F & BeVZV$vzvmiuml<10000,c(1:3,6:9)]
head(BeVZV) # N=2759
range(BeVZV$age)

# Reading B19 data from Belgium (2001-2003)
BeB19<-read.table("data/serodata_vzv_b19_be.csv",header=T,sep=",",na.strings = c("",NA))
BeB19<-BeB19[is.na(BeB19$parvouml)==F,c(1:5,8:9)]
head(BeB19) # N=3098
range(BeB19$age)

# Reading VZV data from England and Wales (2001-2003)
EWVZV<-read.table("data/serodata_vzv_b19_ew.csv",header=T,sep=",",na.strings = c("",NA))
EWVZV<-EWVZV[is.na(EWVZV$vzvmiuml)==F,c(1:3,6:9)]
head(EWVZV) # N=2091
range(EWVZV$age)

##### PANEL 3X2 WITH HISTOGRAMS AND SCATTER PLOTS #####

### VZV Belgium
# Histogram with fixed cut-off point
#svg("figs/Fig0a.svg",width=5.,height=5.)
par(mfrow=c(1,1),mar=c(5,5.5,3,2),mgp=c(3.5,1,0),
    cex.lab=1.5,cex.main=1.5,cex.sub=1.5,cex.axis=1.5,las=1)
hist(log10(BeVZV$vzvmiuml+1),breaks=100,freq=F,xlim=c(0,4),ylim=c(0,1),
     main="",xlab=expression("VZV ["*log[10]*"(OD+1)]"),axes=F)
abline(v=log10(50+1),lty=2,lwd=2)
abline(v=log10(100+1),lty=2,lwd=2)
axis(1,at=seq(0,4,0.5))
axis(2,at=seq(0,1.,0.25))
#dev.off()

#svg("figs/Fig0b.svg",width=5.,height=5.)
# Scatter plot by age with fixed cut-off point
par(mfrow=c(1,1),mar=c(5,5.5,3,2),mgp=c(3.5,1,0),
    cex.lab=1.5,cex.main=1.5,cex.sub=1.5,cex.axis=1.5,las=1)
plot(BeVZV$age,log10(BeVZV$vzvmiuml+1),xlab="Age (in years)",col="darkgray",
     main="",ylab=expression("VZV ["*log[10]*"(OD+1)]"),ylim=c(0,4),axes=F)
abline(h=log10(50+1),lty=2,lwd=2)
abline(h=log10(100+1),lty=2,lwd=2)
axis(1,at=seq(0,80,10))
axis(2,at=seq(0,4,0.5))
#dev.off()

# Proportion seropositive by age according to current status data
BeVZV1<-BeVZV[BeVZV$vzvmiuml<=50 | BeVZV$vzvmiuml>=100,]
prop<-prop.table(table(BeVZV1$vzvres,BeVZV1$age),2)[3,] # age-specific prevalence from cut-off point
#svg("figs/Fig0c.svg",width=5.,height=5.)
par(mfrow=c(1,1),mar=c(5,5.5,3,2),mgp=c(3.5,1,0),
    cex.lab=1.5,cex.main=1.5,cex.sub=1.5,cex.axis=1.5,las=1)
plot(0:79,prop*100,xlab="Age (in years)",ty="p",
     main="",ylab="VZV Seroprevalence (%)",
     cex=0.025*table(BeVZV1$age),ylim=c(0,100),axes=F)
axis(1,at=seq(0,80,10))
axis(2,at=seq(0,100,10))
#dev.off()

### B19
#svg("figs/Fig1a.svg",width=5.,height=5.)
# Histogram with fixed cut-off point
par(mfrow=c(1,1),mar=c(5,5.5,3,2),mgp=c(3.5,1,0),
    cex.lab=1.5,cex.main=1.5,cex.sub=1.5,cex.axis=1.5,las=1)
hist(log10(BeB19$parvouml+1),breaks=100,freq=F,ylim=c(0,1.5),xlim=c(0.25,2.75),
        main="Belgium",xlab=expression("B19 ["*log[10]*"(OD+1)]"),axes=F)
abline(v=log10(20+1),lty=2,lwd=2)
abline(v=log10(24+1),lty=2,lwd=2)
axis(1,at=seq(0,3,0.5))
axis(2,at=seq(0,1.5,0.3))
#dev.off()

#svg("figs/Fig1b.svg",width=5.,height=5.)
# Scatter plot by age with fixed cut-off point
par(mfrow=c(1,1),mar=c(5,5.5,3,2),mgp=c(3.5,1,0),
    cex.lab=1.5,cex.main=1.5,cex.sub=1.5,cex.axis=1.5,las=1)
plot(BeB19$age,log10(BeB19$parvouml+1),xlab="Age (in years)",col="darkgray",
     main="Belgium",ylab=expression("B19 ["*log[10]*"(OD+1)]"),axes=F)
abline(h=log10(20+1),lty=2,lwd=2)
abline(h=log10(24+1),lty=2,lwd=2)
axis(1,at=seq(0,80,10))
axis(2,at=seq(0,3,0.5))
#dev.off()

### VZV
#svg("figs/Fig1c.svg",width=5.,height=5.)
# Histogram with fixed cut-off point
par(mfrow=c(1,1),mar=c(5,5.5,3,2),mgp=c(3.5,1,0),
    cex.lab=1.5,cex.main=1.5,cex.sub=1.5,cex.axis=1.5,las=1)
hist(log10(EWVZV$vzvmiuml+1),breaks=100,freq=F,ylim=c(0,1.5),xlim=c(0.25,2.75),
     main="England and Wales",xlab=expression("VZV ["*log[10]*"(OD+1)]"),axes=F)
abline(v=log10(15+1),lty=2,lwd=2)
abline(v=log10(20+1),lty=2,lwd=2)
axis(1,at=seq(0,3,0.5))
axis(2,at=seq(0,1.5,0.3))
#dev.off()

#svg("figs/Fig1d.svg",width=5.,height=5.)
# Scatter plot by age with fixed cut-off point
par(mfrow=c(1,1),mar=c(5,5.5,3,2),mgp=c(3.5,1,0),
    cex.lab=1.5,cex.main=1.5,cex.sub=1.5,cex.axis=1.5,las=1)
plot(EWVZV$age,log10(EWVZV$vzvmiuml+1),xlab="Age (in years)",col="darkgray",xlim=c(0,20),
     main="England and Wales",ylab=expression("VZV ["*log[10]*"(OD+1)]"),axes=F)
abline(h=log10(15+1),lty=2,lwd=2)
abline(h=log10(20+1),lty=2,lwd=2)
axis(1,at=seq(0,20,5))
axis(2,at=seq(0,3,0.5))
#dev.off()

### MeV
#svg("figs/Fig1e.svg",width=5.,height=5.)
# Histogram with fixed cut-off point
par(mfrow=c(1,1),mar=c(5,5.5,3,2),mgp=c(3.5,1,0),
    cex.lab=1.5,cex.main=1.5,cex.sub=1.5,cex.axis=1.5,las=1)
hist(ItMeV$logod,breaks=100,freq=F,xlim=c(0.5,4.5),ylim=c(0,1.),
     main="Tuscany (Italy)",xlab=expression("MeV ["*log[10]*"(OD+1)]"),axes=F)
abline(v=2.2,lty=2,lwd=2)
abline(v=2.5,lty=2,lwd=2)
axis(1,at=seq(0.5,4.5,0.5))
axis(2,at=seq(0,1.,0.25))
#dev.off()

#svg("figs/Fig1f.svg",width=5.,height=5.)
# Scatter plot by age with fixed cut-off point
par(mfrow=c(1,1),mar=c(5,5.5,3,2),mgp=c(3.5,1,0),
    cex.lab=1.5,cex.main=1.5,cex.sub=1.5,cex.axis=1.5,las=1)
plot(ItMeV$age,ItMeV$logod,xlab="Age (in years)",col="darkgray",xlim=c(0,50),
     main="Tuscany (Italy)",ylab=expression("MeV ["*log[10]*"(OD+1)]"),axes=F)
abline(h=2.2,lty=2,lwd=2)
abline(h=2.5,lty=2,lwd=2)
axis(1,at=seq(0,50,10))
axis(2,at=seq(0.5,4.5,0.5))
#dev.off()



