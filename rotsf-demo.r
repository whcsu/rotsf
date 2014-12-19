rm(list=ls())
library(rpart)
library(MASS)
library(gbm)
library(kernlab)
library(survival)
library(Hmisc)
library(survcomp)
library(survivalROC)
library(randomForestSRC)

source("rotsf-funs.R")

set.seed(123)

#PBC DATA
data(pbc, package = "randomForestSRC")
pbc=na.omit(pbc)

n=dim(pbc)[1]

L=sample(1:n,ceiling(n*0.8))
trset<-pbc[L,]
teset<-pbc[-L,]



rotsf.fit=rotsf(Surv(days, status)~., data=trset)
rotsf_pre=rotsf.predict(rotsf.fit,teset)

ci=concordance.index(rotsf_pre,teset$days,teset$status)
ci[1]

#CANCER DATA
data(cancer, package = "survival")
cancer[,3]=cancer[,3]-1
cancer=cancer[,-c(1)]
cancer=na.omit(cancer)

n2=dim(cancer)[1]

L2=sample(1:n2,ceiling(n2*0.8))
trset2<-cancer[L2,]
teset2<-cancer[-L2,]



rotsf.fit2=rotsf(Surv(time, status)~., data=trset2)
rotsf_pre2=rotsf.predict(rotsf.fit2,teset2)

ci2=concordance.index(rotsf_pre2,teset2$time,teset2$status)
ci2[1]

#STAGEC DATA
data(stagec, package = "rpart")
stagec[,8]=as.numeric(stagec[,8])
stagec=na.omit(stagec)
n3=dim(stagec)[1]


L3=sample(1:n3,ceiling(n3*0.8))
trset3<-stagec[L3,]
teset3<-stagec[-L3,]



rotsf.fit3=rotsf(Surv(pgtime, pgstat)~., data=trset3)
rotsf_pre3=rotsf.predict(rotsf.fit3,teset3)

ci3=concordance.index(rotsf_pre3,teset3$pgtime,teset3$pgstat)
ci3[1]
