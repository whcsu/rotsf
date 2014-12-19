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

set.seed(123)
data(pbc, package = "randomForestSRC")
pbc=na.omit(pbc)
n=dim(pbc)[1]

L=sample(1:n,ceiling(n*0.8))
trset<-pbc[L,]
teset<-pbc[-L,]

source("rotsf-funs.R")

rotsf.fit=rotsf(Surv(days, status)~., data=trset)
rotsf_pre=rotsf.predict(rotsf.fit,teset)

ci=concordance.index(rotsf_pre,teset$days,teset$status)
ci[1]
