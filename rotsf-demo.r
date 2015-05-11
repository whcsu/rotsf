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

Rn=1000
testresult <-data.frame(row.names=1:Rn)

#AUC
ci_rotsf<-c(rep(0,Rn))

ci_rsf<-c(rep(0,Rn))
ci_cox<-c(rep(0,Rn))
ci_gbm<-c(rep(0,Rn))

for (i in 1:Rn)
   {
L=sample(1:n,ceiling(n*0.8))
trset<-pbc[L,]
teset<-pbc[-L,]
rii=c(1,2)
print(i)
# rotsf

rotsf.fit=rotsf(Surv(days, status)~., data=trset,trlength=1000)
rotsf_pre=rotsf.predict(rotsf.fit,teset[,-c(rii)],trlength=1000)
ci1=concordance.index(rotsf_pre,teset$days,teset$status)
ci_rotsf[i]=unlist(ci1[1])



# rfsrc 
rsf=rfsrc(Surv(days, status)~., data = trset) 
rsfpre<-predict(rsf,teset[,-c(rii)])
ci2=concordance.index(rsfpre$predicted,teset$days,teset$status)
ci_rsf[i]=unlist(ci2[1])

#cox
cox.obj <- coxph(Surv(days, status)~., data = trset)
coxpre<-predict(cox.obj,teset[,-c(rii)])
ci3=concordance.index(coxpre,teset$days,teset$status)
ci_cox[i]=unlist(ci3[1])

#gbm
gbm.obj <- gbm(Surv(days, status)~.,distribution="coxph", data = trset,n.trees=1000)
gbmpre<-predict(gbm.obj,teset[,-c(rii)],n.trees=1000)
ci4=concordance.index(gbmpre,teset$days,teset$status)
ci_gbm[i]=unlist(ci4[1])

}

testresult=cbind(testresult,data.frame(unlist(ci_rotsf),unlist(ci_rsf),unlist(ci_cox),unlist(ci_gbm)))

write.csv(testresult,"surv-pbc.csv")
