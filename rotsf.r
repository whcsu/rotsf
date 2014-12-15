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
p=dim(pbc)[2]

testresult <-data.frame(row.names=1:200)

#AUC
auc_rotsf<-c(rep(0,200))
auc_rsf<-c(rep(0,200))
auc_cox<-c(rep(0,200))
auc_gbm<-c(rep(0,200))

ci_rotsf<-c(rep(0,200))
ci_rsf<-c(rep(0,200))
ci_cox<-c(rep(0,200))
ci_gbm<-c(rep(0,200))

for ( kkk in 1:200) {


L=sample(1:n,ceiling(n*0.8))
trset<-pbc[L,]
teset<-pbc[-L,]

k=2
trlength=1000
rotms<- vector(mode = "list", length = trlength)

trees <- vector(mode = "list", length = trlength)


for (i in 1:trlength)
   {
 #create a bagging version for rotation
  
  trainindex=sample(nrow(trset),replace=T)
  trsetold=trset[trainindex,]  
   # rotation for gr i bagging version  
	pp=c(1:p)
	#excluding reponse index:
	ri=c(1,2)
	rii=c(1,2)
	#the rest predictors 
	pp=setdiff(pp,ri)
	#how many groups in the rotation 
	
	rp=p-length(ri)
	gr=(rp)/k
	#initialize the rotation matrix 
	rotationm=matrix(rep(0,(rp)^2),nrow=rp,ncol=rp)
		for (j in 1:gr)
		{
		  trainindex=sample(nrow(trsetold),replace=T)
		  trainb=trsetold[trainindex,]  
		d=sample(pp,k)
		olddata=trsetold[trainindex,d]
		#for pca
		rotj=prcomp(olddata)
		
		#for kernel pca
	#	rotj=kpca(olddata)
	#	kk=dim(pcv(rotj))[2]
		#kk for kernelpca
		#k for pca
		for (jj in 1:k) {
            #for pca
			rotationm[d-length(rii),(j-1)*k+jj]=rotj$rotation[,jj]          		
		#	for kernelpca
		#	rotationm[d-length(rii),(j-1)*kk+jj]=pcv(rotj)[,jj] 
	     } 
		pp=setdiff(pp,d)
		}
    dd=rotationm
   # those features not used in the current pca is set to 1		
	   rotationm[pp-length(rii),]=0
    
	# the final the rotation matrix 
	rotms[[i]]=rotationm
	
	
	# the final bagging training set 
	 rmatrix=as.matrix(trainb[,-c(rii)]) %*% t(rotationm)
	 rmatrixnew=as.data.frame(rmatrix)
	 rmatrixnew=cbind(trainb[,ri],rmatrixnew)
	 names(rmatrixnew)=names(trainb)	 
	
	trees[[i]]=rpart(Surv(days, status)~., data = rmatrixnew)
  
}

   
   # classify the test data   
testpre<-NULL
for (i in 1:trlength) {
#if (oobacc[i]<=avroobacc)
{

# preparing for testing 
  rmatrix2=as.matrix(teset[,-c(rii)]) %*% t(rotms[[i]])
  rmatrixnew2=as.data.frame(rmatrix2)
  rmatrixnew2=cbind(teset[,ri],rmatrixnew2)
  names(rmatrixnew2)=names(teset)	

	 predicts<-predict(trees[[i]],rmatrixnew2[,-c(rii)])

testpre<-cbind(predicts,testpre)
}
}
ensemble_predictions<-rowMeans(testpre)

#cindex=rcorr.cens(ensemble_predictions,Surv(teset$days,teset$status))
#print(cindex[1])

 
 ci=concordance.index(ensemble_predictions,teset$days,teset$status)
 ci_rotsf[kkk]=ci[1]
 
 
# rfsrc 
rsf=rfsrc(Surv(days, status)~., data = trset) 
rsfpre<-predict(rsf,teset[,-c(rii)])
ci=concordance.index(rsfpre$predicted,teset$days,teset$status)
ci_rsf[kkk]=ci[1]

cox.obj <- coxph(Surv(days, status)~., data = trset)
coxpre<-predict(cox.obj,teset[,-c(rii)])
ci=concordance.index(coxpre,teset$days,teset$status)
 ci_cox[kkk]=ci[1]

 # bag.obj <- bagging(Surv(days, status)~., data = trset)
# bagpre<-predict(bag.obj,newdata=teset[,-c(rii)])
# ci=concordance.index(bagpre,teset$days,teset$status)
# print("bag result:")
 # print(ci[1])
 
 gbm.obj <- gbm(Surv(days, status)~.,distribution="coxph", data = trset,n.trees=1000)
gbmpre<-predict(gbm.obj,teset[,-c(rii)],n.trees=1000)
ci=concordance.index(gbmpre,teset$days,teset$status)
ci_gbm[kkk]=unlist(ci[1])


print(200-kkk)

}


testresult=cbind(testresult,data.frame(auc_rotsf,auc_rsf,auc_cox,auc_gbm,unlist(ci_rotsf),unlist(ci_rsf),unlist(ci_cox),unlist(ci_gbm)))



write.csv(testresult,"survivalresult200.csv")
