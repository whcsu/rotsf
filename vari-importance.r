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

testresult <-data.frame(row.names=1:30)



L=sample(1:n,ceiling(n*0.8))
trset<-pbc[L,]
teset<-pbc[-L,]

k=3
trlength=1000
rotms<- vector(mode = "list", length = trlength)

trees <- vector(mode = "list", length = trlength)
varimp<-NULL

for (i in 1:trlength)
   {
 #create a bagging version for rotation
  trainindex=sample(nrow(trset),replace=T)
  trsetold=trset[trainindex,]  
    
  
    
    # oobset  
    train_posp<-1:nrow(trset) %in% trainindex
    
    oobset=trset[!train_posp,]
	
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
  
	#initialize varimportance
	varimport=c(rep(0,rp))  
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
	   rotationm[pp-length(rii),]=1
    
	# the final the rotation matrix 
	rotms[[i]]=rotationm
	
	
	# the final bagging training set 
	 rmatrix=as.matrix(trainb[,-c(rii)]) %*% t(rotationm)
	 rmatrixnew=as.data.frame(rmatrix)
	 rmatrixnew=cbind(trainb[,ri],rmatrixnew)
	 names(rmatrixnew)=names(trainb)	 
	
	trees[[i]]=rpart(Surv(days, status)~., data = rmatrixnew)
  
	# preparing for predicting oob before permutation
  
	rmatrix0=as.matrix(oobset[,-c(rii)]) %*% t(rotationm)
	rmatrixnew0=as.data.frame(rmatrix0)
	rmatrixnew0=cbind(oobset[,ri],rmatrixnew0)
	names(rmatrixnew0)=names(oobset)	 
    
	predict_oob<-predict(trees[[i]],rmatrixnew0[,-c(rii)])
	
	ci_before=concordance.index(predict_oob,oobset$time,oobset$status)
  
  
	# CAL c-index DECREASE
	
	#test first covariate
	X=cbind(sample(oobset[,3]),oobset[,c(4:p)])
	

	rmatrix1=as.matrix(X) %*% t(rotationm)
	rmatrixnew0=as.data.frame(rmatrix1)
	rmatrixnew0=cbind(oobset[,ri],rmatrixnew0)
	names(rmatrixnew0)=names(oobset)		
	
  predict_oob1<-predict(trees[[i]],rmatrixnew0[,-c(rii)])
	
	ci_after=concordance.index(predict_oob1,oobset$time,oobset$status)
	
	
	varimport[1]=unlist(ci_before[1])-unlist(ci_after[1]) 
	
	#test last column
	X=cbind(oobset[,3:(p-1)],sample(oobset[,p]))
  
	rmatrix1=as.matrix(X) %*% t(rotationm)
	rmatrixnew0=as.data.frame(rmatrix1)
	rmatrixnew0=cbind(oobset[,ri],rmatrixnew0)
	names(rmatrixnew0)=names(oobset)		
	
	predict_oob1<-predict(trees[[i]],rmatrixnew0[,-c(rii)])
	
	ci_after=concordance.index(predict_oob1,oobset$time,oobset$status)
	
	
	varimport[p-2]=unlist(ci_before[1])-unlist(ci_after[1]) 

	
	#test others
	for (jjj in 4:(p-1)) {
	  X=cbind(oobset[,3:(jjj-1)],sample(oobset[,jjj]),oobset[,(jjj+1):p])
	 
	  rmatrix1=as.matrix(X) %*% t(rotationm)
	  rmatrixnew0=as.data.frame(rmatrix1)
	  rmatrixnew0=cbind(oobset[,ri],rmatrixnew0)
	  names(rmatrixnew0)=names(oobset)		
	  
	  predict_oob1<-predict(trees[[i]],rmatrixnew0[,-c(rii)])
	  
	  ci_after=concordance.index(predict_oob1,oobset$time,oobset$status)
	  
	  
	  varimport[jjj-2]=unlist(ci_before[1])-unlist(ci_after[1]) 
	}
 print(i)
	varimp=cbind(varimp,varimport)	 
  
  
}

write.csv(varimp,"var.csv" )
dd=rowMeans(varimp)
names(dd)=names(pbc[3:19])
dd=dd*100
write.csv(dd,"var2.csv")
dd=read.csv("var2.csv")
dd=dd[order(-dd[,2]),]
y=dd
#barplot(y[1:20,2],beside = TRUE,xlim=c(0,200),xlab="Importance",main="Variable Importance by AUC Decrease",names.arg=y[1:20,1],horiz=TRUE,cex.names = 0.8,las=2,)
par(mar=c(5,17,2,2))
barplot(y[1:10,2],beside = TRUE,xlim=c(0,10),xlab="Importance",names.arg=y[1:10,1],horiz=TRUE,cex.names = 0.8,las=2,)

