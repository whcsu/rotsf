rotsf<-function (formula, data, trlength=500,m=2,control = control, na.action =  na.omit,vari_status=FALSE) 
{
  Call <- match.call()
  
  data=na.omit(data)
  
  mf <- model.frame(formula, data)
  
  x <- model.matrix(attr(mf, "terms"), data = mf)
  
  y <- model.response(mf)
  
   
  if (!inherits(y, "Surv"))
    stop("Response must be a 'survival' object - use the 'Surv()' function")
  
  ny <- ncol(y)
  n <- nrow(y)
  
  status <- y[, ny]
  survtime=y[, 1L]
  

  
  if (any(survtime <= 0)) stop("Observation time must be > 0")
  
  if (all(status == 0)) stop("No deaths in training data set")
  
  
  if (!missing(control)) 
    controls[names(control)] <- control
  

  #names(data)=c("time","status",names(x))
  

  p=dim(mf)[2]
  
  rotms<- vector(mode = "list", length = trlength)
  
  trees <- vector(mode = "list", length = trlength)
  varimp<-NULL
  
  for (i in 1:trlength)
  {
    #create a bagging version for rotation
    
    trainindex=sample(nrow(mf),replace=T)
    
    trsetold=mf[trainindex,]      
    
    # oobset  
    train_posp<-1:nrow(mf) %in% trainindex
    
    oobset=trset[!train_posp,]      
    
    # rotation for gr i bagging version  
    pp=c(1:p)
    
    #excluding reponse index:          
    rii=c(1)
    
    #the rest predictors      
    pp=setdiff(pp,rii)
    #how many groups in the rotation       
    rp=p-length(rii)
    gr=(rp)/m
    
    #initialize varimportance
    varimport=c(rep(0,rp)) 
    
    
    #initialize the rotation matrix 
    rotationm=matrix(rep(0,(rp)^2),nrow=rp,ncol=rp)
    for (j in 1:gr)
    {
      trainindex=sample(nrow(trsetold),replace=T)
      trainb=trsetold[trainindex,]  
      d=sample(pp,m)        
      olddata=trsetold[trainindex,d]
      
      #for pca
      rotj=prcomp(olddata)
      
      
      for (jj in 1:m) {
        #for pca
        rotationm[d-length(rii),(j-1)*m+jj]=rotj$rotation[,jj]            	
        
      } 
      pp=setdiff(pp,d)
    }
    dd=rotationm
    # those features not used in the current pca is set to 0		
    rotationm[pp-length(rii),]=0
    
    # the final the rotation matrix 
    rotms[[i]]=rotationm
    
    
    # the final bagging training set 
    rmatrix=as.matrix(trainb[,-c(rii)]) %*% t(rotationm)      
    rmatrixnew=as.data.frame(rmatrix)      
    rmatrixnew=cbind(trainb[,rii],rmatrixnew)      
    names(rmatrixnew)=names(trainb)	   
    
    trees[[i]]=rpart(rmatrixnew,control = control)
  
if (vari_status)    
{    
    # preparing for predicting oob before permutation      
    rmatrix0=as.matrix(oobset[,-c(rii)]) %*% t(rotationm)
    rmatrixnew0=as.data.frame(rmatrix0)
    rmatrixnew0=cbind(oobset[,rii],rmatrixnew0)
    names(rmatrixnew0)=names(oobset)	       
    predict_oob<-predict(trees[[i]],rmatrixnew0[,-c(rii)])      
    ci_before=concordance.index(predict_oob,oobset$time,oobset$status)
    
    
    # CAL c-index DECREASE
    
    #test first covariate
    X=cbind(sample(oobset[,3]),oobset[,c(4:p)])       
    rmatrix1=as.matrix(X) %*% t(rotationm)
    rmatrixnew0=as.data.frame(rmatrix1)
    rmatrixnew0=cbind(oobset[,rii],rmatrixnew0)
    names(rmatrixnew0)=names(oobset)		      
    predict_oob1<-predict(trees[[i]],rmatrixnew0[,-c(rii)])      
    ci_after=concordance.index(predict_oob1,oobset$time,oobset$status)      
    varimport[1]=unlist(ci_before[1])-unlist(ci_after[1]) 
    
    #test last column
    X=cbind(oobset[,3:(p-1)],sample(oobset[,p]))      
    rmatrix1=as.matrix(X) %*% t(rotationm)
    rmatrixnew0=as.data.frame(rmatrix1)
    rmatrixnew0=cbind(oobset[,rii],rmatrixnew0)
    names(rmatrixnew0)=names(oobset)      
    predict_oob1<-predict(trees[[i]],rmatrixnew0[,-c(rii)])      
    ci_after=concordance.index(predict_oob1,oobset$time,oobset$status)      
    varimport[p-2]=unlist(ci_before[1])-unlist(ci_after[1]) 
    
    
    #test others
    for (jjj in 4:(p-1)) {
      X=cbind(oobset[,3:(jjj-1)],sample(oobset[,jjj]),oobset[,(jjj+1):p])
      
      rmatrix1=as.matrix(X) %*% t(rotationm)
      rmatrixnew0=as.data.frame(rmatrix1)
      rmatrixnew0=cbind(oobset[,rii],rmatrixnew0)
      names(rmatrixnew0)=names(oobset)		        
      predict_oob1<-predict(trees[[i]],rmatrixnew0[,-c(rii)])        
      ci_after=concordance.index(predict_oob1,oobset$time,oobset$status)        
      varimport[jjj-2]=unlist(ci_before[1])-unlist(ci_after[1])        
    }
    if (vari_status)  {  varimp=cbind(varimp,varimport)  }
}
  }
  
  fit=trees
if (vari_status)   
  {
  vari=rowMeans(varimp)
  #make the difference apparent
  vari=vari*trlength
  
  names(vari)=names(mf[,-c(rii)])
}
  class(fit) <- "rotsf"  
if (vari_status)  
  return(list(trees=trees,rotms=rotms,rii=rii,vari=vari))
else
  return(list(trees=trees,rotms=rotms,rii=rii))
  
}

rotsf.predict<-function(rotsffit,newdata,trlength=500){
  
  trees=rotsffit$trees
  rotms=rotsffit$rotms
  rii=rotsffit$rii
  
  if (trlength>length(rotsf.fit$rotms)) 
    stop("Number of Trees for prediction should not be more than Number of Trees Fitted")
  
  # classify the test data   
  testpre<-NULL
  for (i in 1:trlength) {
    #if (oobacc[i]<=avroobacc)
{
  
  # preparing for testing 
  rmatrix2=as.matrix(newdata) %*% t(rotms[[i]])
  rmatrixnew2=as.data.frame(rmatrix2)

  names(rmatrixnew2)=names(newdata)  
  
  predicts<-predict(trees[[i]],rmatrixnew2)
  
  testpre<-cbind(predicts,testpre)
}
  }

ensemble_predictions<-rowMeans(testpre)

return(ensemble_predictions)

}
