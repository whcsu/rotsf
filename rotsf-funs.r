rotsf<-function (formula, data, trlength=500,k=2,control = control, na.action =  na.omit) 
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
  
  if (all(status == 0)) stop("No deaths in data set")
  

  if (!missing(control)) 
    controls[names(control)] <- control
  
  p=dim(data)[2]
  
    rotms<- vector(mode = "list", length = trlength)
    
    trees <- vector(mode = "list", length = trlength)
    
    
    for (i in 1:trlength)
    {
      #create a bagging version for rotation
      
      trainindex=sample(nrow(data),replace=T)
      
      trsetold=data[trainindex,]  
      
      # rotation for gr i bagging version  
      pp=c(1:p)
      
      #excluding reponse index:
          
      rii=c(1,2)
      
      #the rest predictors
      
      pp=setdiff(pp,rii)
      #how many groups in the rotation 
      
      rp=p-length(rii)
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
        
  
        for (jj in 1:k) {
          #for pca
          rotationm[d-length(rii),(j-1)*k+jj]=rotj$rotation[,jj]          		
     
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
      
      trees[[i]]=rpart(formula, data = rmatrixnew,control = control)
      
    }
   fit=trees
  
  class(fit) <- "rotsf"
  
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
  rmatrix2=as.matrix(newdata[,-c(rii)]) %*% t(rotms[[i]])
  rmatrixnew2=as.data.frame(rmatrix2)
  rmatrixnew2=cbind(newdata[,rii],rmatrixnew2)
  names(rmatrixnew2)=names(newdata)  
  
  predicts<-predict(trees[[i]],rmatrixnew2[,-c(rii)])
  
  testpre<-cbind(predicts,testpre)
}
}

ensemble_predictions<-rowMeans(testpre)

return(ensemble_predictions)

}
