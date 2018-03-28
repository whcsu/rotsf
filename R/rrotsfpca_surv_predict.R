##' @export
rrotsfpca_surv.predict<-function(rrotsfpcafit,newdata,uniquetimes,trlength=rrotsfpcafit$trlength){

  trees=rrotsfpcafit$pectrees
  rotms=rrotsfpcafit$rotms

  colindexes=rrotsfpcafit$colindexes


  if (trlength>length(rrotsfpcafit$rotms))
    stop("Number of Trees for prediction should not be more than Number of Trees Fitted")

 
  # Preparing testpre dataframe 
  testpre <- matrix(0,nrow = dim(newdata)[1],  ncol= length(uniquetimes))
  colnames(testpre)=paste0(uniquetimes)
  for (i in 1:trlength) {
    #if (oobacc[i]<=avroobacc)
{

  # preparing for testing
  testdata=newdata[,colindexes[[i]]]
  rmatrix2=as.matrix(testdata) %*% t(rotms[[i]])
  rmatrixnew2=as.data.frame(rmatrix2)

  colnames(rmatrixnew2)=colnames(testdata)

      predicts=predictSurvProb(trees[[i]],rmatrixnew2,uniquetimes)
      #convert all NA into zero
      predicts[is.na(predicts)]=0
      colnames(predicts)=paste0(uniquetimes)
      #print((predicts[1,100]))
      testpre<-testpre+predicts
}
  }

 return(testpre/trlength)

}
