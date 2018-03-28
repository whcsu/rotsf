##' @export
rrotsfpca.predict<-function(rrotsfpcafit,newdata,trlength=500){

  trees=rrotsfpcafit$pectrees
  rotms=rrotsfpcafit$rotms
  colindexes=rrotsfpcafit$colindexes


  if (trlength>length(rrotsfpcafit$rotms))
    stop("Number of Trees for prediction should not be more than Number of Trees Fitted")

  # classify the test data
  testpre<-NULL
  for (i in 1:trlength) {
    #if (oobacc[i]<=avroobacc)
{

  # preparing for testing
  testdata=newdata[,colindexes[[i]]]
  rmatrix2=as.matrix(testdata) %*% t(rotms[[i]])
  rmatrixnew2=as.data.frame(rmatrix2)

  colnames(rmatrixnew2)=colnames(testdata)

  predicts<-predict(trees[[i]]$rpart,rmatrixnew2)

  testpre<-cbind(predicts,testpre)
}
  }

ensemble_predictions<-rowMeans(testpre)

return(ensemble_predictions)

}
