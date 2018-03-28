##' Prediction with rrotsfspls and return a list of with mean hazard function 
##' @title rotsf predict.rrotsfspls
##' @param object  An object that inherits from class rrotsfspls.
##' @param testx  A data frame in which to look for variables with which to predict. 
##' @param trlength  Number of based models used for prediction, shouble be less than and equal to the number for training.
##' @param ... Additional arguments.
##' @return produces a vector of predictions or a matrix of predictions
##' @seealso \code{pec} R package
##' @author Hong Wang and Lifeng Zhou
##' @references
##' \itemize{
##'   \item Zhou L, Xu Q, Wang H. (2015) Rotation survival forest for right censored data. PeerJ 3:e1009 https://doi.org/10.7717/peerj.1009.
##'   \item Zhou, L., Wang, H., & Xu, Q. (2016). Random rotation survival forest for high dimensional censored data. SpringerPlus, 5(1), 1425.
##'  }
##' @examples
##' set.seed(123)
##' require(rotsf)
##' require(survival)
##' ## Survival Forest with PLS  with default settings
##' #Lung DATA
##' data(lung)
##' lung=na.omit(lung)
##' lung[,3]=lung[,3]-1
##' n=dim(lung)[1]
##' L=sample(1:n,ceiling(n*0.5))
##' trset<-lung[L,]
##' teset<-lung[-L,]
##' rii=c(2,3)
##' plssurvmodel=rrotsfspls(x=trset[,-rii],y=Surv(trset[,rii[1]], trset[,rii[2]]))
##' mypre=predict(plssurvmodel, teset[,-rii])
##' @export
predict.rrotsfspls<-function(object,testx,trlength=object$trlength,...){

  trees=object$pectrees
  rotms=object$rotms
  colindexes=object$colindexes
  nplscomp=object$nplscomp

  if (trlength>object$trlength)
    stop("Number of Trees for prediction should not be more than Number of Trees Fitted")

  # classify the test data
  testpre<-NULL
  for (i in 1:trlength) {
    #if (oobacc[i]<=avroobacc)
{

  # preparing for testing
  
  testdata=testx[,colindexes[[i]]]
  rmatrix2=as.matrix(testdata) %*% (rotms[[i]])
  rmatrixnew2=as.data.frame(rmatrix2)

  colnames(rmatrixnew2)=colnames(testdata)[1:nplscomp]

  predicts<-predict(trees[[i]]$rpart,rmatrixnew2,...)

  testpre<-cbind(predicts,testpre)
}
  }

ensemble_predictions<-rowMeans(testpre)

return(ensemble_predictions)

}
