##' Grow a survival forest with partial least squares(PLS)
##' @title rotsf rrotsfspls
##' @param x  The covariates(predictor variables) of training data.
##' @param y  Survival time and censored status of training data. Must be a Surv  \code{survival} object
##' @param mtry   The number of covariates(predictor variables) used in each base tree model. Default is the square root of the number of all avaibable covariates.
##' @param trlength  The ensemle size (the number of base  survival trees). Default is 100.
##' @param vari_status  Whether or not calculate variables importance scores. Default is "FALSE".
##' @param ... Additional arguments for the base decision tree, see the \code{rpart} package for details.
##' @return Object of class \code{rrotsfspls} with elements
##'   \tabular{ll}{
##'       \code{pectrees}    \tab  A list of base models \code{pecRpart} in \code{pec} R package of size \code{trlength}. To retrieve a particular base model: use  pectrees[[i]], where i takes values between 1 and \code{trlength} \cr
##'       \code{colindexes} \tab A list of covaraite subspace index for each base tree. \cr
##'       \code{trlength} \tab Number of bases models trained. \cr
##'       \code{rotms} \tab A list of PLS weights  matrix.To retrieve a particular weight matrix: use  rotms[[i]], where i takes values between 1 and \code{trlength} \cr
##'       \code{varimp} \tab If  \code{vari_status=FALSE}, return a Matrix of variable importance scores. \cr
##'   }
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
##' # Get the 1th base model
##' firstbasemodel=plssurvmodel$pectrees[[1]]
##' #second PLS weight matrix
##' secondweigmatrix=plssurvmodel$rotms[[2]]
##' plssurvmodel2=rrotsfspls(x=trset[,-rii],y=Surv(trset[,rii[1]], trset[,rii[2]]),vari_status=TRUE)
##' #variable importance
##' varimpscores=plssurvmodel2$varimp
##' @export
rrotsfspls <-
function (x, y, trlength=100,mtry=floor(sqrt(ncol(x))),impute=TRUE, nplscomp=floor(sqrt(ncol(x))),vari_status=FALSE,...)
{
  Call <- match.call()


  if (!inherits(y, "Surv"))
    stop("Response must be a 'survival' object - use the 'Surv()' function")

  ny <- ncol(y)
  n <- nrow(y)

  status <- y[, ny]
  survtime=y[, 1L]



  if (any(survtime <= 0)) stop("Observation time must be > 0")

  if (all(status == 0)) stop("No deaths in training data set")



  #names(data)=c("time","status",names(x))


  rotms<- vector(mode = "list", length = trlength)
  colindexes <- vector(mode = "list", length = trlength)

  pectrees <- vector(mode = "list", length = trlength)
  ptm <- proc.time()
  unusedtr=0
  varimp=Matrix(0, nrow = ncol(x), ncol = 1, sparse = TRUE)
  for (i in 1:trlength)
  {
   

   #create a bagging version for rotation

    colindex=sample(ncol(x),size=mtry)
    colindexes[[i]]=colindex
	
	mf=data.frame(survtime,status,x[,colindex])
    
	p=dim(mf)[2]-2
    
	#kval=max(p+1,nrow(mf))
	
    trainindex=sample(nrow(mf),replace=T)

    bagtrset=mf[trainindex,]
	
	rii=c(1,2)
    colnames(bagtrset)[rii]=c("time","status")
	colnames(bagtrset)[-rii]=colnames(x)[colindex]
	#buckely-james imputation
	
	#imputey
 if(impute){
    newy=bjimpute(y=bagtrset[,1], cen=bagtrset[,2], x=bagtrset[,-c(1,2)],inibeta=NULL)
    
	#replace the old y with new y
    trsetold=data.frame(newy,bagtrset[,-c(1)])
	}else{
	trsetold=bagtrset
	}
    colnames(trsetold)=colnames(bagtrset)
		
    #initialize the rotation matrix
		
    plsres=pls::simpls.fit(as.matrix(trsetold[,-c(1,2)]),trsetold[,1],ncomp=nplscomp)
	
    if (vari_status)
	{
		SS <- c(plsres$Yloadings)^2 * colSums(plsres$scores^2)
		Wnorm2 <- colSums(plsres$projection^2)
		SSW <- sweep(plsres$projection^2, 2, SS / Wnorm2, "*")
		vip= sqrt(nrow(SSW) * apply(SSW, 1, cumsum) / cumsum(SS))
		vipscore=colMeans(vip)	
	    varimp[colindex]=varimp[colindex]+vipscore
	}
	
	
	# the final the rotation matrix
	 rotms[[i]]=plsres$projection
    
    #print(dim(rotms[[i]]))
	#print(dim(trsetold[,-c(rii)]))
    # the final bagging training set
    rmatrix=as.matrix(trsetold[,-c(rii)]) %*% (rotms[[i]])
    rmatrixnew=as.data.frame(rmatrix)
    
    rmatrixnew=data.frame(trsetold[,rii][,1],trsetold[,rii][,2],rmatrixnew)

    colnames(rmatrixnew)=colnames(bagtrset)[1:nplscomp]

	#colnames(rmatrixnew)[rii]=c("time","status")
    
  
      #trees[[i]]=rpart(rmatrixnew,control = control)
    #print(head(rmatrixnew))
    fixed_col_names=paste(colnames(rmatrixnew)[3:nplscomp], collapse='+')
    pectrees[[i]]=pecRpart(as.formula(paste("Surv(time,status)","~",fixed_col_names)),data=rmatrixnew)
    #pectrees[[i]]=pecRpart(Surv(time,status)~.,data=rmatrixnew)
    
    
	# if ((i%%20==0) && (i<trlength))	{
	 # runned=unlist((proc.time() - ptm)[1])	
	 # totaltime=runned*trlength/i
	 # #add cat to print out
	 # if (totaltime>5)
	 # cat(sprintf("\t %.0f%% compeleted %.0f minutes remaining \n ",i/trlength*100,  (totaltime-runned)/60))
	 #}
	 
}
  

  fit=pectrees

  class(fit) <- "rrotsfspls"

  if (vari_status)
  { 
    #normalize the scores
    varimp=varimp/trlength*100
    rownames(varimp)=colnames(x)
}


	fit <- list()
    fit$pectrees=pectrees
	fit$rotms=rotms
	fit$colindexes=colindexes
	fit$trlength=trlength
	fit$nplscomp=nplscomp 


if (vari_status){
   	fit$varimp=varimp
	class(fit) <- "rrotsfspls"
    fit
  }
else{
	class(fit) <- "rrotsfspls"
    fit
	}
}
