##' Grow a rotation survival random forest
##' @title rotsf rostf_pca
##' @param x  The covariates(predictor variables) of training data.
##' @param y  Survival time and censored status of training data. Must be a Surv  \code{survival} object
##' @param testx  The covariates(predictor variables) of test data.
##' @param mtry   The number of covariates(predictor variables) used in each base ELM model. Default is the square root of the number of all avaibable covariates.
##' @param trlength  The ensemle size (the number of base ELM survival models). Default is 100.
##' @param Regularization_coefficient  Ridge or Tikhonov regularization parameter. Default is 10000. Also known as \eqn{C} in the ELM paper.
##' @param Kernel_type Type of kernel matrix. Currently four options avaibable. "RBF_kernel",a RBF kernel;"lin_kernel" , a linear kernel;poly_kernel ,a polynomial kernel;sigmoid_kernel, a sigmoid kernel. Default is "lin_kernel".
##' @param Kernel_para Parameters for different types of kernels. A single value for RBF and linear kernels. A vector for polynomial and sigmoid kernels and progam stops if only a single value is supplied. However, if the vector of values is supplied in the cases of RBF and liner kernels, only the first value will be used. Default is a vector value "c(2,1)"
##' @return Object of class \code{ELMSurvEN} with elements
##'   \tabular{ll}{
##'       \code{elmsurvfit}    \tab  A list of base models \code{elm_surv} of size \code{trlength}. To retrieve a particular base model: use  elmsurvfit[[i]], where i takes values between 1 and \code{trlength} \cr
##'       \code{precitedtime} \tab Esitmated survival times of test data. \cr
##'   }
##' @seealso \code{\link{elm_surv}}
##' @author Hong Wang
##' @references
##' \itemize{
##'   \item Zhou L, Xu Q, Wang H. (2015) Rotation survival forest for right censored data. PeerJ 3:e1009 https://doi.org/10.7717/peerj.1009.
##'  }
##' @examples
##' set.seed(123)
##' require(ELMSurv)
##' require(survival)
##' ## Survival Ensemble of ELM  with default settings
##' #Lung DATA
##' data(lung)
##' lung=na.omit(lung)
##' lung[,3]=lung[,3]-1
##' n=dim(lung)[1]
##' L=sample(1:n,ceiling(n*0.5))
##' trset<-lung[L,]
##' teset<-lung[-L,]
##' rii=c(2,3)
##' elmsurvmodel=ELMSurvEN(x=trset[,-rii],y=Surv(trset[,rii[1]], trset[,rii[2]]),testx=teset[,-c(rii)])
##' testpretimes=elmsurvmodel$precitedtime
##' #The predicted survival times on the first test example
##' head(testpretimes[1,])
##' #The predicted survival times of all test examples by the third model
##' head(testpretimes[,3])
##' # Get the 1th base model
##' firstbasemodel=elmsurvmodel$elmsurvfit[[1]]
##' @export
rrotsfpca<-function (x, y, trlength=500,m=2,mtry=floor(sqrt(ncol(x))), control = control, na.action =  na.omit,vari_status=FALSE)
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


  if (!missing(control))
    controls[names(control)] <- control


  #names(data)=c("time","status",names(x))




  rotms<- vector(mode = "list", length = trlength)

  pectrees <- vector(mode = "list", length = trlength)

  colindexes <- vector(mode = "list", length = trlength)

  varimp<-NULL

  for (i in 1:trlength)
  {
    #create a bagging version for rotation

    colindex=sample(ncol(x),size=mtry)
    colindexes[[i]]=colindex
    mf=data.frame(y,x[,colindex])
    p=dim(mf)[2]

    trainindex=sample(nrow(mf),replace=T)

    trsetold=mf[trainindex,]

    # oobset
    train_posp<-1:nrow(mf) %in% trainindex

    oobset=mf[!train_posp,]

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
    rmatrix=as.matrix(trsetold[,-c(rii)]) %*% t(rotationm)
    rmatrixnew=as.data.frame(rmatrix)
	
    rmatrixnew=data.frame(trsetold[,rii][,1],trsetold[,rii][,2],rmatrixnew)
    
	colnames(rmatrixnew)[c(1,2)]=c("time","status")
    colnames(rmatrixnew)[-c(1,2)]=colnames(x[,colindex])
   
      #trees[[i]]=rpart(rmatrixnew,control = control)
    #print(head(rmatrixnew))
   pectrees[[i]]=pecRpart(Surv(time,status)~.,data=rmatrixnew)

  }

  fit=pectrees

  class(fit) <- "rotsf"

  return(list(pectrees=pectrees,rotms=rotms,colindexes=colindexes,trlength=trlength))

}

