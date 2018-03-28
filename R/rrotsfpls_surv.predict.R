# -------------------------------------------------------------------------------
#   This file is part of Rotsf.
#
# Rotsf is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Rotsf is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Rotsf. If not, see <http://www.gnu.org/licenses/>.
#
# Written by:
#
#   Hong Wang
# 
# http://www.csu.edu.cn
# wh@csu.edu.cn
# -------------------------------------------------------------------------------

##' Prediction with new data and return a saved forest with mean hazard function 
##' @export
rrotsfpls_surv.predict<-function(object,newdata,uniquetimes,trlength=object$trlength){

  trees=object$pectrees
  rotms=object$rotms
  colindexes=object$colindexes
  nplscomp=object$nplscomp

  if (trlength>object$trlength)
    stop("Number of Trees for prediction should not be more than Number of Trees Fitted")

  # classify the test data
   # Preparing testpre dataframe 
  testpre <- matrix(0,nrow = dim(newdata)[1],  ncol= length(uniquetimes))
  colnames(testpre)=paste0(uniquetimes)
  
  for (i in 1:trlength) {
    #if (oobacc[i]<=avroobacc)
{

  # preparing for testing
  
  testdata=newdata[,colindexes[[i]]]
  rmatrix2=as.matrix(testdata) %*% (rotms[[i]])
  rmatrixnew2=as.data.frame(rmatrix2)

  #colnames(rmatrixnew2)=colnames(testdata)
  colnames(rmatrixnew2)=colnames(testdata)[1:nplscomp]
  
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
