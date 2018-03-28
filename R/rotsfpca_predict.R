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
rotsfpca.predict <-
function(rotsffit,newdata,trlength=500){

  trees=rotsffit$pectrees
  rotms=rotsffit$rotms
  rii=rotsffit$rii

  if (trlength>length(rotsffit$rotms))
    stop("Number of Trees for prediction should not be more than Number of Trees Fitted")

  # classify the test data
  testpre<-NULL
  for (i in 1:trlength) {
    #if (oobacc[i]<=avroobacc)
{

  # preparing for testing
  rmatrix2=as.matrix(newdata) %*% t(rotms[[i]])
  rmatrixnew2=as.data.frame(rmatrix2)

  colnames(rmatrixnew2)=colnames(newdata)

  predicts<-predict(trees[[i]]$rpart,rmatrixnew2)

  testpre<-cbind(predicts,testpre)
}
  }

ensemble_predictions<-rowMeans(testpre)

return(ensemble_predictions)

}
