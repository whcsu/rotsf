\name{rotsf.predict}
\alias{rotsf.predict}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prediction with new data and return a saved forest with mean hazard function
}
\description{
Prediction with new data and return a saved forest with mean hazard function
}
\usage{
rotsf.predict(rotsffit, newdata, trlength = 500)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rotsffit}{
%%     ~~Describe \code{rotsffit} here~~
}
  \item{newdata}{
%%     ~~Describe \code{newdata} here~~
}
  \item{trlength}{
%%     ~~Describe \code{trlength} here~~
}
}
\details{
Prediction with new data and return a saved forest with mean hazard function
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Zhou L, Xu Q, Wang H. (2015) Rotation survival forest for right censored data. PeerJ 3:e1009 https://doi.org/10.7717/peerj.1009
}
\author{
Hong Wang
}
\note{
Zhou L, Xu Q, Wang H. (2015) Rotation survival forest for right censored data. PeerJ 3:e1009 https://doi.org/10.7717/peerj.1009
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
Zhou L, Xu Q, Wang H. (2015) Rotation survival forest for right censored data. PeerJ 3:e1009 https://doi.org/10.7717/peerj.1009
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (rotsffit, newdata, trlength = 500) 
{
    trees = rotsffit$pectrees
    rotms = rotsffit$rotms
    rii = rotsffit$rii
    if (trlength > length(rotsf.fit$rotms)) 
        stop("Number of Trees for prediction should not be more than Number of Trees Fitted")
    testpre <- NULL
    for (i in 1:trlength) {
        {
            rmatrix2 = as.matrix(newdata) \%*\% t(rotms[[i]])
            rmatrixnew2 = as.data.frame(rmatrix2)
            colnames(rmatrixnew2) = colnames(newdata)
            predicts <- predict(trees[[i]]$rpart, rmatrixnew2)
            testpre <- cbind(predicts, testpre)
        }
    }
    ensemble_predictions <- rowMeans(testpre)
    return(ensemble_predictions)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
