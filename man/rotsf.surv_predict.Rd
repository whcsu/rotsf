\name{rotsf.surv_predict}
\alias{rotsf.surv_predict}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prediction with new data and return a saved forest with survival proability at Each unique time points
}
\description{
Prediction with new data and return a saved forest with survival proability at Each unique time points
}
\usage{
rotsf.surv_predict(rotsffit, newdata, uniquetimes, trlength = 500)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rotsffit}{
%%     ~~Describe \code{rotsffit} here~~
}
  \item{newdata}{
%%     ~~Describe \code{newdata} here~~
}
  \item{uniquetimes}{
%%     ~~Describe \code{uniquetimes} here~~
}
  \item{trlength}{
%%     ~~Describe \code{trlength} here~~
}
}
\details{
Prediction with new data and return a saved forest with survival proability at Each unique time points
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
function (rotsffit, newdata, uniquetimes, trlength = 500) 
{
    trees = rotsffit$pectrees
    rotms = rotsffit$rotms
    rii = rotsffit$rii
    if (trlength > length(rotsf.fit$rotms)) 
        stop("Number of Trees for prediction should not be more than Number of Trees Fitted")
    testpre <- matrix(0, nrow = dim(newdata)[1], ncol = length(uniquetimes))
    colnames(testpre) = paste0(uniquetimes)
    for (i in 1:trlength) {
        {
            rmatrix2 = as.matrix(newdata) \%*\% t(rotms[[i]])
            rmatrixnew2 = as.data.frame(rmatrix2)
            colnames(rmatrixnew2) = colnames(newdata)
            predicts = predictSurvProb(trees[[i]], rmatrixnew2, 
                uniquetimes)
            predicts[is.na(predicts)] = 0
            colnames(predicts) = paste0(uniquetimes)
            testpre <- testpre + predicts
        }
    }
    return(testpre/trlength)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
