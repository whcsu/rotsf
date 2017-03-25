\name{rotsf}
\alias{rotsf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Grow a rotation survival random forest
}
\description{
Grow a rotation survival random forest
}
\usage{
rotsf(formula, data, trlength = 500, m = 2, control = control, na.action = na.omit, vari_status = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
%%     ~~Describe \code{formula} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{trlength}{
%%     ~~Describe \code{trlength} here~~
}
  \item{m}{
%%     ~~Describe \code{m} here~~
}
  \item{control}{
%%     ~~Describe \code{control} here~~
}
  \item{na.action}{
%%     ~~Describe \code{na.action} here~~
}
  \item{vari_status}{
%%     ~~Describe \code{vari_status} here~~
}
}
\details{
Grow a rotation survival random forest
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
function (formula, data, trlength = 500, m = 2, control = control, 
    na.action = na.omit, vari_status = FALSE) 
{
    Call <- match.call()
    data = na.omit(data)
    mf <- model.frame(formula, data)
    x <- model.matrix(attr(mf, "terms"), data = mf)
    y <- model.response(mf)
    if (!inherits(y, "Surv")) 
        stop("Response must be a 'survival' object - use the 'Surv()' function")
    ny <- ncol(y)
    n <- nrow(y)
    status <- y[, ny]
    survtime = y[, 1L]
    if (any(survtime <= 0)) 
        stop("Observation time must be > 0")
    if (all(status == 0)) 
        stop("No deaths in training data set")
    if (!missing(control)) 
        controls[names(control)] <- control
    p = dim(mf)[2]
    rotms <- vector(mode = "list", length = trlength)
    pectrees <- vector(mode = "list", length = trlength)
    varimp <- NULL
    for (i in 1:trlength) {
        trainindex = sample(nrow(mf), replace = T)
        trsetold = mf[trainindex, ]
        train_posp <- 1:nrow(mf) \%in\% trainindex
        oobset = trset[!train_posp, ]
        pp = c(1:p)
        rii = c(1)
        pp = setdiff(pp, rii)
        rp = p - length(rii)
        gr = (rp)/m
        varimport = c(rep(0, rp))
        rotationm = matrix(rep(0, (rp)^2), nrow = rp, ncol = rp)
        for (j in 1:gr) {
            trainindex = sample(nrow(trsetold), replace = T)
            trainb = trsetold[trainindex, ]
            d = sample(pp, m)
            olddata = trsetold[trainindex, d]
            rotj = prcomp(olddata)
            for (jj in 1:m) {
                rotationm[d - length(rii), (j - 1) * m + jj] = rotj$rotation[, 
                  jj]
            }
            pp = setdiff(pp, d)
        }
        dd = rotationm
        rotationm[pp - length(rii), ] = 0
        rotms[[i]] = rotationm
        rmatrix = as.matrix(trsetold[, -c(rii)]) \%*\% t(rotationm)
        rmatrixnew = as.data.frame(rmatrix)
        rmatrixnew = data.frame(trsetold[, rii][, 1], trsetold[, 
            rii][, 2], rmatrixnew)
        colnames(rmatrixnew)[c(1, 2)] = c("time", "status")
        colnames(rmatrixnew)[-c(1, 2)] = colnames(trsetold)[-c(1)]
        pectrees[[i]] = pecRpart(formula, data = rmatrixnew)
        if (vari_status) {
            rmatrix0 = as.matrix(oobset[, -c(rii)]) \%*\% t(rotationm)
            rmatrixnew0 = as.data.frame(rmatrix0)
            rmatrixnew0 = cbind(oobset[, rii], rmatrixnew0)
            colnames(rmatrixnew0) = colnames(oobset)
            predict_oob <- predict(trees[[i]], rmatrixnew0[, 
                -c(rii)])
            ci_before = concordance.index(predict_oob, oobset$time, 
                oobset$status)
            X = cbind(sample(oobset[, 3]), oobset[, c(4:p)])
            rmatrix1 = as.matrix(X) \%*\% t(rotationm)
            rmatrixnew0 = as.data.frame(rmatrix1)
            rmatrixnew0 = cbind(oobset[, rii], rmatrixnew0)
            colnames(rmatrixnew0) = colnames(oobset)
            predict_oob1 <- predict(trees[[i]], rmatrixnew0[, 
                -c(rii)])
            ci_after = concordance.index(predict_oob1, oobset$time, 
                oobset$status)
            varimport[1] = unlist(ci_before[1]) - unlist(ci_after[1])
            X = cbind(oobset[, 3:(p - 1)], sample(oobset[, p]))
            rmatrix1 = as.matrix(X) \%*\% t(rotationm)
            rmatrixnew0 = as.data.frame(rmatrix1)
            rmatrixnew0 = cbind(oobset[, rii], rmatrixnew0)
            colnames(rmatrixnew0) = colnames(oobset)
            predict_oob1 <- predict(trees[[i]], rmatrixnew0[, 
                -c(rii)])
            ci_after = concordance.index(predict_oob1, oobset$time, 
                oobset$status)
            varimport[p - 2] = unlist(ci_before[1]) - unlist(ci_after[1])
            for (jjj in 4:(p - 1)) {
                X = cbind(oobset[, 3:(jjj - 1)], sample(oobset[, 
                  jjj]), oobset[, (jjj + 1):p])
                rmatrix1 = as.matrix(X) \%*\% t(rotationm)
                rmatrixnew0 = as.data.frame(rmatrix1)
                rmatrixnew0 = cbind(oobset[, rii], rmatrixnew0)
                colnames(rmatrixnew0) = colnames(oobset)
                predict_oob1 <- predict(trees[[i]], rmatrixnew0[, 
                  -c(rii)])
                ci_after = concordance.index(predict_oob1, oobset$time, 
                  oobset$status)
                varimport[jjj - 2] = unlist(ci_before[1]) - unlist(ci_after[1])
            }
            if (vari_status) {
                varimp = cbind(varimp, varimport)
            }
        }
    }
    fit = pectrees
    if (vari_status) {
        vari = rowMeans(varimp)
        vari = vari * trlength
        colnames(vari) = colnames(mf[, -c(rii)])
    }
    class(fit) <- "rotsf"
    if (vari_status) 
        return(list(pectrees = pectrees, rotms = rotms, rii = rii, 
            vari = vari))
    else return(list(pectrees = pectrees, rotms = rotms, rii = rii))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
