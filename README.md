# Survival Forest with PCA and Partial Least Squares

We investigate the plausibility of extending a rotation forest, originally proposed for classification purpose, to survival analysis. Currently, PCA(Principal components analysis) , SPLS(Sparse partial least squares) rotation are supported.

For a detailed information, see the papers by Zhou L, Xu Q, Wang H. (2015) Rotation survival forest for right censored data. <https://doi.org/10.7717/peerj.1009> and Zhou L, Xu Q, Wang H. (2016), Random rotation survival forest for high dimensional censored data. <https://doi.org/10.1186/s40064-016-3113-5>. Zhou L,  Wang H. Xu Q, (2018) Survival Forest with Partial Least Squares For High Dimensional Censored Data https://doi.org/10.1016/j.chemolab.2018.05.005  

R version >= 3.1 and the latest new Rtools toolchain need to be installed to compile the package.

Download the package via github

install.packages("devtools") # if you have not installed "devtools" package

devtools::install_github("whcsu/rotsf")
