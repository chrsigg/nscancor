library(MASS)
library(glmnet)
data(nutrimouse, package="CCA")

set.seed(1)

### Unconstrained CCA, identical to cancor(nutrimouse$gene[ , 1:10], nutrimouse$lipid)

ypredict = function(X, y, cc) {
  return(ginv(X)%*%y)
}
xpredict = function(Y, x, cc) {
  return(ginv(Y)%*%x)
} 
nscancor(nutrimouse$gene[ , 1:10], nutrimouse$lipid, xpredict=xpredict, 
         ypredict=ypredict)$cor


### Non-negative sparse CCA using glmnet

ypredict <- function(X, y, cc) {
    en <- glmnet(X, y, alpha=0.5, intercept=FALSE, dfmax=5, lower.limits=0)
    W <- coef(en)
    return(W[2:nrow(W), ncol(W)])
}
xpredict <- function(Y, x, cc) {
    en <- glmnet(Y, x, alpha=0.5, intercept=FALSE, dfmax=3, lower.limits=0)
    V <- coef(en)
    return(V[2:nrow(V), ncol(V)])
}
nscancor(nutrimouse$gene, nutrimouse$lipid, xpredict=xpredict, ncomp=5,
         ypredict=ypredict, verbosity=2)