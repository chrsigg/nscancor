library(MASS)
library(glmnet)
data(nutrimouse, package="CCA")

set.seed(1)

### Unconstrained CCA, produces identical results to calling 
# cancor(nutrimouse$gene[ , 1:10], nutrimouse$lipid)

ypredict <- function(X, yv, pp) {
  return(ginv(X)%*%yv)
}
xpredict <- function(Y, xv, pp) {
  return(ginv(Y)%*%xv)
} 
nscancor(nutrimouse$gene[ , 1:10], nutrimouse$lipid, xpredict=xpredict, 
         ypredict=ypredict)


### Non-negative sparse CCA using glmnet() as the regression function

ypredict <- function(X, yv, pp) {
    en <- glmnet(X, yv, alpha=0.5, intercept=FALSE, dfmax=5, lower.limits=0)
    W <- coef(en)
    return(W[2:nrow(W), ncol(W)])
}
xpredict <- function(Y, xv, pp) {
    en <- glmnet(Y, xv, alpha=0.5, intercept=FALSE, dfmax=3, lower.limits=0)
    V <- coef(en)
    return(V[2:nrow(V), ncol(V)])
}
nscancor(nutrimouse$gene, nutrimouse$lipid, xpredict=xpredict, npairs=5,
         ypredict=ypredict, verbosity=2)