library(MASS)
library(glmnet)
data(nutrimouse, package="CCA")

set.seed(1)

### Unconstrained CCA, produces identical results to calling 
# cancor(nutrimouse$gene[ , 1:10], nutrimouse$lipid)

ypredict <- function(x, yv, vv) {
  return(ginv(x)%*%yv)
}
xpredict <- function(y, xv, vv) {
  return(ginv(y)%*%xv)
} 
cc <- nscancor(nutrimouse$gene[ , 1:10], nutrimouse$lipid, xpredict=xpredict, 
         ypredict=ypredict)


### Non-negative sparse CCA using glmnet() as the regression function

ypredict <- function(x, yv, vv) {
    en <- glmnet(x, yv, alpha=0.5, intercept=FALSE, dfmax=5, lower.limits=0)
    W <- coef(en)
    return(W[2:nrow(W), ncol(W)])
}
xpredict <- function(y, xv, vv) {
    en <- glmnet(y, xv, alpha=0.5, intercept=FALSE, dfmax=3, lower.limits=0)
    V <- coef(en)
    return(V[2:nrow(V), ncol(V)])
}
nscc <- nscancor(nutrimouse$gene, nutrimouse$lipid, xpredict=xpredict, nvar=5,
         ypredict=ypredict)
