library(MASS)
library(glmnet)
data(nutrimouse, package="CCA")

set.seed(1)

### 
# Unconstrained CCA, produces identical results to calling 
# cancor(nutrimouse$gene[ , 1:10], nutrimouse$lipid)

ypredict <- function(x, yv, vv) {
  return(ginv(x)%*%yv)
}
xpredict <- function(y, xv, vv) {
  return(ginv(y)%*%xv)
} 
cc <- nscancor(nutrimouse$gene[ , 1:10], nutrimouse$lipid, xpredict=xpredict, 
               ypredict=ypredict)


### 
# Non-negative sparse CCA using glmnet() as the regression function, where
# different regularisers are enforced on the different data domains and pairs of 
# canonical variables.

dfmax_w <- c(40, 15, 10, 10)
ypredict <- function(x, yv, vv) {
  en <- glmnet(x, yv, alpha=0.5, intercept=FALSE, dfmax=dfmax_w[vv], lower.limits=0)
  W <- coef(en)
  return(W[2:nrow(W), ncol(W)])
}
dfmax_v <- c(7, 5, 5, 3)
xpredict <- function(y, xv, vv) {
  en <- glmnet(y, xv, alpha=0.5, intercept=FALSE, dfmax=dfmax_v[vv])
  V <- coef(en)
  return(V[2:nrow(V), ncol(V)])
}
nscc <- nscancor(nutrimouse$gene, nutrimouse$lipid, nvar=3,
                 xpredict=xpredict, ypredict=ypredict)

# continue the computation from a partial model
nscc <- nscancor(nutrimouse$gene, nutrimouse$lipid, nvar=4,
                 xpredict=xpredict, ypredict=ypredict,
                 partial_model=nscc)
