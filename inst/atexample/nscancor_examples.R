\donttest{
\dontshow{
if (requireNamespace("CCA", quietly = TRUE) &&
    requireNamespace("glmnet", quietly = TRUE) &&
    requireNamespace("MASS", quietly = TRUE)) \{
}
data(nutrimouse, package = "CCA")

###
# Unconstrained CCA, produces results close to calling
# cancor(nutrimouse$gene[ , 1:5], nutrimouse$lipid)

ypredict <- function(x, yc, cc) {
  return(MASS::ginv(x)%*%yc)
}
xpredict <- function(y, xc, cc) {
  return(MASS::ginv(y)%*%xc)
}
cc <- nscancor(nutrimouse$gene[ , 1:5], nutrimouse$lipid, xpredict = xpredict,
               ypredict = ypredict)
print(cc$cor)

###
# Non-negative sparse CCA using glmnet() as the regression function, where
# different regularizers are enforced on the different data domains and pairs
# of canonical variables.

dfmax_w <- c(40, 15, 10)
ypredict <- function(x, yc, cc) {
  en <- glmnet::glmnet(x, yc, alpha = 0.5, intercept = FALSE,
                       dfmax = dfmax_w[cc], lower.limits = 0)
  W <- coef(en)
  return(W[2:nrow(W), ncol(W)])
}
dfmax_v <- c(7, 5, 5)
xpredict <- function(y, xc, cc) {
  en <- glmnet::glmnet(y, xc, alpha = 0.5, intercept = FALSE,
                       dfmax = dfmax_v[cc])
  V <- coef(en)
  return(V[2:nrow(V), ncol(V)])
}
nscc <- nscancor(nutrimouse$gene, nutrimouse$lipid, nvar = 2,
                 xpredict = xpredict, ypredict = ypredict)

# continue the computation of canonical variables from a partial model
nscc <- nscancor(nutrimouse$gene, nutrimouse$lipid, nvar = 3,
                 xpredict = xpredict, ypredict = ypredict,
                 partial_model = nscc)
print(nscc$cor)
print(nscc$xcoef)
print(nscc$ycoef)
\dontshow{\}}
}
