library(glmnet)
data(breastdata, package="PMA")

set.seed(1)

# Three data domains: genes and CGH spots for the first and second chromosome
x <- with(breastdata, 
     list(t(rna), t(dna)[ , chrom==1], t(dna)[ , chrom==2])
)

# Sparse regression functions with different cardinalities for different domains
generate_predict <- function(dfmax) {
  force(dfmax)
  return(
    function(x, sv, vv) {
      en <- glmnet(x, sv, alpha=0.05, intercept=FALSE, dfmax=dfmax)
      W <- coef(en)
      return(W[2:nrow(W), ncol(W)])
    }
  )
}
predict <- lapply(c(20, 10, 10), generate_predict)

# Compute three canonical variables per domain
mcc <- mcancor(x, predict=predict, nvar=3, verbosity = 2)
