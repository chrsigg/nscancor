\donttest{
if (requireNamespace("glmnet", quietly = TRUE) &&
    requireNamespace("PMA", quietly = TRUE)) {

  data(breastdata, package="PMA")

  set.seed(1)

  # Three data domains: a subset of genes, and CGH spots for the first and
  # second chromosome
  x <- with(breastdata,
            list(t(rna)[ , 1:100], t(dna)[ , chrom == 1], t(dna)[ , chrom == 2])
  )

  # Sparse regression functions with different cardinalities for different domains
  generate_predict <- function(dfmax) {
    force(dfmax)
    return(
      function(x, sc, cc) {
        en <- glmnet::glmnet(x, sc, alpha = 0.05, intercept = FALSE, dfmax = dfmax)
        W <- coef(en)
        return(W[2:nrow(W), ncol(W)])
      }
    )
  }
  predict <- lapply(c(20, 10, 10), generate_predict)

  # Compute two canonical variables per domain
  mcc <- mcancor(x, predict = predict, nvar = 2)

  # Compute another canonical variable for each domain
  mcc <- mcancor(x, predict = predict, nvar = 3, partial_model = mcc)
}
}
