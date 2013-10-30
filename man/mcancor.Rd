\name{mcancor}
\alias{mcancor}
\title{Non-Negative and Sparse Multi-Domain CCA}
\usage{
  mcancor(x, center = TRUE, scale_ = FALSE, nvar = NULL,
    predict, cor_tol = NULL, nrestart = 10,
    iter_tol = 0.001, iter_max = 30, verbosity = 0)
}
\arguments{
  \item{x}{a list of numeric matrices which contain the
  data from the different domains}

  \item{center}{a list of logical values indicating whether
  the empirical mean of (each column of) the corresponding
  data matrix should be subtracted.  Alternatively, a list
  of vectors can be supplied, where each vector specifies
  the mean to be subtracted from the corresponding data
  matrix.  Each list element is passed to
  \code{\link{scale}}.}

  \item{scale_}{a list of logical values indicating whether
  the columns of the corresponding data matrix should be
  scaled to have unit variance before the analysis takes
  place. The default is \code{FALSE} for consistency with
  \code{nscancor}. Alternatively, a list of vectors can be
  supplied, where each vector specifies the standard
  deviations used to rescale the columns of the
  corresponding data matrix. Each list element is passed to
  \code{\link{scale}}.}

  \item{nvar}{the number of canonical variables to be
  computed for each domain.  With the default setting,
  canonical variables are computed until at least one data
  matrix is fully deflated.}

  \item{predict}{a list of regression functions to predict
  the sum of the canonical variables of all other domains.
  The formal arguments for each regression function are the
  design matrix \code{x} corresponding to the data from the
  current domain, the regression target \code{sv} as the
  sum of the canonical variables for all other domains, and
  \code{vv} as a counter of which canonical variable is
  currently computed (e.g. for enforcing different
  constraints for subsequent canonical vectors of a given
  domain).  See the examples for an illustration.}

  \item{cor_tol}{a threshold indicating the magnitude below
  which canonical variables should be omitted. Variables
  are omitted if the sum of all their correlations are less
  than or equal to \code{cor_tol} times the sum of all
  correlations of the first canonical variables of all
  domains. With the default \code{NULL} setting, no
  variables are omitted.}

  \item{nrestart}{the number of random restarts for
  computing the canonical variables via iterated regression
  steps. The solution achieving maximum explained
  correlation over all random restarts is kept. A value
  greater than one can help to avoid poor local maxima.}

  \item{iter_tol}{If the relative change of the objective
  is less than \code{iter_tol} between iterations, the
  procedure is asssumed to have converged to a local
  optimum.}

  \item{iter_max}{the maximum number of iterations to be
  performed. The procedure is terminated if either the
  \code{iter_tol} or the \code{iter_max} criterion is
  satisfied.}

  \item{verbosity}{an integer specifying the verbosity
  level. Greater values result in more output, the default
  is to be quiet.}
}
\value{
  \code{mcancor} returns a list with the following
  elements: \item{cor}{a multi-dimensional array containing
  the additional correlations explained by each pair of
  canonical variables. The first two dimensions correspond
  to the domains, and the third dimension corresponds to
  the different canonical variables per domain (see also
  \code{\link{macor}}).} \item{coef}{a list of matrices
  containing the canonical vectors related to each data
  domain. The canonical vectors are stored as the columns
  of each matrix.} \item{center}{the list of empirical
  means used to center the data matrices} \item{xscale}{the
  list of empirical standard deviations used to scale the
  data matrices}
}
\description{
  Performs a canonical correlation analysis (CCA) on
  multiple data domains, where constraints such as
  non-negativity or sparsity are enforced on the canonical
  vectors.
}
\details{
  \code{mcancor} generalizes \code{\link{nscancor}} to the
  case where more than two data domains are available for
  an analysis. Its objective is to maximize the sum of all
  pairwise correlations of the canonical variables.
}
\note{
  Deflating the data matrices accumulates numerical errors
  over successive canonical vectors.
}
\examples{
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
}
\seealso{
  \code{\link{macor}}, \code{\link{nscancor}},
  \code{\link{scale}}
}

