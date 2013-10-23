\name{nscancor}
\alias{nscancor}
\title{Non-Negative and Sparse CCA}
\usage{
  nscancor(x, y, xcenter = TRUE, ycenter = TRUE,
    xscale = FALSE, yscale = FALSE,
    npairs = min(dim(x), dim(y)), xpredict, ypredict,
    cor_tol = NULL, nrestart = 10, iter_tol = 0.001,
    iter_max = 30, verbosity = 0)
}
\arguments{
  \item{x}{a numeric matrix which provides the data from
  the first domain}

  \item{y}{a numeric matrix which provides the data from
  the second domain}

  \item{xcenter}{a logical value indicating whether the
  empirical mean of (each column of) \code{x} should be
  subtracted. Alternatively, a vector of length equal the
  number of columns of \code{x} can be supplied. The value
  is passed to \code{\link{scale}}.}

  \item{ycenter}{analogous to \code{xcenter}}

  \item{xscale}{a logical value indicating whether the
  columns of \code{x} should be scaled to have unit
  variance before the analysis takes place. The default is
  \code{FALSE} for consistency with \code{cancor}.
  Alternatively, a vector of length equal the number of
  columns of \code{x} can be supplied.  The value is passed
  to \code{\link{scale}}.}

  \item{yscale}{analogous to \code{xscale}}

  \item{npairs}{the number of pairs of canonical variables
  to be computed. With the default setting, pairs of
  variables are computed until either \code{x} or \code{y}
  is fully deflated.}

  \item{xpredict}{the regression function to predict the
  canonical variable for \code{x}, given \code{y}. The
  formal arguments are the design matrix \code{y}, the
  regression target \code{xv} as the current canonical
  variable for \code{x}, and \code{pp} as a counter of the
  current pair of canonical variables (e.g. for enforcing
  different constraints for different canonical vectors).
  See the examples for an illustration.}

  \item{ypredict}{analogous to \code{xpredict}}

  \item{cor_tol}{a threshold indicating the magnitude below
  which canonical variables should be omitted. Variables
  are omitted if their explained correlations are less than
  or equal to \code{cor_tol} times the correlation of the
  first pair of canonical variables. With the default
  \code{NULL} setting, no variables are omitted.}

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
  \code{nscancor} returns a list with the following
  elements: \item{cor}{the additional correlation explained
  by each pair of canonical variables, see
  \code{\link{acor}}.} \item{xcoef}{the matrix containing
  the canonical vectors related to \code{x} as its columns}
  \item{ycoef}{analogous to \code{xcoef}} \item{xcenter}{if
  \code{xcenter} is \code{TRUE} the centering vector, else
  the zero vector (in accordance with \code{cancor})}
  \item{ycenter}{analogous to \code{xcenter}}
  \item{xscale}{if \code{xscale} is \code{TRUE} the scaling
  vector, else FALSE } \item{yscale}{analogous to
  \code{xscale}}
}
\description{
  Performs a canonical correlation analysis(CCA) where
  constraints such as non-negativity or sparsity are
  enforced on the canonical vectors. The result of the
  analysis is returned as a list with the same elements as
  the list returned by \code{cancor}.
}
\details{
  \code{nscancor} computes the canonical vectors (called
  \code{xcoef} and \code{ycoef}) using iterated regression
  steps, where the constraints suitable for each domain are
  enforced by choosing the appropriate regression method.
  See Sigg et al. (2007) for an early application of the
  principle (not yet including generalized deflation).

  Because constrained canonical vectors no longer
  correspond to true eigenvectors of the cross-covariance
  matrix and are usually not pairwise conjugate (i.e. the
  canonical variables are not uncorrelated), special
  attention needs to be paid when computing more than a
  single pair of canonical vectors. \code{nscancor}
  implements a generalized deflation (GD) scheme which
  builds on GD for PCA as proposed by Mackey (2009). For
  each domain, a basis of the space spanned by the previous
  canonical variables is computed. Then, the correlation of
  the current pair of canonical variables is maximized
  after projecting each current canonical vector to the
  ortho-complement space of its respective basis. This
  procedure maximizes the additional correlation not
  explained by previous canonical variables, and is
  identical to standard CCA if the canonical vectors are
  the eigenvectors of the cross-covariance matrix.

  See the references for further details.
}
\note{
  Deflating the data matrices accumulates numerical errors
  over successive canonical vectors.
}
\examples{
library(MASS)
library(glmnet)
data(nutrimouse, package="CCA")

set.seed(1)

### Unconstrained CCA, produces identical results to calling 
# cancor(nutrimouse$gene[ , 1:10], nutrimouse$lipid)

ypredict = function(X, yv, vv) {
  return(ginv(X)\%*\%yv)
}
xpredict = function(Y, xv, vv) {
  return(ginv(Y)\%*\%xv)
} 
nscancor(nutrimouse$gene[ , 1:10], nutrimouse$lipid, xpredict=xpredict, 
         ypredict=ypredict)


### Non-negative sparse CCA using glmnet() as the regression function

ypredict <- function(X, yv, vv) {
    en <- glmnet(X, yv, alpha=0.5, intercept=FALSE, dfmax=5, lower.limits=0)
    W <- coef(en)
    return(W[2:nrow(W), ncol(W)])
}
xpredict <- function(Y, xv, vv) {
    en <- glmnet(Y, xv, alpha=0.5, intercept=FALSE, dfmax=3, lower.limits=0)
    V <- coef(en)
    return(V[2:nrow(V), ncol(V)])
}
nscancor(nutrimouse$gene, nutrimouse$lipid, xpredict=xpredict, npairs=5,
         ypredict=ypredict, verbosity=2)
}
\references{
  Sigg, C. and Fischer, B. and Ommer, B. and Roth, V. and
  Buhmann, J. (2007) Nonnegative CCA for Audiovisual Source
  Separation. In \emph{Proceedings of the 2007 IEEE
  Workshop on Machine Learning for Signal Processing} (pp.
  253--258).

  Mackey, L. (2009) Deflation Methods for Sparse PCA. In
  \emph{Advances in Neural Information Processing Systems}
  (pp. 1017--1024).
}
\seealso{
  \code{\link{acor}}, \code{\link{cancor}},
  \code{\link{scale}}
}

