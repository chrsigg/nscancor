\name{nscancor}
\alias{nscancor}
\title{Non-Negative Sparse CCA}
\usage{
  nscancor(x, y, xcenter = TRUE, ycenter = TRUE,
    xscale = FALSE, yscale = FALSE,
    ncomp = min(dim(x), dim(y)), xpredict, ypredict,
    cor_tol = NULL, nrestart = 10, iter_tol = 0.001,
    iter_max = 30, verbosity = 0)
}
\arguments{
  \item{x}{a numeric matrix which provides the data from
  the first domain}

  \item{y}{a numeric matrix which provides the data from
  the second domain}

  \item{xcenter}{a logical value indicating whether the
  empirical mean of \code{X} should be subtracted.
  Alternatively, a vector of length equal the number of
  columns of \code{X} can be supplied.  The value is passed
  to \code{\link{scale}}.}

  \item{ycenter}{see \code{xcenter}}

  \item{xscale}{a logical value indicating whether the
  columns of \code{X} should be scaled to have unit
  variance before the analysis takes place. The default is
  \code{FALSE} for consistency with \code{cancor}.
  Alternatively, a vector of length equal the number of
  columns of \code{X} can be supplied.  The value is passed
  to \code{\link{scale}}.}

  \item{yscale}{see \code{xscale}}

  \item{ncomp}{the number of canonical components to be
  computed. With the default setting, components are
  computed until either \code{x} or \code{y} is fully
  deflated.}

  \item{xpredict}{the regression function to predict the
  canonical variable of \code{x} given \code{y}. The formal
  arguments of the function are the design matrix \code{y},
  the current canonical variable of \code{x} as the target
  and the current component \code{cc}, e.g. for specifying
  different constraints for each component.}

  \item{ypredict}{see \code{xpredict}}

  \item{cor_tol}{a threshold indicating the magnitude below
  which components should be omitted. Components are
  omitted if their explained correlations are less than or
  equal to \code{cor_tol} times the correlation of the
  first component.  With the default \code{NULL} setting,
  no components are omitted.  With \code{cor_tol = 0} or
  \code{cor_tol = sqrt(.Machine$double.eps)}, essentially
  constant components are omitted.}

  \item{nrestart}{the number of random restarts for
  computing the projection pairs via iterated regression.
  The solution achieving maximum explained correlation over
  all random restarts is kept. A value greater than one can
  help to avoid bad local maxima.}

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
  \code{nscancor} returns a list containing the following
  elements: \item{cor}{the additional correlation explained
  by each component, see \code{\link{acor}}.}
  \item{xcoef}{the matrix containing the projection vectors
  related to \code{x} as its columns} \item{ycoef}{the
  matrix containing the projection vectors related to
  \code{y} as its columns} \item{xcenter}{if \code{xcenter}
  is \code{TRUE} the centering vector, else the zero vector
  } \item{xscale}{if \code{xscale} is \code{TRUE} the
  scaling vector, else FALSE } \item{ycenter}{see
  \code{xcenter}} \item{yscale}{see \code{xscale}}
}
\description{
  Performs a constrained canonical correlation analysis,
  where non-negativity, sparsity or other constraints are
  enforced on the projection vectors. The result of the
  analysis is returned as a list with the same elements as
  the list returned by \code{cancor}.
}
\details{
  \code{nscancor} computes the projection vectors using
  iterated regression steps, where the constraints suitable
  for each domain are enforced by choosing the appropriate
  regression method. See Sigg et al. (2007) for an early
  application of the principle (not yet including
  generalized deflation).

  Because constrained projection vectors no longer
  correspond to true eigenvectors of the cross-covariance
  matrix and are usually not pairwise conjugate, special
  attention needs to be paid when computing more than a
  single component. The algorithm implements a generalized
  deflation (GD) scheme which builds on GD for PCA as
  proposed by Mackey (2009). Given a basis of the space
  spanned by the canonical variables, the correlation of
  the component is maximized after projecting the current
  projection vector to the ortho-complement space of the
  basis. This procedure maximizes the additional
  correlation not explained by previous components, and is
  identical to standard CCA if no additional constraints
  are enforced.

  See the references for further details.
}
\note{
  Deflating the data matrices accumulates numerical errors
  over successive components.
}
\examples{
library(MASS)
library(glmnet)
data(nutrimouse, package="CCA")

set.seed(1)

### Unconstrained CCA, identical to cancor(nutrimouse$gene[ , 1:10], nutrimouse$lipid)

ypredict = function(X, y, cc) {
  return(ginv(X)\%*\%y)
}
xpredict = function(Y, x, cc) {
  return(ginv(Y)\%*\%x)
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

