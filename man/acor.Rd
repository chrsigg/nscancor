\name{acor}
\alias{acor}
\title{Additional Explained Correlation}
\usage{
  acor(X, W, Y, V, xcenter = TRUE, ycenter = TRUE,
    xscale = FALSE, yscale = FALSE)
}
\arguments{
  \item{X}{a numeric matrix which provides the data from
  the first domain}

  \item{W}{a numeric data matrix with the projection
  vectors related to \code{X} as columns (returned as
  \code{xcoef} from \code{\link{nscancor}})}

  \item{Y}{a numeric matrix which provides the data from
  the second domain}

  \item{V}{a numeric data matrix with the projection
  vectors related to \code{Y} as columns (returned as
  \code{ycoef} from \code{\link{nscancor}})}

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
}
\description{
  \code{acor} computes the additional standard correlation
  explained by each canonical component, taking into
  account the possible non-orthogonality of
  \eqn{\mathbf{W}}{W} and \eqn{\mathbf{V}}{V}.
}
\details{
  The additional correlation of a component is measured
  after projecting the corresponding projection vectors to
  the ortho-complement space spanned by the previous
  canonical variables. This procedure ensures that the
  correlation explained by non-orthogonal projection
  vectors is not counted multiple times.

  See Mackey (2009) for a presentation of generalized
  deflation in the context of principal component analysis
  (PCA), which was adapted here to CCA.
}
\references{
  Mackey, L. (2009) Deflation Methods for Sparse PCA. In
  \emph{Advances in Neural Information Processing Systems}
  (pp. 1017--1024).
}

