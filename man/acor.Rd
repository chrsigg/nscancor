\name{acor}
\alias{acor}
\title{Additional Explained Correlation}
\usage{
  acor(x, w, y, v, xcenter = TRUE, ycenter = TRUE,
    xscale = FALSE, yscale = FALSE)
}
\arguments{
  \item{x}{a numeric matrix which provides the data from
  the first domain}

  \item{w}{a numeric data matrix with the canonical vectors
  related to \code{x} as columns (returned as \code{xcoef}
  from \code{\link{nscancor}})}

  \item{y}{a numeric matrix which provides the data from
  the second domain}

  \item{v}{a numeric data matrix with the canonical vectors
  related to \code{y} as columns (returned as \code{ycoef}
  from \code{\link{nscancor}})}

  \item{xcenter}{a logical value indicating whether the
  empirical mean of (each column of) \code{x} should be
  subtracted. Alternatively, a vector of length equal to
  the number of columns of \code{x} can be supplied. The
  value is passed to \code{\link{scale}}.}

  \item{ycenter}{analogous to \code{xcenter}}

  \item{xscale}{a logical value indicating whether the
  columns of \code{x} should be scaled to have unit
  variance before the analysis takes place. The default is
  \code{FALSE} for consistency with \code{cancor}.
  Alternatively, a vector of length equal to the number of
  columns of \code{x} can be supplied.  The value is passed
  to \code{\link{scale}}.}

  \item{yscale}{analogous to \code{xscale}}
}
\description{
  \code{acor} computes the additional standard correlation
  explained by each canonical variable, taking into account
  the possible non-conjugacy of the canonical vectors in
  \eqn{\mathbf{W}}{W} and \eqn{\mathbf{V}}{V}.
}
\details{
  The additional correlation is measured after projecting
  the corresponding canonical vectors to the
  ortho-complement space spanned by the previous canonical
  variables. This procedure ensures that the correlation
  explained by non-conjugate canonical vectors is not
  counted multiple times.

  See Mackey (2009) for a presentation of generalized
  deflation in the context of principal component analysis
  (PCA), which was adapted here to CCA.
}
\references{
  Mackey, L. (2009) Deflation Methods for Sparse PCA. In
  \emph{Advances in Neural Information Processing Systems}
  (pp. 1017--1024).
}

