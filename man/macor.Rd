\name{macor}
\alias{macor}
\title{Multi-Domain Additional Explained Correlation}
\usage{
  macor(x, w, center = TRUE, scale_ = FALSE)
}
\arguments{
  \item{x}{a list of numeric matrices which contain the
  data from the different domains}

  \item{w}{a list of matrices containing the canonical
  vectors related to each data domain. The canonical
  vectors must be stored as the columns of each matrix.}

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
  \code{acor}. Alternatively, a list of vectors can be
  supplied, where each vector specifies the standard
  deviations used to rescale the columns of the
  corresponding data matrix. Each list element is passed to
  \code{\link{scale}}.}
}
\value{
  \code{macor} returns a multi-dimensional array containing
  the additional correlations explained by each pair of
  canonical variables. The first two dimensions correspond
  to the domains, and the third dimension corresponds to
  the different canonical variables per domain.
}
\description{
  \code{macor} generalizes \code{\link{acor}} to the case
  of more than two data domains.
}

