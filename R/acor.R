#  Copyright 2013 Christian Sigg
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

#' Additional Explained Correlation
#' 
#' \code{acor} computes the additional
#' standard correlation explained by each canonical component, taking into account 
#' the possible non-orthogonality of \eqn{\mathbf{W}}{W} and \eqn{\mathbf{V}}{V}.
#' 
#' The additional correlation of a component is measured after projecting the 
#' corresponding projection vectors to the ortho-complement space spanned by the 
#' previous canonical variables. This procedure ensures that the correlation explained
#' by non-orthogonal projection vectors is not counted multiple times. 
#' 
#' See Mackey (2009) for a presentation of generalized deflation in the context
#' of principal component analysis (PCA), which was adapted here to CCA.
#' 
#' @references Mackey, L. (2009) Deflation Methods for Sparse PCA. In 
#'   \emph{Advances in Neural Information Processing Systems} (pp. 1017--1024).
#'
#' @export
#' @param X a numeric matrix which provides the data from the first domain
#' @param W a numeric data matrix with the projection vectors related to
#'   \code{X} as columns (returned as \code{xcoef} from \code{\link{nscancor}})
#' @param Y a numeric matrix which provides the data from the second domain
#' @param V a numeric data matrix with the projection vectors related to
#'   \code{Y} as columns (returned as \code{ycoef} from \code{\link{nscancor}})
#' @param xcenter a logical value indicating whether the empirical mean of \code{X}
#'   should be subtracted. Alternatively, a vector of
#'   length equal the number of columns of \code{X} can be supplied.
#'   The value is passed to \code{\link{scale}}.
#' @param ycenter see \code{xcenter}
#' @param xscale a logical value indicating whether the columns of \code{X} should
#'   be scaled to have unit variance before the analysis takes
#'   place. The default is \code{FALSE} for consistency with \code{cancor}.
#'   Alternatively, a vector of length
#'   equal the number of columns of \code{X} can be supplied.  The
#'   value is passed to \code{\link{scale}}.
#' @param yscale see \code{xscale}
acor <- function(X, W, Y, V,  xcenter = TRUE, ycenter = TRUE, 
                 xscale = FALSE, yscale = FALSE) {
  
  dx <- ncol(X)
  dy <- ncol(Y)
  ncomp <- ncol(W)
  
  X <- scale(X, center = xcenter, scale = xscale)
  Y <- scale(Y, center = ycenter, scale = yscale)
  
  corr <- rep(0, ncomp)  # additional explained correlation
  Qx <- matrix(0, dx, ncomp)  # orthonormal basis spanned by the canonical variables X%*%W
  Qy <- matrix(0, dy, ncomp)  # orthonormal basis spanned by the canonical variables Y%*%V
  Xp <- X  # X projected to the orthocomplement space spanned by Qx
  Yp <- Y  # Y projected to the orthocomplement space spanned by Qy
  for (cc in seq_len(ncomp)) {
    w <- W[ ,cc]
    v <- V[ ,cc]
    corr[cc] <- cor(Xp%*%w, Yp%*%v)
    
    # update Qx and Qy
    XtXw <- t(X)%*%(X%*%w)
    YtYv <- t(Y)%*%(Y%*%v)
    if (cc > 1) {
      qx <- XtXw - Qx[ , 1:(cc-1)]%*%(t(Qx[ , 1:(cc-1)])%*%XtXw) 
      qy <- YtYv - Qy[ , 1:(cc-1)]%*%(t(Qy[ , 1:(cc-1)])%*%YtYv) 
    } else {
      qx <- XtXw
      qy <- YtYv
    }
    qx <- qx/normv(qx)
    qy <- qy/normv(qy)
    Qx[ , cc] <- qx
    Qy[ , cc] <- qy
    
    # deflate data matrices
    Xp <- Xp - Xp%*%qx%*%t(qx)  
    Yp <- Yp - Yp%*%qy%*%t(qy)
  }
  
  return(corr)
}