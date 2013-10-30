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

#' Multi-Domain Additional Explained Correlation
#' 
#' \code{macor} generalizes \code{\link{acor}} to the case of more than two
#' data domains.
#' 
#' @export
#' @param x a list of numeric matrices which contain the data from the different
#'   domains
#' @param w a list of matrices containing the canonical vectors related to 
#'   each data domain. The canonical vectors must be stored as the columns of each 
#'   matrix.
#' @param center a list of logical values indicating whether the empirical mean 
#'   of (each column of) the corresponding data matrix should be subtracted. 
#'   Alternatively, a list of vectors can be supplied, where each vector 
#'   specifies the mean to be subtracted from the corresponding data matrix. 
#'   Each list element is passed to \code{\link{scale}}.
#' @param scale_ a list of logical values indicating whether the columns of the 
#'   corresponding data matrix should be scaled to have unit variance before the
#'   analysis takes place. The default is \code{FALSE} for consistency with 
#'   \code{acor}. Alternatively, a list of vectors can be supplied, where each
#'   vector specifies the standard deviations used to rescale the columns of the
#'   corresponding data matrix. Each list element is passed to 
#'   \code{\link{scale}}.
#' @return \code{macor} returns a multi-dimensional array containing the
#'   additional correlations explained by each pair of canonical variables. The
#'   first two dimensions correspond to the domains, and the third dimension
#'   corresponds to the different canonical variables per domain.
macor <- function(x, w, center = TRUE, scale_ = FALSE) {
  
  X <- x
  W <- w
  m <- length(X)  # number of domains
  n <- nrow(X[[1]])  # number of observations
  nvar <- ncol(W[[1]])  # number of canonical variables for each domain
  
  dx <- numeric(m)  # dimensionality of data domain
  for (mm in 1:m) {
    X[[mm]] <- scale(X[[mm]], 
                     if (is.list(center)) center[[mm]] else center, 
                     if (is.list(scale_)) scale_[[mm]] else scale_
    )
    dx[mm] <- ncol(X[[mm]])
  }
  
  corr <- array(NA, dim = c(m, m, nvar))  # additional explained correlation
  Q <- list()  # orthonormal basis spanned by the canonical variables
  for (mm in 1:m) {
    Q[[mm]] <- matrix(NA, dx[mm], nvar)
  }
  Xp <- X  # X projected to the orthocomplement space spanned by Q
  
  for (pp in seq(nvar)) {

    XpW <- matrix(NA, n, m)  
    for (mm in 1:m) {
      w <- W[[mm]][ , pp]
      
      XpW[ , mm] <- Xp[[mm]]%*%w
      
      # update Q 
      XtXw <- t(X[[mm]])%*%(X[[mm]]%*%w)
      if (pp > 1) {
        q <- XtXw - Q[[mm]][ , 1:(pp-1)]%*%(t(Q[[mm]][ , 1:(pp-1)])%*%XtXw) 
      } else {
        q <- XtXw
      }
      q <- q/normv(q)
      Q[[mm]][ , pp] <- q
      
      # deflate data matrix
      Xp[[mm]] <- Xp[[mm]] - Xp[[mm]]%*%q%*%t(q)  
    }
    
    corr[ , , pp] <- cor(XpW, XpW)
  }
    
  return(corr)
}