#  Copyright 2013 Christian Sigg
#  Copyright 1995-2013 The R Core Team
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

#' Non-Negative and Sparse Multi-Domain CCA
#' 
#' Performs a canonical correlation analysis (CCA) on multiple data domains, 
#' where constraints such as non-negativity or sparsity are enforced on the 
#' canonical vectors.
#' 
#' \code{mcancor} generalizes \code{\link{nscancor}} to the case where more than
#' two data domains are available for an analysis. Its objective is to maximize 
#' the sum of all pairwise correlations of the canonical variables.
#' 
#' @export
#' @param x a list of numeric matrices which contain the data from the different
#'   domains
#' @param center a list of logical values indicating whether the empirical mean 
#'   of (each column of) the corresponding data matrix should be subtracted. 
#'   Alternatively, a list of vectors can be supplied, where each vector 
#'   specifies the mean to be subtracted from the corresponding data matrix. 
#'   Each list element is passed to \code{\link{scale}}.
#' @param scale_ a list of logical values indicating whether the columns of the 
#'   corresponding data matrix should be scaled to have unit variance before the
#'   analysis takes place. The default is \code{FALSE} for consistency with 
#'   \code{nscancor}. Alternatively, a list of vectors can be supplied, where 
#'   each vector specifies the standard deviations used to rescale the columns 
#'   of the corresponding data matrix. Each list element is passed to 
#'   \code{\link{scale}}.
#' @param nvar the number of canonical variables to be computed for each domain.
#'   With the default setting, canonical variables are computed until at least 
#'   one data matrix is fully deflated.
#' @param predict a list of regression functions to predict the sum of the 
#'   canonical variables of all other domains. The formal arguments for each 
#'   regression function are the design matrix \code{x} corresponding to the 
#'   data from the current domain, the regression target \code{sv} as the sum of
#'   the canonical variables for all other domains, and \code{vv} as a counter 
#'   of which canonical variable is currently computed (e.g. for enforcing 
#'   different constraints for subsequent canonical vectors of a given domain). 
#'   See the examples for an illustration.
#' @param cor_tol a threshold indicating the magnitude below which canonical 
#'   variables should be omitted. Variables are omitted if the sum of all their 
#'   correlations are less than or equal to \code{cor_tol} times the sum of all 
#'   correlations of the first canonical variables of all domains. With the 
#'   default \code{NULL} setting, no variables are omitted.
#' @param nrestart the number of random restarts for computing the canonical 
#'   variables via iterated regression steps. The solution achieving maximum 
#'   explained correlation over all random restarts is kept. A value greater 
#'   than one can help to avoid poor local maxima.
#' @param iter_tol If the relative change of the objective is less than 
#'   \code{iter_tol} between iterations, the procedure is asssumed to have 
#'   converged to a local optimum.
#' @param iter_max the maximum number of iterations to be performed. The 
#'   procedure is terminated if either the \code{iter_tol} or the 
#'   \code{iter_max} criterion is satisfied.
#' @param verbosity an integer specifying the verbosity level. Greater values 
#'   result in more output, the default is to be quiet.
#'   
#' @return \code{mcancor} returns a list with the following elements: 
#'   \item{cor}{a multi-dimensional array containing the additional correlations
#'   explained by each pair of canonical variables. The first two dimensions
#'   correspond to the domains, and the third dimension corresponds to the
#'   different canonical variables per domain (see also \code{\link{macor}}).} 
#'   \item{coef}{a list of matrices containing the canonical vectors related to 
#'   each data domain. The canonical vectors are stored as the columns of each 
#'   matrix.} \item{center}{the list of empirical means used to center the data
#'   matrices} \item{xscale}{the list of empirical standard deviations used to
#'   scale the data matrices}
#'   
#' @note Deflating the data matrices accumulates numerical errors over 
#'   successive canonical vectors.
#'   
#' @seealso  \code{\link{macor}}, \code{\link{nscancor}}, \code{\link{scale}}
#'   
#' @example inst/mcancor_examples.R
mcancor <- function (x, center = TRUE, scale_ = FALSE, nvar = NULL, predict,  
                     cor_tol = NULL, nrestart = 10, iter_tol = 1e-3, iter_max = 30,
                     verbosity = 0) {
  
  m <- length(x)  # number of domains
  
  X <- list();  # data sets
  cen <- list(); sc <- list();  # centering and scaling
  dx <- numeric(m)  # dimensionality of data domain
  for (mm in 1:m) {
    X[[mm]] <- scale(as.matrix(x[[mm]]), 
                     if (is.list(center)) center[[mm]] else center, 
                     if (is.list(scale_)) scale_[[mm]] else scale_
    )
    dx[mm] <- ncol(X[[mm]])

    cent <- attr(X[[mm]], "scaled:center")
    cen[[mm]] <- if(is.null(cent)) rep.int(0, dx[mm]) else cent  # follows cancor convention

    scal <- attr(X[[mm]], "scaled:scale")
    if(any(scal == 0))
      stop(paste("cannot rescale a constant column to unit variance in domain", mm))
    sc[[mm]] <- if(is.null(scal)) FALSE else scal
  }
  
  n <- nrow(X[[1]])  # number of observations
  if (is.null(nvar))
    nvar <- min(n, dx)
  
  corr <- array(NA, dim = c(m, m, nvar))  # additional explained correlation
  W <- list()  # canonical vectors 
  Q <- list()  # orthonormal basis spanned by the canonical variables
  for (mm in 1:m) {
    W[[mm]] <- matrix(NA, dx[mm], nvar)
    rownames(W[[mm]]) <- colnames(X[[mm]])
    Q[[mm]] <- matrix(NA, dx[mm], nvar)
  }
  Xp <- X  # X projected to the orthocomplement space spanned by Q
  
  for (vv in seq(nvar)) {
    obj_opt <- -Inf
    for (rr in seq(nrestart)) {
      
      res <- mcc_inner(X, Xp, dx, predict, vv, iter_tol, iter_max, verbosity)               
      
      # keep solution with maximum objective
      if (res$obj > obj_opt) {
        obj_opt <- res$obj
        w_opt <- res$w
        XpW <- res$XpW
      }
      if (verbosity > 0) {
        print(paste("canonical variable ", vv, ": ",
                    "maximum objective is ", format(obj_opt, digits = 4),
                    " at random restart ", rr-1, sep = ""))
      }
    }
    
    corr[ , , vv] <- cor(XpW, XpW)
    
    for (mm in 1:m) {
      w <- w_opt[[mm]]
      W[[mm]][ , vv] <- w
      rownames(W[[mm]]) <- colnames(X[[mm]])
      
      # update Q 
      XtXw <- t(X[[mm]])%*%(X[[mm]]%*%w)
      if (vv > 1) {
        q <- XtXw - Q[[mm]][ , 1:(vv-1)]%*%(t(Q[[mm]][ , 1:(vv-1)])%*%XtXw) 
      } else {
        q <- XtXw
      }
      q <- q/normv(q)
      Q[[mm]][ , vv] <- q
      
      # deflate data matrix
      Xp[[mm]] <- Xp[[mm]] - Xp[[mm]]%*%q%*%t(q)  
    }
    
    # current additionally explained correlation is below threshold cor_tol
    sum_corr_pp <- sum(corr[ , , vv] - diag(m))/2
    sum_corr_1 <- sum(corr[ , , 1] - diag(m))/2
    if (!is.null(cor_tol) && sum_corr_pp < cor_tol*sum_corr_1) {
      corr <- corr[ , , 1:(vv-1), drop=FALSE]
      for (mm in 1:m) {
        W[[mm]] <- W[[mm]][ , 1:(vv-1), drop=FALSE]  
      }
      break
    }    
    # at least one data matrix is fully deflated
    deflated <- logical(m)
    for (mm in 1:m) {
      deflated[mm] <- all(abs(Xp[[mm]]) < 1e-14)
    }
    if (vv < nvar && any(deflated)) { 
      if (verbosity > 0) {
        print("at least one data matrix is fully deflated, less than 'nvar' pairs of variables could be computed")            
      }
      corr <- corr[ , , 1:vv, drop=FALSE]
      for (mm in 1:m) {
        W[[mm]] <- W[[mm]][ , 1:vv, drop=FALSE]  
      }
      break
    }
  }
  
  return(list(cor = corr, coef = W, center = cen, scale = sc))
}

mcc_inner <- function(X, Xp, dx, predict, vv, iter_tol, iter_max, verbosity) {
  
  m <- length(X)  # number of domains
  n <- nrow(X[[1]])  # number of observations
  
  # initialize w as a random unit vector, non-negative if the prediction
  # function returns a non-negative vector
  w <- list()
  for (mm in 1:m) {
    v <- rnorm(dx[mm]);  
    if (all(predict[[mm]](Xp[[mm]], Xp[[mm]][ , 1], vv) >= 0)) {
      v <- abs(v)
    }
    w[[mm]] <- v/normv(Xp[[mm]]%*%v)   
  }    
  
  obj_old <- -Inf
  ii <- 0
  XpW <- matrix(NA, n, m)
  while(ii < iter_max) { 
    
    for (mm in 1:m) {
      XpW[ , mm] <- Xp[[mm]]%*%w[[mm]]
    }
    obj <- sum(t(XpW)%*%(XpW))
    if (abs(obj - obj_old)/obj < iter_tol) {
      break
    }
    obj_old <- obj
    
    for (mm in 1:m) {
      v <- predict[[m]](Xp[[mm]], rowSums(XpW[ , -mm, drop=FALSE]), vv)  
      if (all(v == 0))
        stop("w collapsed to the zero vector, try relaxing the constraints")
      
      w[[mm]] <- v/normv(Xp[[mm]]%*%v)
      XpW[ , mm] <- Xp[[mm]]%*%w[[mm]]
    }
    
    ii <- ii + 1
  }
  if ((verbosity > 0 ) && (ii == iter_max))
    print("maximum number of iterations reached before convergence")
  
  for (mm in 1:m) {
    w[[mm]] <- w[[mm]]/normv(X[[mm]]%*%w[[mm]])       
  }
  
  return(list(obj = obj, w = w, XpW = XpW))
}