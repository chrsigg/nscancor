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

#' Non-Negative and Sparse CCA
#' 
#' Performs a canonical correlation analysis(CCA) where constraints such as 
#' non-negativity or  sparsity are enforced on the canonical vectors. The result
#' of the analysis is returned as a list with the same elements as the list 
#' returned by \code{cancor}.
#' 
#' \code{nscancor} computes the canonical vectors (called \code{xcoef} and 
#' \code{ycoef}) using iterated regression steps, where the constraints suitable
#' for each domain are enforced by choosing the appropriate regression method. 
#' See Sigg et al. (2007) for an early application of the principle (not yet 
#' including generalized deflation).
#' 
#' Because constrained canonical vectors no longer correspond to true 
#' eigenvectors of the cross-covariance matrix and are usually not pairwise 
#' conjugate (i.e. the canonical variables are not uncorrelated), special 
#' attention needs to be paid when computing more than a single pair of 
#' canonical vectors. \code{nscancor} implements a generalized deflation (GD) 
#' scheme which builds on GD for PCA as proposed by Mackey (2009). For each 
#' domain, a basis of the space spanned by the previous canonical variables is 
#' computed. Then, the correlation of the current pair of canonical variables is
#' maximized after projecting each current canonical vector to the 
#' ortho-complement space of its respective basis. This procedure maximizes the 
#' additional correlation not explained by previous canonical variables, and is 
#' identical to standard CCA if the canonical vectors are the eigenvectors of 
#' the cross-covariance matrix.
#' 
#' See the references for further details.
#' 
#' @export nscancor
#' @param x a numeric matrix which provides the data from the first domain
#' @param y a numeric matrix which provides the data from the second domain
#' @param xcenter a logical value indicating whether the empirical mean of (each
#'   column of) \code{x} should be subtracted. Alternatively, a vector of length
#'   equal the number of columns of \code{x} can be supplied. The value is
#'   passed to \code{\link{scale}}.
#' @param ycenter analogous to \code{xcenter}
#' @param xscale a logical value indicating whether the columns of \code{x} 
#'   should be scaled to have unit variance before the analysis takes place. The
#'   default is \code{FALSE} for consistency with \code{cancor}. Alternatively, 
#'   a vector of length equal the number of columns of \code{x} can be supplied.
#'   The value is passed to \code{\link{scale}}.
#' @param yscale analogous to \code{xscale}
#' @param npairs the number of pairs of canonical variables to be computed. With
#'   the default setting, pairs of variables are computed until either \code{x} or
#'   \code{y} is fully deflated.
#' @param xpredict the regression function to predict the canonical variable for
#'   \code{x}, given \code{y}. The formal arguments are the design matrix 
#'   \code{y}, the regression target \code{xv} as the current canonical variable
#'   for \code{x}, and \code{pp} as a counter of the current pair of canonical 
#'   variables (e.g. for enforcing different constraints for different canonical
#'   vectors). See the examples for an illustration.
#' @param ypredict analogous to \code{xpredict}
#' @param cor_tol a threshold indicating the magnitude below which canonical 
#'   variables should be omitted. Variables are omitted if their explained 
#'   correlations are less than or equal to \code{cor_tol} times the correlation
#'   of the first pair of canonical variables. With the default \code{NULL} 
#'   setting, no variables are omitted.
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
#' @return \code{nscancor} returns a list with the following elements: 
#'   \item{cor}{the additional correlation explained by each pair of canonical 
#'   variables, see \code{\link{acor}}.} \item{xcoef}{the matrix containing the 
#'   canonical vectors related to \code{x} as its columns} 
#'   \item{ycoef}{analogous to \code{xcoef}} \item{xcenter}{if \code{xcenter} is
#'   \code{TRUE} the centering vector, else the zero vector (in accordance with 
#'   \code{cancor})} \item{ycenter}{analogous to \code{xcenter}} 
#'   \item{xscale}{if \code{xscale} is \code{TRUE} the scaling vector, else 
#'   FALSE } \item{yscale}{analogous to \code{xscale}}
#'   
#' @note Deflating the data matrices accumulates numerical errors over 
#'   successive canonical vectors.
#'   
#' @references Sigg, C. and Fischer, B. and Ommer, B. and Roth, V. and Buhmann, 
#'   J. (2007) Nonnegative CCA for Audiovisual Source Separation. In 
#'   \emph{Proceedings of the 2007 IEEE Workshop on Machine Learning for Signal 
#'   Processing} (pp. 253--258).
#' @references Mackey, L. (2009) Deflation Methods for Sparse PCA. In 
#'   \emph{Advances in Neural Information Processing Systems} (pp. 1017--1024).
#'   
#' @seealso  \code{\link{acor}}, \code{\link{cancor}}, \code{\link{scale}}
#'   
#' @example inst/nscancor_examples.R
nscancor <- function (x, y, xcenter = TRUE, ycenter = TRUE, 
                      xscale = FALSE, yscale = FALSE, npairs = min(dim(x), dim(y)),
                      xpredict, ypredict, 
                      cor_tol = NULL, nrestart = 10, iter_tol = 1e-3, iter_max = 30,
                      verbosity = 0) {
    
    n <- nrow(x)
    dx <- ncol(x)
    dy <- ncol(y)
    
    X <- as.matrix(x)
    X <- scale(X, center = xcenter, scale = xscale)
    xcen <- attr(X, "scaled:center")
    xsc <- attr(X, "scaled:scale")
    Y <- as.matrix(y)
    Y <- scale(Y, center = ycenter, scale = yscale)
    ycen <- attr(Y, "scaled:center")
    ysc <- attr(Y, "scaled:scale")
    if(any(xsc == 0) || any(ysc == 0))
        stop("cannot rescale a constant/zero column to unit variance")
    
    corr <- rep(0, npairs)  # additional explained correlation
    W <- matrix(0, dx, npairs)  # canonical vectors for X
    V <- matrix(0, dy, npairs)  # canonical vectors for Y
    Qx <- matrix(0, dx, npairs)  # orthonormal basis spanned by the canonical variables X%*%W
    Qy <- matrix(0, dy, npairs)  # orthonormal basis spanned by the canonical variables Y%*%V
    Xp <- X  # X projected to the orthocomplement space spanned by Qx
    Yp <- Y  # Y projected to the orthocomplement space spanned by Qy
    
    for (pp in seq(npairs)) {
        obj_opt <- -Inf
        for (rr in seq(nrestart)) {
            
            res <- emcca(X, Xp, Y, Yp, xpredict, ypredict, pp, iter_tol, 
                         iter_max, verbosity)               
            
            # keep solution with maximum objective
            if (res$obj > obj_opt) {
                obj_opt <- res$obj
                w_opt <- res$w
                v_opt <- res$v
            }
            if (verbosity > 0) {
                print(paste("variable pair ", pp, ": ",
                            "maximum objective is ", format(obj_opt, digits = 4),
                            " at random restart ", rr-1, sep = ""))
            }
        }
        w <- w_opt  
        v <- v_opt
        W[ ,pp] <- w
        V[ ,pp] <- v
        corr[pp] <- cor(Xp%*%w, Yp%*%v)
        
        # update Qx and Qy
        XtXw <- t(X)%*%(X%*%w)
        YtYv <- t(Y)%*%(Y%*%v)
        if (pp > 1) {
            qx <- XtXw - Qx[ , 1:(pp-1)]%*%(t(Qx[ , 1:(pp-1)])%*%XtXw) 
            qy <- YtYv - Qy[ , 1:(pp-1)]%*%(t(Qy[ , 1:(pp-1)])%*%YtYv) 
        } else {
            qx <- XtXw
            qy <- YtYv
        }
        qx <- qx/normv(qx)
        qy <- qy/normv(qy)
        Qx[ , pp] <- qx
        Qy[ , pp] <- qy
        
        # deflate data matrices
        Xp <- Xp - Xp%*%qx%*%t(qx)  
        Yp <- Yp - Yp%*%qy%*%t(qy)  
        
        # current additionally explained correlation is below threshold cor_tol
        if (!is.null(cor_tol) && corr[pp] < cor_tol*corr[1]) {
            corr <- corr[1:(pp-1)]
            W <- W[ , 1:(pp-1), drop=FALSE]
            V <- V[ , 1:(pp-1), drop=FALSE]
            break
        }    
        # at least one data matrix is fully deflated
        else if (pp < npairs && (all(abs(Xp) < 1e-14) || all(abs(Yp) < 1e-14))) { 
            if (verbosity > 0) {
                print("at least one data matrix is fully deflated, less than 'npairs' pairs of variables could be computed")            
            }
            corr <- corr[1:pp]
            W <- W[ , 1:pp, drop=FALSE]
            V <- V[ , 1:pp, drop=FALSE]
            break
        }
    }
    
    rownames(W) <- colnames(X)
    rownames(V) <- colnames(Y)
    return(list(cor = corr, xcoef = W, ycoef = V,
                xcenter = if(is.null(xcen)) rep.int(0, dx) else xcen,  # return value follows cancor interface
                xscale = if(is.null(xsc)) FALSE else xsc,
                ycenter = if(is.null(ycen)) rep.int(0, dy) else ycen,
                yscale = if(is.null(ysc)) FALSE else ysc))
}

emcca <- function(X, Xp, Y, Yp, xpredict, ypredict, pp, iter_tol, iter_max, 
                  verbosity) {
    
    n <- nrow(X)
    dx <- ncol(X)
    dy <- ncol(Y)
    
    # initialize w and v as random unit vectors, non-negative if the prediction
    # functions return non-negative vectors
    w <- rnorm(dx);
    v <- rnorm(dy);
    if (all(ypredict(Xp, Yp[ , 1], pp) >= 0)) {
        w <- abs(w)
    }
    if (all(xpredict(Yp, Xp[ , 1], pp) >= 0)) {
        v <- abs(v)
    }
    w <- w/normv(Xp%*%w)     
    v <- v/normv(Yp%*%v)     
    
    obj_old <- -Inf
    ii <- 0
    while(ii < iter_max) { 
        obj <- t(Xp%*%w)%*%(Yp%*%v)
        if (abs(obj - obj_old)/obj < iter_tol) {
            break
        }
        obj_old <- obj
        
        w <- ypredict(Xp, Yp%*%v, pp)
        if (all(w == 0))
            stop("w collapsed to the zero vector, try relaxing the constraints")
        w <- w/normv(Xp%*%w)
        
        v <- xpredict(Yp, Xp%*%w, pp)
        if (all(v == 0))
            stop("v collapsed to the zero vector, try relaxing the constraints")
        v <- v/normv(Yp%*%v)
        
        ii <- ii + 1
    }
    if ((verbosity > 0 ) && (ii == iter_max))
        print("maximum number of iterations reached before convergence")
    
    w <- w/normv(X%*%w)     
    v <- v/normv(Y%*%v) 
    return(list(obj=obj, w=w, v=v))
}