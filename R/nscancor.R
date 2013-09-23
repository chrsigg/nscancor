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

#' Non-Negative Sparse CCA
#' 
#' Performs a constrained canonical correlation analysis,
#' where non-negativity, sparsity or other constraints are enforced on the 
#' projection vectors. The result of the analysis is returned as a list with 
#' the same elements as the list returned by \code{cancor}.
#' 
#' \code{nscancor} computes the projection vectors using iterated regression
#' steps, where the constraints suitable for each domain are enforced by 
#' choosing the appropriate regression method. See Sigg et al. (2007) for an
#' early application of the principle (not yet including generalized deflation). 
#' 
#' Because constrained projection vectors no longer correspond to true eigenvectors 
#' of the cross-covariance matrix and are usually not pairwise conjugate, special
#' attention needs to be paid when computing more than a single component. The
#' algorithm implements a generalized deflation (GD) scheme which builds on
#' GD for PCA as proposed by Mackey (2009). Given a basis of the space spanned by 
#' the canonical variables, the
#' correlation of the component is maximized after projecting the current 
#' projection vector to the
#' ortho-complement space of the basis. This procedure maximizes the
#' additional correlation not explained by previous components, and is
#' identical to standard CCA if no additional constraints
#' are enforced.
#' 
#' See the references for further details. 
#' 
#' @export nscancor
#' @param x a numeric matrix which provides the data from the first domain
#' @param y a numeric matrix which provides the data from the second domain
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
#' @param ncomp the number of canonical components to be computed. With the 
#'   default setting, components are 
#'   computed until either \code{x} or \code{y} is fully deflated. 
#' @param xpredict the regression function to predict the canonical variable of 
#'   \code{x} given \code{y}. The formal arguments of the function are
#'   the design matrix \code{y}, the current canonical variable of \code{x} as
#'   the target and the current component \code{cc}, e.g. for specifying 
#'   different constraints for each component.
#' @param ypredict see \code{xpredict}
#' @param cor_tol a threshold indicating the magnitude below which components
#'   should be omitted. Components are omitted if their explained
#'   correlations are less than or equal to \code{cor_tol} times the
#'   correlation of the first component.
#'   With the default \code{NULL} setting, no components
#'   are omitted.  With \code{cor_tol = 0} or \code{cor_tol = sqrt(.Machine$double.eps)}, 
#'   essentially constant components are omitted.
#' @param nrestart the number of random restarts for computing
#'   the projection pairs via iterated regression. The solution 
#'   achieving maximum explained correlation over all random restarts is kept. A 
#'   value greater than one can help to avoid bad local maxima.
#' @param iter_tol If the relative change
#'   of the objective is less than \code{iter_tol} between iterations, 
#'   the procedure is asssumed to have converged to a local optimum.
#' @param iter_max the maximum number of iterations to be performed. The
#'   procedure is terminated if either the \code{iter_tol} or the \code{iter_max}
#'   criterion is satisfied.
#' @param verbosity an integer specifying the verbosity level. Greater values
#'   result in more output, the default is to be quiet.
#' 
#' @return \code{nscancor} returns a list containing the following elements:
#' \item{cor}{the additional correlation explained by each component,
#'   see \code{\link{acor}}.}
#' \item{xcoef}{the matrix containing the projection vectors related to 
#'   \code{x} as its columns}
#' \item{ycoef}{the matrix containing the projection vectors related to 
#'   \code{y} as its columns}
#' \item{xcenter}{if \code{xcenter} is \code{TRUE} the centering vector, else
#'   the zero vector }
#' \item{xscale}{if \code{xscale} is \code{TRUE} the scaling vector, else
#'   FALSE }
#' \item{ycenter}{see \code{xcenter}}
#' \item{yscale}{see \code{xscale}}
#' 
#' @note Deflating the data matrices accumulates numerical errors over successive
#' components.
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
                      xscale = FALSE, yscale = FALSE, ncomp = min(dim(x), dim(y)),
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
    
    corr <- rep(0, ncomp)  # additional explained correlation
    W <- matrix(0, dx, ncomp)  # canonical axes for X
    V <- matrix(0, dy, ncomp)  # canonical axes for Y
    Qx <- matrix(0, dx, ncomp)  # orthonormal basis spanned by the canonical variables X%*%W
    Qy <- matrix(0, dy, ncomp)  # orthonormal basis spanned by the canonical variables Y%*%V
    Xp <- X  # X projected to the orthocomplement space spanned by Qx
    Yp <- Y  # Y projected to the orthocomplement space spanned by Qy
    
    for (cc in seq(ncomp)) {
        obj_opt <- -Inf
        for (rr in seq(nrestart)) {
            
            res <- emcca(X, Xp, Y, Yp, xpredict, ypredict, cc, iter_tol, 
                         iter_max, verbosity)               
            
            # keep solution with maximum objective
            if (res$obj > obj_opt) {
                obj_opt <- res$obj
                w_opt <- res$w
                v_opt <- res$v
            }
            if (verbosity > 0) {
                print(paste("component ", cc, ": ",
                            "maximum objective is ", format(obj_opt, digits = 4),
                            " at random restart ", rr-1, sep = ""))
            }
        }
        w <- w_opt  
        v <- v_opt
        W[ ,cc] <- w
        V[ ,cc] <- v
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
        
        # current additionally explained correlation is below threshold cor_tol
        if (!is.null(cor_tol) && corr[cc] < cor_tol*corr[1]) {
            corr <- corr[1:(cc-1)]
            W <- W[ , 1:(cc-1), drop=FALSE]
            V <- V[ , 1:(cc-1), drop=FALSE]
            break
        }    
        # at least one data matrix is fully deflated
        else if (cc < ncomp && (all(abs(Xp) < 1e-14) || all(abs(Yp) < 1e-14))) { 
            if (verbosity > 0) {
                print("at least one data matrix is fully deflated, less than 'ncomp' components could be computed")            
            }
            corr <- corr[1:cc]
            W <- W[ , 1:cc, drop=FALSE]
            V <- V[ , 1:cc, drop=FALSE]
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

emcca <- function(X, Xp, Y, Yp, xpredict, ypredict, cc, iter_tol, iter_max, 
                  verbosity) {
    
    n <- nrow(X)
    dx <- ncol(X)
    dy <- ncol(Y)
    
    # initialize w and v as random unit vectors, non-negative if the prediction
    # functions return non-negative vectors
    w <- rnorm(dx);
    v <- rnorm(dy);
    if (all(ypredict(Xp, Yp[ , 1], cc) >= 0)) {
        w <- abs(w)
    }
    if (all(xpredict(Yp, Xp[ , 1], cc) >= 0)) {
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
        
        w <- ypredict(Xp, Yp%*%v, cc)
        if (all(w == 0))
            stop("w collapsed to the zero vector, try relaxing the constraints")
        w <- w/normv(Xp%*%w)
        
        v <- xpredict(Yp, Xp%*%w, cc)
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