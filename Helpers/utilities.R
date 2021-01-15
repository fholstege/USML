#' Binary Search for Penalization Parameter
#' 
#' Performs binary search to find a value of the shrinkage parameter
#' lambda which obtains a required L1-norm value after soft-thresholding and
#' normalization to length one. The value of lambda found is the smallest value 
#' that attains the required constraint value.
#' 
#' @details Note that in the \code{\link{PMA}} package (version 1.2), the binary search will find the largest 
#' value of lambda that imposes the required limit on the L1-norm after soft-thresholding in case
#' this norm has regions of lambda for which it is constant. This leads to sparser solutions
#' than what is necessary.
#' 
#' @param x Vector to be soft-thresholded and normalized
#' @param c Desired value of the L1-norm after soft-thresholding and normalization
#' @param maxit Maximum number of iterations for the binary search
#' 
#' @author Pieter C. Schoonees
#' 
#' @examples 
#' 
#' ## Random vector to soft-threshold and normalize
#' set.seed(1)
#' (z <- rnorm(25))
#' 
#' ## Binary search for L1-norm of 1 (returning lambda)
#' (lambda <- binary_search(z, 1))
#' 
#' ## Graphical illustration: Progress of estimate
#' pts <- seq(from = 0, to = max(abs(z)), length.out = 100)
#' plot(pts, sapply(pts, function(a) sum(abs(soft_l2norm(z, a)))),
#'      ylab = "L1-norm", xlab = expression(lambda))
#' abline(h = 1, lty = 2)
#' abline(v = binary_search(z, 1, maxit = 1), col = 2, lty = 2)
#' abline(v = binary_search(z, 1, maxit = 2), col = 2, lty = 2)
#' abline(v = binary_search(z, 1, maxit = 3), col = 2, lty = 2)
#' abline(v = binary_search(z, 1, maxit = 4), col = 2, lty = 2)
#' abline(v = binary_search(z, 1, maxit = 15), col = 2, lty = 2)
#' 
#' ## More tests: Check the L2-norm is the required one
#' sum(abs(soft_l2norm(z, lambda)))
#' sum(abs(soft_l2norm(z, binary_search(z, 2))))
#' sum(abs(soft_l2norm(z, binary_search(z, 3))))
#' 
#' @export
binary_search <- function(x, c, maxit = 100) {
  
  ## Check if x already meets the condition, or whether c is zero, and return 0 in either case
  x_l1 <- sum(abs(x))
  if (x_l1 <= c || c == 0) {
    return(0)
  }
  
  ## Initialize search boundaries and iteration counter
  lo <- 0
  hi <- max(abs(x))
  i <- 0
  
  ## Iterate
  while(i < maxit) {
    
    ## Increase iteration count
    i <- i + 1
    
    ## Get mean of current bounds, and evaluate L1 norm there
    m <- lo + (hi - lo)/2
    m_l1 <- sum(abs(soft_l2norm(x, m)))
    
    ## Update either lower or upper threshold
    if (m_l1 <=  c) { ## Note: must be <= otherwise lambda will be too large (as in PMA)
      hi <- m
    } else {
      lo <- m
    }
    
    ## Break if lo and hi have converged
    if (abs(lo - hi) <= sqrt(.Machine$double.eps)) {
      break
    }
  }
  
  return(lo + (hi - lo)/2)
}

#' Soft-thresholding Function
#' 
#' Simple vectorized soft-thresholding utility function.
#' 
#' @param x Vector to be soft-thresholded
#' @param lambda Threshold value
#' 
#' @examples 
#' set.seed(1)
#' (z <- rnorm(10))
#' soft_thresh(z, 0.5)
#' 
#' y <- seq(from = -2, to = 2, by = 0.1)
#' plot(y, soft_thresh(y, 0.4), type = "l")
#' lines(y, y, lty = 3)
#' 
soft_thresh <- function(x, lambda) {
  sign(x) * (pmax(abs(x) - lambda, 0))
}
#' Soft-threshold and normalize
#' 
#' Soft-threshold, and then normalize the result to L2-norm of 1. Care is taken with
#' vectors of L2-norm empirically equal to zero.
#' 
#' @param x Vector to be soft-thresholded
#' @param lambda Threshold value
#' 
#' @examples 
#' ## Check that L2-norm is correct after normalization
#' set.seed(1)
#' (z <- rnorm(10))
#' soft_l2norm(z, 0.5)
#' sqrt(sum(soft_l2norm(z, 0.5)^2))
#' 
#' (z <- rnorm(100))
#' soft_l2norm(z, 0.75)
#' sqrt(sum(soft_l2norm(z, 0.5)^2))
#' 
soft_l2norm <- function(x, lambda) {
  x <- soft_thresh(x = x, lambda = lambda)
  x_l2 <- sqrt(sum(x^2))
  if (abs(x_l2) > .Machine$double.eps^0.5) {
    x <- x / x_l2
  }
  x
}
