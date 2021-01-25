#################
# Authors: Floris Holstege, Yuchou Peng, Chao Liang, Eva Mynott, Gabrielle Mignoli
# Purpose: Code for Unsupervised Machine Learning Assignment 2 - applying MDS to a basket of co-purchased goods
# Sections: 
#           A) Load packages and pre-process the data
#           B) Define our functions for MDS
#           C) Implementation SMACOF algorithm
#           D) Compare our implementation with package
#################



################################################################################
# A) Load packages and pre-process the data
################################################################################

# empty global environment
rm(list=ls())

# Packes required for subsequent analysis. P_load ensures these will be installed and loaded. 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(smacof,
               ggplot2,
               rlist) 

# load the data
load("Data/basket.Rdata")

# set seed
set.seed(123)


################################################################################
# B) Define our functions for MDS
################################################################################


### calc_mDis: Create dissimilarity matrix
# 
# Arguments:
#     mX: nxn matrix 
# Output:
#     mDis: nxn dissimilarity matrix
#

calc_mDis <- function(mX){
  
  # scale similarities (covariances -> correlations)
  mX_transformed <- solve(diag(sqrt(diag(mX)))) %*% mX %*% solve(diag(sqrt(diag(mX))))
  
  # obtain dissimilarities from scaled similarities
  mDis <- 1-mX_transformed 
  
  # make sure diagonal is exactly 0 (not rounded)
  diag(mDis) <- 0
  
  return(mDis)
}


### calc_mDis_Ord: Turn dissimilarities to ordinal dissimilarities 
# 
# Arguments:
#     mDis: nxn dissimilarity matrix 
# Output:
#     mDis_ord: nxn dissimilarity matrix, iwith ordinal rank
#
calc_mDis_Ord <- function(mDis) {
  
  # initialize new matrix of ranks
  mOrd <- matrix(0L, nrow=nrow(mDis), ncol=ncol(mDis))
  
  # replace upper and lower matrix with rank values
  mOrd[lower.tri(mOrd, diag = FALSE)] <- rank(mDis[lower.tri(mDis, diag = FALSE)])
  mOrd[upper.tri(mOrd, diag = FALSE)] <- rank(mDis[upper.tri(mDis, diag = FALSE)])
  
  # ensure values are integers
  mDis_ord <- mapply(mOrd, FUN=as.integer)
  mDis_ord <- matrix(data=mOrd, ncol=ncol(mDis), nrow=nrow(mDis))
  
  return(mDis_ord)
}



###
# euc_dist_M: calculates euclidean distance of an nxp matrix for each row with all other rows
#
# Arguments:
#       mX: nxp matrix object
# 
# Output: 
#       result: nxn matrix of euclidean distances

euc_dist_M <- function(mX){
  
  # define M as XX', and C as diag(M)
  M <- mX %*% t(mX)
  C <- diag(M)
  
  # calculate euclidean distance matrix
  result <- sqrt(outer(C, C, '+') - 2*M)
  
  return(result)
}

###
# LowerT_sum: takes sum of lower triangular of a matrix
# 
# Arguments: 
#     mX: matrix 
#
# Output:
#     constant, sum of lower triangular 
# 
LowerT_sum <- function(mX){sum(lower.tri(mX), diag = FALSE)}


### 
# calc_mB: calculates matrix B
#
# Arguments: 
#     mDis: nxn dissimilarity matrix 
#     mEucD: nxn euclidean distance matrix
# 
# Output:
#     mBx: nxn matrix B, used for calculating stress
# 
calc_mB <- function(mDis, mEucD){
  
  # create matrix F
  mF <- mDis/mEucD
  diag(mF) <- 0
  
  # create Bx matrix
  mBx <- diag(rowSums(mF)) - mF
  
  return(mBx)
}



### calc_stress: Calculate (normalized) stress value
#
# Arguments; 
#     mDis: nxn dissimilarity matrix
#     X: 
#     tX.X:
#     B
#     normalized: boolean, if true then normalize stress
# 
#   Output:
#      stress: variable with raw or normalized stress
#     


calc_stress <- function(mDis, X, tX.X, B, normalized=TRUE) {
  
  # n of variables in 
  n <- nrow(X)
  
  # constant eta, eta squared and rho make up the raw stress
  Eta <- LowerT_sum(mDis^2)
  Eta2 <- n*sum(diag(tX.X))
  rho <- sum(diag(t(X) %*% B %*% X))
  
  # calculate raw stress
  raw.stress <- Eta + Eta2 - 2*rho

  if(normalized){
    
    # normalize the raw stress 
    norm.stress <- raw.stress/Eta
    stress <- norm.stress
    
  }else{
    stress <- raw.stress
  }
  return(stress)
}


################################################################################
# C) Smacof algorithm
################################################################################

SMACOF <- function(mDis, config = NULL, eps = 1e-06) {
  
  # Number of variables
  n <- nrow(mDis)
  
  # Use random or initial Torgerson configuration
  if (is.null(config) == T) {
    X <- matrix(data = rnorm(n*2,), nrow = n, ncol = 2)
  } else {
    X <- config
  }

  # Set Z to be equal to the initial configuration
  Z <- list(X)
  
  # Get initial euclidean distance
  DZ <- euc_dist_M(X)

  # Get initial B
  BZ <- calc_mB(mDis, DZ)
  
  # Pre-calculate stress value
  stress <- c(0)
  stress <- c(stress,calc_stress(mDis, X, t(X)%*%X, BZ))
  
  # Initialize counter
  k <- 1
  
  # Start loop
  while (k == 1 | stress[k]-stress[k+1] > eps){
    k <- k+1
    
    # Update X
    X <- (1/n) * BZ %*% Z[[k-1]]
    Z <- list.append(Z,X)
    
    # Obtain new distances
    DZ <- euc_dist_M(Z[[k]])
    
    # Calculate updated B
    BZ <- calc_mB(mDis, DZ)
    
    # Compute updated stress
    stress <- c(stress, calc_stress(mDis, Z[[k]], t(Z[[k]])%*%Z[[k]], BZ))
    
    # Check whether stress goes down
    if (stress[k+1]-stress[k] > 0) {
      cat('\n\nSomething goes wrong! Stress is increasing by ', stress[k+1]-stress[k])
    }
    
  } # End loop
  
  cat('\n\nIterations:', k) # Total number of iterations
  
  # Output a final stress and configuration
  cat('\n\nFinal stress:', stress[k])
  cat('\n\nFinal configuration:\n')
  print(data.frame(Z[[k]]))
  
  # return list with results
  output <- list(conf = data.frame(Z[[k]]), stress = stress[k], niter = k)
  
  return(output)
}


create_coordinate_plot <- function(result_MDS){
  
  
  # plot of package results
  dfResult <- data.frame(result_MDS$conf)
  rownames(dfResult) <- labels
  coordinate_plot <- ggplot(data = dfResult, aes(x = D1, y = D2))+
    geom_point(size = 1.5) + ggtitle('MDS co-purchases') + 
    ylab("") + xlab("") + theme_bw() + 
    coord_cartesian(ylim = c(-1, 1), xlim = c(-1, 1)) + 
    geom_text(aes(label=labels),hjust=0.3, vjust=1.3, size=4)
  
  return(coordinate_plot)
  
}


################################################################################
# D) Compare our implementation with package
################################################################################


labels <- colnames(basket)
mDis <- calc_mDis(as.matrix(basket))
mDis_Ord <- calc_mDis_Ord(mDis)

# obtain initial configeration using torgerson algorithm
initial <- mds(mDis, init='torgerson', itmax=1)
initial_conf <- initial$conf


# own implementation, with initial configuration
own.result <- SMACOF(mDis, config = initial_conf, eps=1e-6) 
own.result$stress
own.result$conf
own.result$niter


# package results
pack.result <- mds(mDis, itmax = 1000, init=initial_conf, eps=1e-06)
pack.result$stress
pack.result$conf
pack.result$niter

# plots of both implementations
create_coordinate_plot(own.result)
create_coordinate_plot(pack.result)


# check the ordinal result
pack.result_ord <- mds(mDis, itmax = 1000, init=initial_conf, eps=1e-06,  type = 'ordinal')
create_coordinate_plot(pack.result_ord)

