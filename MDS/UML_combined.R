#################
# Authors: Floris Holstege, Yuchou Peng, Chao Liang, Eva Mynott, Gabrielle Mignoli
# Purpose: Code & Text for Unsupervised Machine Learning Assignment 2 - applying MDS to a basket of co-purchased goods
# Sections: 
#           A) Load packages and data
#           B) Define our functions for SMACOF + dissimilarity matrix + visualization
#           C) Implementation SMACOF algorithm
#           D) Explicit answers to questions in the exercise 
#################


################################################################################
# A) Load packages
################################################################################

# empty global environment
rm(list=ls())

# Packes required for subsequent analysis. P_load ensures these will be installed and loaded. 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(smacof,
               ggplot2,
               rlist,
               fBasics) 

# load the data
load("Data/basket.Rdata")

# set seed
set.seed(123)


################################################################################
# B) Define our functions for SMACOF + dissimilarity matrix + visualization
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
# UpperT_sum: takes sum of lower triangular of a matrix
# 
# Arguments: 
#     mX: matrix 
#
# Output:
#     constant, sum of lower triangular 
# 
UpperT_sum <- function(mX){sum(upper.tri(mX), diag = FALSE)}


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

### calc_stress: Calculate normalized stress value
#
# Arguments; 
#     mDis: nxn dissimilarity matrix
#     X: nxp input matrix
#     tX.X: X'X, pxp matrix
#     B: matrix B, calculated with calc_mB
# 
#   Output:
#      stress: float, normalized stress
#     
calc_stress <- function(mDis, mX, tX.X, mB, normalized=TRUE) {
  
  # n of variables in 
  n <- nrow(mX)
  
  # constant eta, eta squared and rho make up the raw stress
  Eta <- (sum(lower.tri(mDis, diag = FALSE) * (mDis ^ 2 )))
  Eta <- UpperT_sum(mDis^2)
  Eta2 <- n*sum(diag(tX.X))
  rho <- sum(diag(t(X) %*% B %*% X))
  
  # calculate raw stress
  raw_stress <- Eta + Eta2 - 2*rho

  # normalize the stress and return 
  stress <- raw_stress/Eta

  return(stress)
}


### create_coordinate_plot: creates plot of coordinates from MDS result
#
# Arguments; 
#     result_MDS: result object from mds() function or SMACOF, with conf attribute
#
#   Output:
#      coordinate_plot: ggplot object with MDS coordinates
#     
create_coordinate_plot <- function(result_MDS){
  
  # plot of mds results
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
# C) Smacof algorithm
################################################################################


### SMACOF: implements SMACOF algorithm to a dissimilarity matrix
#
# Arguments; 
#     mDis: nxn dissimilarity matrix
#    
#   Output:
#      output: list, with final configuration, n of iterations, and stress of final configuration
#  

SMACOF <- function(mDis, config = NULL, eps = 1e-06) {
  
  # Number of variables
  n <- nrow(mDis)
  p <- 2
  
  # Use random or initial Torgerson configuration
  if (is.null(config) == T) {
    mX<- matrix(data = rnorm(n*p,), nrow = n, ncol = p)
  } else {
    mX<- config
  }

  # Set Z to be equal to the initial configuration
  mZ <- list(mX)
  
  # Get initial euclidean distance
  mDZ <- euc_dist_M(mX)

  # Get initial B
  mBZ <- calc_mB(mDis, mDZ)
  
  # Pre-calculate stress value
  stress <- c(0)
  stress <- c(stress,calc_stress(mDis, mX, t(mX)%*%mX, mBZ))
  
  # Initialize counter
  k <- 1
  
  # Start loop
  while (k == 1 | stress[k]-stress[k+1] > eps){
    k <- k+1
    
    # Update X
    mX<- (1/n) * mBZ %*% mZ[[k-1]]
    mZ <- list.append(mZ,mX)
    
    # Obtain new distances
    mDZ <- euc_dist_M(mZ[[k]])
    
    # Calculate updated B
    mBZ <- calc_mB(mDis, mDZ)
    
    # Compute updated stress
    stress <- c(stress, calc_stress(mDis, mZ[[k]], t(mZ[[k]]) %*% mZ[[k]], mBZ))
    
    # Check whether stress goes down
    if (stress[k+1]-stress[k] > 0) {
      cat('\n\nSomething goes wrong! Stress is increasing by ', stress[k+1]-stress[k])
    }
    
  } # End loop
  
  cat('\n\nIterations:', k) # Total number of iterations
  
  # Output a final stress and configuration
  cat('\n\nFinal stress:', stress[k])
  cat('\n\nFinal configuration:\n')
  print(data.frame(mZ[[k]]))
  
  # return list with results
  output <- list(conf = data.frame(mZ[[k]]), stress = stress[k], niter = k)
  
  return(output)
}


################################################################################
# D) Explicit answers to questions in the exercise 
################################################################################

### A) 

# create dissimilarity matrix and ordinal dissimilarity matrix
labels <- colnames(basket)
mDis <- calc_mDis(as.matrix(basket))
mDis_Ord <- calc_mDis_Ord(mDis)

### B) The function for the euclidean distance matrix is euc_dist_M, found in line 83-100

### C) see the SMACOF() function for our implemenation, section C of the code. We implement it below a

# obtain initial configeration using torgerson algorithm
initial <- mds(mDis, init='torgerson', itmax=1)
initial_conf <- initial$conf

# own implementation, with initial configuration
own.result <- SMACOF(mDis, config = initial_conf, eps=1e-6) 
own.result$stress
own.result$conf
own.result$niter

### D)  There is a difference of 0.083 between our implemenation and the package in terms of raw stress
### Why this difference? 1/n = 1/12 = 0.083 - so likely we are somewhere not scaling the stress accurately, which influences our configuration

# package results
pack.result <- mds(mDis, itmax = 1000, init=initial_conf, eps=1e-06)
pack.result$stress
pack.result$conf
pack.result$niter

# plots of both implementations
create_coordinate_plot(own.result)
create_coordinate_plot(pack.result)

### E) The ordinal solution finds a much more 'clustered' solution - the products are often much closer or further away from each other, 
# instead of an roughly equal distance in the previous implemenation. This is likely because of our first dissimilarity matrix has many similar values (range of 0.6-0.8)
# but in the ordinal dissimilarity matrix, the differences become bigger, since its only considered with the relative rank

# check the ordinal result
pack.result_ord <- mds(mDis, itmax = 1000, init=initial_conf, eps=1e-06,  type = 'ordinal')
create_coordinate_plot(pack.result_ord)

