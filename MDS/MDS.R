
# Packes required for subsequent analysis. P_load ensures these will be installed and loaded. 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(smacof) 




# load the data
load("Data/basket.Rdata")




###
# mEuc_dist_M: calculates euclidean distance of an nxp matrix for each row with all other rows
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
#LowerT_sum: takes sum of lower triangular of a matrix
# 
#
#
LowerT_sum <- function(mX){sum(lower.tri(mX),  na.rm=T)}


### 
# stand_dissim: standardizes dissimilarity matrix
#   
# Arguments: 
#       mDis: nxn dissimilarity matrix
#
# Output: 

stand_dissim <- function(mDis){
  
  diag_mDis <- diag(mDis)
  
  mDis_stand <- 1 - (mDis/diag_mDis)
  
  return(mDis_stand)
}

### 
# calc_stress: calculates stress
#
#
#


calc_stress <- function(mX, mDis,mEucD, normalize = TRUE){
  
  # define X'X
  mXtX <- t(mX) %*% mX

  eta2 <- nrow(mX) * sum(diag((mXtX)))
  
  mBx <- calc_Bx(mDis, mEucD)
  
  rho <- sum(diag(t(mX) %*% mBx %*% mX))
  
  stress <- eta2 - 2* rho
  
  if(normalize){
    
    stress <- stress / mEucD^2
    
  }
    
  return(stress)
}

### 
# calc_Bx: calculates matrix Bx
#
#
#


calc_Bx <- function(mDis, mEucD){
  
  # create matrix F
  mF <- mDis/mEucD
  diag(mF) <- 0
  
  # create Bx matrix
  mBx <- diag(colSums(mF)) - mF
  
  return(mBx)
}


### 
# SMACOF: implements SMACOF algorithm
#
#
#


SMACOF <- function(mDis, mX = "random", eps = 0.0000005){
  
  # set n and k 
  n = nrow(mDis)
  p = 2
  k = 0 
  
  if(mX == "random"){
    # create initial random matrix with dimension nx2 
    mZ <-matrix(rexp(n*p, rate=.1), ncol=p)
  }else{
    
    mZ <- mX
  }

  # get initial euclidean distances
  mEucDZ <- euc_dist_M(mZ)
  
  # define Bx with dissimilarities and euclidean matrix
  mBz <- calc_Bx(mDis, mEucDZ)
  
  # intiialize variables to store previous stress and current stress
  stress_prev <- 0
  stress_k <- LowerT_sum(calc_stress(mZ, mDis,mEucDZ))
  
  # initialize while loop 
  while(k == 0 | stress_prev - stress_k > eps){
    
    # to continue iterations
    k = k + 1
    
    # update mX_k 
    mX_k <- (1/n) * mBz %*% mZ

    # update Z and the euclidean distance of Z
    mZ <- mX_k
    mEucDZ <- euc_dist_M(mZ)
    
    # set previous stress to current, and update current stress
    stress_prev <- stress_k
    stress_k <- sum(calc_stress(mZ,mDis, mEucDZ), na.rm = T)
    
    # check if constantly decreasing: 
    print(paste0("Iteration : ", k," Current stress: ", stress_k))
    
    
  }
}


# Create dissimilarity matrix
df <- 1-basket/diag(as.matrix(basket))

# labels 
labels <- colnames(basket)
mDis_stand <- stand_dissim(as.matrix(basket))
nrow(mDis_stand)

# obtain initial configeration using classical torgerson
mX_0 <- mds(df, ndim=2, itmax=1, eps=1-06)

# package results
pack.result <- mds(df, itmax = 10000000, init='torgerson', eps=1e-06)
pack.result$stress
pack.result$conf
pack.result$niter

mX_0$conf

SMACOF(mDis_stand, mX = "random")

