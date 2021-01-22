


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
# stand_dissim: standardizes dissimilarity matrix
#
#
#

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


calc_stress <- function(mX, mDis,mEucD){
  
  mXtX <- t(mX) %*% mX

  eta2 <- nrow(mX) * sum(diag((mXtX)))
  
  mBx <- calc_Bx(mDis, mEucD)
  
  stress <- eta2 - 2* mBx
  
  return(stress)
}

### 
# calc_Bx: calculates matrix Bx
#
#
#


calc_Bx <- function(mDis, mEucD){
  
  mF <- mDis/mEucD
  diag(mF) <- 0
  
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
  k = 0 
  
  # create initial random matrix with dimension nx2 
  mZ <-matrix(rexp(n*2, rate=.1), ncol=2)

  # get initial euclidean distances
  mEucDZ <- euc_dist_M(mZ)
  
  # define Bx with dissimilarities and euclidean matrix
  mBz <- calc_Bx(mDis, mEucDZ)
  
  # intiialize variables to store previous stress and current stress
  stress_prev <- 0
  stress_k <- sum(calc_stress(mZ, mDis,mEucDZ), na.rm=T)
  
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

mDis_stand <- stand_dissim(as.matrix(basket))
nrow(mDis_stand)

SMACOF(mDis_stand)

