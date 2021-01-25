# packages ----------------------------------------------------------------

library(smacof)

# functions ---------------------------------------------------------------

fnEuclDist <- function(mX){
  
  #calculate cross product
  mTCross <- tcrossprod(mX)
  
  #calculate diagonals
  vC <- diag(mTCross)
  
  #calculate distance matrix
  mDist <- outer(vC, vC, FUN = "+") - 2 * mTCross
  
  # take square root
  mDist <- sqrt(mDist)
  return(mDist)
}

fnSMACOF<- function(mDiss, mX = NULL, dEps = 1e-6, verbose = TRUE){
  # get size
  iN <- nrow(mDiss)
  
  #fill mX with random numbers
  if (is.null(mX)){
    mX = matrix(runif(2*iN)*max(mDiss), ncol = 2, nrow = iN)
  }
  
  #set counter
  iK <- 0
  
  #initialize Z and compute distances, stress
  mZ <- mX
  mDist <- fnEuclDist(mZ)
  dSumDiss <- sum(mDiss^2)
  dStress <- 0.5*sum((mDist - mDiss)^2)/dSumDiss
  dTol <- dEps + 1
  
  while (iK == 0 || dTol > dEps){
    #update k
    iK <- iK + 1
    
    #find Bx
    mF <- ifelse(mDist > 0, mDiss/mDist, 0)
    mF1 <- diag(rowSums(mF))
    mBx <- mF1 - mF
    
    #update X, z + calculate distances
    mX = (1/iN)*(mBx %*% mZ)
    mZ = mX
    mDist <- fnEuclDist(mZ)
    
    #caclulate stress and find toll
    dStressOld <- dStress
    dStress <- 0.5*sum((mDist - mDiss)^2)/dSumDiss
    dTol <- dStressOld - dStress
    
    #update
    if(verbose){
      cat("iteration ", formatC(iK), " Obj value ", formatC(dStress, digits = 5),
          " improvement ", formatC(dTol), "\n")
    }
  }
  return(list(dStress = dStress, mX = mX))
}

#set seed and load data
set.seed(1234)
load("~/UML/basket.Rdata")

#set basket to matrix and calculate dissimilarities (eucl dist)
mBasket <- as.matrix(basket)
iN <- ncol(mBasket)
mDiss <- matrix(0, ncol = iN, nrow = iN)#fnEuclDist(mBasket)
for( i in 1:iN){
  for (j in 1:iN){
    mDiss[i,j] <- 0.5*(1 - mBasket[i,j]/mBasket[i,i]) + 0.5*(1- mBasket[i,j]/mBasket[j,j])
  }
}
#lResults <-} fnSMACOF(mDiss = mDiss)
#lBestResults <- lResults


#for (i in 1:1000){
#  print(i)
#  lResults <- fnSMACOF(mDiss = mDiss, verbose = FALSE)
#  if(lResults$dStress<= lBestResults$dStress){
#    lBestResults <- lResults
# }
#}

#plot(x = lBestResults$mX[,1], y = lBestResults$mX[,2], xlab = "Axis 1", ylab = "Axis 2")
#text(x =  lBestResults$mX[,1], y = lBestResults$mX[,2], colnames(basket), col = "red")


#do 1 iteration to get initial matrix
lResSmacof1 <- mds(delta = mDiss, type = "ratio", itmax = 1)

#apply MDS
lResSmacof <- mds(delta = mDiss, type = "ratio", init = lResSmacof1$init)
lResOwn <- fnSMACOF(mDiss, mX = lResSmacof$init)

#plot own results
plot(x = lResOwn$mX[,1], y = lResOwn$mX[,2], xlab = "Axis 1", ylab = "Axis 2")
text(x =  lResOwn$mX[,1], y = lResOwn$mX[,2], colnames(basket), col = "red")

# plot smacof results
plot(x = lResSmacof$conf[,1], y = lResSmacof$conf[,2], xlab = "Axis1", ylab = "Axis 2")
text(x = lResSmacof$conf[,1], y = lResSmacof$conf[,2], colnames(basket), col = "red")

#ordinal results
lResSmacofOrd <- mds(delta = mDiss, type = "ordinal")

plot(x = lResSmacofOrd$conf[,1], y = lResSmacofOrd$conf[,2], xlab = "Axis1", ylab = "Axis 2")
text(x = lResSmacofOrd$conf[,1], y = lResSmacofOrd$conf[,2], colnames(basket), col = "red")
