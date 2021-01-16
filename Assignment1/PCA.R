

####
# svd_components: extracts the principal components and their explained variaance from a matrix X
# arguments: 
#     mX: a matrix object
# output: 
#     result: a list object with the components and their explained variance
svd_components <- function(mX){
  
  # parameters; k columns, n observations, vector i of ones
  k = ncol(mX)
  n <- dim(mX)[1]
  i <- rep(1,n)
  
  # use this matrix H to center matrix X
  mH <- (diag(n) -  (1/n) * i %*% t(i))
  mX_center <-mH %*% mX
  
  # get use single value decomposition (svd)
  lSvd_mX <- svd(mX_center)
  
  # extract Principal components and explained variance from single value decomposition
  mPComponents <- lSvd_mX$v[,1:k]
  vExplained_var <- 1/(n-1)*(lSvd_mX$d^2)[1:k]
  
  # create list of results, with components and explained var for each
  result <- list(components=mPComponents,     
                 explained_variance=vExplained_var,
                 d = lSvd_mX$d)
  
  return(result)
}


####
# sparce_svd: applies sparce svd to a matrix X
# arguments: 
#     mX: a matrix object
#     c1: a constant, influences penalty term by constraining L1 norm of vector v 
#     c2: a constant, influences penalty term by constraining L1 norm of vector u
#     iMax: default = 100, maximum of iterations
#     v: optional, pass on which to take
# output: 
#     result: a list object with the vector v, u, and sigma 

sparse_svd<- function(mX, c1, c2,iMax=100, v=NULL){
  
  if(is.null(v)){
    v <- as.matrix(svd_components(mX)$components[,1])
  }
  
  # start while loop
  i = 0
  while(i <= iMax){
    
    # update each iteration
    i <- i + 1
    
    # define matrix Xv
    mXv <- mX %*% v
    
    lambda1 <- binary_search(mXv, c1)
    
    # get u, given V
    u  <- soft_l2norm(mXv, lambda1)
    
    # define matrix Xu
    mXu <- t(mX) %*% u
    
    # then, determine lambda 2
    lambda2 <- binary_search(mXu, c2)
    
    # redefine v
    v <- soft_l2norm(mXu, lambda2)
    
  }
  
  # calculate sigma
  sigma = t(u) %*% mX %*% v
  
  # return the results
  results <- list(v = v, u = u, sigma = sigma)
  return(results)
}



####
# sparce_PCA: applies sparce PCA to a matrix X, returns rank K
# arguments: 
#     mX: a matrix object
#     K : Number of principal components to return
#     c1: a constant, influences penalty term by constraining L1 norm of vector v 
#     c2: a constant, influences penalty term by constraining L1 norm of vector u
# output: 
#     result: a list object with the vector v, u, and sigma 

sparse_PCA <- function(mX, c1, c2, K){
  
  mR <- mX
  
  v.init <- svd_components(mX)$components[,1:K]
  v <- matrix(0, nrow = ncol(mX), ncol = K)
  u <- matrix(0, nrow = nrow(mX), ncol= K)
 

  for(k in 1:K){
    
    result_k <- sparse_svd(mR, c1,c2,v=matrix(v.init[,k],ncol = 1))
    
    v_k <- result_k$v
    u_k <- result_k$u
    sigma_k <- result_k$sigma
    
    v[, k] <- v_k
    u[, k] <- u_k
    
    
    mR <- mR - sigma_k[1,1]*u_k %*% t(v_k)

  }
  
  result = list(v = v, u = u, v.init=v.init)
  
  return(result)
  
}





####
# create_screePlot: creates screeplot for results of PCA
# arguments: 
#     result_PCA: list, created by svd_components or sparse_PCA
# 
# output: 
#     screePlot: a ggplot object with the screeplot

create_screePlot <- function(result_PCA){
  
  # create screeplot
  screePlot <- ggplot(data=data.frame(result_PCA), aes(x=seq(1,length(explained_variance),1), y=explained_variance))+
                      geom_line(linetype="dashed",color="red", size=2)+
                      geom_point(color="blue", size=4)+
                      labs(x = "Component", y="Eigenvalue")+
                      theme_minimal() +
                      theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 16),
                            axis.text.y = element_text(size = 10), axis.title.y = element_text(size = 16))  +
                      geom_hline(yintercept = 1, linetype="dotted", 
                                  color = "black", size=1.5)
  return(screePlot)
  
  
}

####
# create_screePlot: creates plot that shows the % of variance explained for a PCA
# arguments: 
#     result_PCA: list, created by svd_components or sparse_PCA
# 
# output: 
#     varPLot: a ggplot object with the plot showing the cumulative variance per added component


create_varPlot <- function(result_PCA){
  
  # create screeplot
  varPlot <-  ggplot(data=data.frame(result_PCA), aes(x=seq(1,length(explained_variance),1), y=cumsum(explained_variance/sum(explained_variance))))+
    geom_line(linetype="dashed",color="red", size=2)+
    geom_point(color="blue", size=4)+
    labs(x = "Component", y="% Variance Explained")+
    theme_minimal() +
    theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 16),
          axis.text.y = element_text(size = 10), axis.title.y = element_text(size = 16))+
    scale_y_continuous(labels = function(x) paste0(x*100, "%"))
  
  return(varPlot)
  
  
}

create_LoadingTable <- function(result_PCA, K){
  
  lVariables <-colnames(mX_attributes)
  
  dfLoadings <- data.frame(lVariables,result_PCA$components[,1:K])
  
  colnames(dfLoadings) <- c("Variable", paste0("PC ", seq(1,K,1)))
  
  xtable(dfLoadings, type="latex")
  
  
}





