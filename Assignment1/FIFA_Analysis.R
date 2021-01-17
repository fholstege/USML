#################
# Authors: Floris Holstege, Yuchou Peng, Chao Liang, Eva Mynott, Gabrielle Mignoli
# Purpose: Code for Unsupervised Machine Learning Assignment 1 - applying PCA and sparse PCA to a FIFA dataset
# Sections: 
#           A) Load packages and pre-process the data
#           B) Define our functions for PCA and sparse PCA
#           C) Quickly compare our functions to implementations by standard R packages
#           D) Define functions for visualization ((scree plot, explained variance, loadings, biplot))
#           E) Analyse FIFA data with PCA and sparse PCA, and create relevant plots (scree, explained variance, loadings, biplot)
#################


#######
# A) Load packages and pre-process the data
#######

# Packes required for subsequent analysis. P_load ensures these will be installed and loaded. 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(PMA, 
               ggplot2,
               factoextra,
               xtable,
               tidyverse)

# load the data
load("../Data/FIFA2017_NL.RData")

# load helper functions
source("../Helpers/utilities.R")

# set seed to ensure reproducability
set.seed(123)

# create dataframe with the player attributes
df_fifa_attributes <- fifa[,1:32] %>% filter(Position != "Gk")

# put attributes (without player name, club etc. ) in matrix, and scale these
mX_attributes <- scale(as.matrix(df_fifa_attributes[,4:ncol(df_fifa_attributes)]))


#######
# B) Define our functions for PCA and sparse PCA
#######



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
  
  # define initial matrix R - to be changed as we add v_k, u_k, sigma_k
  mR <- mX
  
  # get initial v
  v.init <- svd_components(mX)$components[,1:K]
  
  # matrices to be filled with values we find per k
  v <- matrix(0, nrow = ncol(mX), ncol = K)
  u <- matrix(0, nrow = nrow(mX), ncol= K)
  sigma <- numeric(0)
  
  # go over k = 1, 2,..K
  for(k in 1:K){
    
    # apply rank-1 algorithm to current version of R
    result_k <- sparse_svd(mR, c1,c2,v=matrix(v.init[,k],ncol = 1))
    
    # get results for this k for v, u, sigma
    v_k <- result_k$v
    u_k <- result_k$u
    sigma_k <- result_k$sigma[1,1]
    
    # add results to matrices to be returned when done
    v[, k] <- v_k
    u[, k] <- u_k
    sigma[k] <- sigma_k
    
    # update matrix R for next k
    mR <- mR - sigma_k*u_k %*% t(v_k)
    
  }
  
  # calculate the explained variance
  VExplained_var <- 1/(nrow(mX)-1)*(sigma^2)[1:k]
  
  # define result and return
  result = list(components = v, u = u,explained_variance=VExplained_var)
  
  return(result)
  
}



#######
# C) Quickly compare our functions to implementations by standard R packages
#######


### Compare our PCA function to the prcomp() function in R
# We get the exact same results as prcomp 
result_prcomp <- prcomp(mX_attributes)
result_PCA <- svd_components(mX_attributes)
result_prcomp$rotation
result_PCA$components

### Compare our sparse PCA function to the PMD() function in R, for rank-1
## the results appear very, very similar - when sumabsu, sumabsv = 1, then difference of 1-e16. No difference when other values are picked.
# our sparse PCA
result_sparse_svd <- sparse_svd(mX_attributes, c1=1.5, c2=1.5)
result_sparse_svd$v
result_sparse_svd$u

# sparce PCA from PMA package
result_pmd <- PMD(mX_attributes, sumabsu = 1.5, sumabsv = 1.5, center = FALSE)
result_pmd$v
result_pmd$u

### Compare our sparse PCA function to the PMD() function with a higher rank 
## not exactly the same when ranks are higher - but this is likely because neither methods find a global maximum
# our sparse PCA for rank > 1
result_sparse_PCA <- sparse_PCA(mX_attributes, c1= sqrt(nrow(mX_attributes)),c2=2,K=5)
result_sparse_PCA$components
result_sparse_PCA$u

# PMD implementation for rank > 1
result_pmd_rank <- PMD(mX_attributes, sumabsu = 1.5, sumabsv = 1.5, K=5)
result_pmd_rank$v
result_pmd_rank$u


#######
# D) Define functions for visualization ((scree plot, explained variance, loadings, biplot))
#######


####
# get_df_var: quick helper function for create_screePlot & create_varPlot, gets variance from PCA result and makes it fit for ggplot
# arguments: 
#     type: string, either PMD or SPC
# 
# output: 
#     df_var: dataframe used in ggplot viz

get_df_var <- function(type, result_PCA){
  if(type == "PMD"){
    
    df_var <- data.frame(explained_variance = result_PCA$explained_variance)
    
  }else if(type == "SPC"){
    
    n <- nrow(result_PCA$u)
    p <- nrow(result_PCA$v)

    d <- result_PCA$d
    explained_variance <- 1/(n-1)*(d^2)[1:p]
    df_var <- data.frame(explained_variance = explained_variance)
    
  }
  
  return(df_var)
}

####
# create_screePlot: creates screeplot for results of PCA
# arguments: 
#     result_PCA: list, created by svd_components or sparse_PCA
# 
# output: 
#     screePlot: a ggplot object with the screeplot

create_screePlot <- function(result_PCA, type="PMD"){
  
  df_var <- get_df_var(type, result_PCA)
  
  # create screeplot
  screePlot <- ggplot(data=df_var,aes(x=seq(1,length(explained_variance),1), y=explained_variance))+
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

create_varPlot <- function(result_PCA, type = "PMD"){
  
  df_var <- get_df_var(type, result_PCA)
  
  # create screeplot
  varPlot <-  ggplot(data=df_var, aes(x=seq(1,length(explained_variance),1), y=cumsum(explained_variance/sum(explained_variance))))+
    geom_line(linetype="dashed",color="red", size=2)+
    geom_point(color="blue", size=4)+
    labs(x = "Component", y="% Variance Explained")+
    theme_minimal() +
    theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 16),
          axis.text.y = element_text(size = 10), axis.title.y = element_text(size = 16))+
    scale_y_continuous(labels = function(x) paste0(x*100, "%"))
  
  return(varPlot)
  
}

####
# create_loadingTable: creates latex code for loadings of a PCA result
# arguments: 
#     result_PCA: list, created by svd_components or sparse_PCA
# 
# output: 
#     varPLot: a ggplot object with the plot showing the cumulative variance per added component

create_LoadingTable <- function(mX,result_PCA, K, type="PMD", latex=TRUE){
  
  # get variable names
  lVariables <-colnames(mX)
  
  # get components based on type
  if(type=="PMD"){
    components <- result_PCA$components[,1:K]
  }else if (type=="SPC"){
    components <- result_PCA$v[,1:K]
  }
  
  # create dataframe with variable names and loadings
  dfLoadings <- data.frame(lVariables,components)
  colnames(dfLoadings) <- c("Variable", paste0("PC ", seq(1,K,1)))
  
  # create latex table
  if(latex){
    xtable(dfLoadings, type="latex")
  }

  return(dfLoadings)
  
}

PCbiplot <- function(PC, x="PC1", y="PC2",type="PMD",mX, rownames = 1:nrow(mX)) {
  
  varnames <- colnames(mX)
  
  if(type=="PMD"){
    data <- data.frame(obsnames=rownames, PC$x)
    components <- PC$rotation

  }else if(type=="SPC"){
    mPC <- data.frame(mX %*% PC$v)
    colnames(mPC) <- c(paste0("PC", seq(1,ncol(mPC),1)))
    
    data <- data.frame(obsnames=rownames, mPC)
    components <- data.frame(PC$v)
    colnames(components) <- c(paste0("PC", seq(1,ncol(mPC),1)))
  }
  


  plot <- ggplot(data, aes_string(x=x, y=y)) + geom_text(alpha=.4, size=3, aes(label=obsnames))
  plot <- plot + geom_hline(yintercept=0, size=.2) + geom_vline(xintercept=0, size=.2)
  datapc <- data.frame(varnames=varnames, components)
  mult <- min(
    (max(data[,y]) - min(data[,y])/(max(datapc[,y])-min(datapc[,y]))),
    (max(data[,x]) - min(data[,x])/(max(datapc[,x])-min(datapc[,x])))
  )
  datapc <- transform(datapc,
                      v1 = .7 * mult * (get(x)),
                      v2 = .7 * mult * (get(y)) 
  )%>% filter(!(v2 == 0 & v1 ==0))
  plot <- plot + coord_equal() + geom_text(data=datapc, aes(x=v1, y=v2, label=varnames), size = 5, vjust=1, color="red")
  plot <- plot + geom_segment(data=datapc, aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="red")
  plot <- plot + theme_minimal()
  plot
}




#######
# E) Analyse FIFA data with PCA and sparse PCA, and create relevant plots (scree, explained variance, loadings, biplot)
#######

### First, for standard PCA - we use our own implementation, since identical

### create an scree plot
screePlot_PCA <- create_screePlot(result_PCA, type="PMD")
screePlot_PCA

## show how variance improves with each component added
varPlot_PCA <- create_varPlot(result_PCA, type="PMD")
varPlot_PCA

## create table of the loadings 
create_LoadingTable(mX_attributes, result_PCA, K=5, type="PMD")

## create biplot
PCbiplot(result_prcomp, mX=mX_attributes, rownames=df_fifa_attributes$name)

## results for sparse PCA - proceed with result from package, using SPC() which is in line with Hui Zou, Trevor Hastie, and Robert Tibshirani 

# First, use cross validation to determine optimal penalty for sparse PCA
result_cv <- SPC.cv(mX_attributes, sumabsvs = seq(1,sqrt(ncol(mX_attributes)), 0.5))
result_cv

result_sparse_PCA_SPC <- SPC(mX_attributes, sumabsv = 1.5, K = 29)

### create an scree plot
screePlot <- create_screePlot(result_sparse_PCA_SPC,type="SPC")
screePlot

## show how variance improves with each component added
varPlot <- create_varPlot(result_sparse_PCA_SPC, type="SPC")
varPlot

## create table of the loadings 
create_LoadingTable(mX_attributes,result_sparse_PCA_SPC, K=6, type="SPC")

## create biplot
PCbiplot(result_sparse_PCA_SPC, mX = mX_attributes,type="SPC",rownames=df_fifa_attributes$name)

