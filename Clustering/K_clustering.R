#################
# Authors: Floris Holstege, Yuchou Peng, Chao Liang, Eva Mynott, Gabriele Mingoli
# Purpose: Code for Unsupervised Machine Learning Assignment 2 - applying K-means and K-medioids
# Sections: 
#           A) Load packages and data
#           B) Helper functions
#           C) K-means implementation
#           D) Compare to package  
#################


################################################################################
# A) Load packages & Data
################################################################################

# empty global environment
rm(list=ls())

# Packes required for subsequent analysis. P_load ensures these will be installed and loaded. 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(pdist,
               tidyverse) 

# load the data
load("Data/cityweather.Rdata")

# set seed
set.seed(123)

###############################################################################
# B) Helper Functions
###############################################################################

sq.euc.dist <- function(x1, x2){sum((x1 - x2)^2)}

calc_sq_euc_dist_m <- function(mX, centroids){ as.matrix(pdist(mX,centroids))^2}

calc_within_ss <- function(mDist, clusters){ 
  
  # get the minimum distance per observation, combined with cluster
  mX_dist_c <- cbind(apply(mDist, 1, FUN=min), clusters)
  
  # calculate within sum of squares
  within_SS <- aggregate(mX_dist_c,by= list(mX_dist_c[,'clusters']), sum)
  
  return(within_SS)
}

###############################################################################
# C) K-means
##############################################################################

K_means <- function(mX, K, n_iter = 10, n_random_centroids=10, eps = 1e-6){
  
  # parameters
  n <- nrow(mX)
  
  # counter for number of random points to pick
  j = 1
  
  # total sum of square
  curr_tot_within_SS = Inf
  
  # placeholder for saving centroids
  fin_centroids <- matrix(0, K, ncol(mX))
  
  fin_within_SS <- matrix(0,1,K)
  

  while(j <= n_random_centroids){
    
    # update
    j = j + 1
    
    # pick the initial centroids
    centroids <- mX[sample(n, K), ]

    # create list of 0's for clusters to continue checking
    clusters <- rep(0,n)
    
    # counter for iterations once random points picked
    i = 1

    while(i <= n_iter){
      
      # update
      i = i + 1
      
      # create distance matrix between centroids and all observations
      mDist <- calc_sq_euc_dist_m(mX, centroids)
      
      # get the clusters for each centroid
      new_clusters <- apply(mDist, 1, FUN=which.min)

      # check - if none of the clusters changed, stop
      if(any(clusters != new_clusters)){
        clusters <- new_clusters
      }else{break}
      
      # check cluster per observation
      mX_cluster <- cbind(mX, clusters)
      
      # get new centroids 
      centroids <- aggregate(mX_cluster,by= list(mX_cluster[,'clusters']), mean)[,-1]
    }
    
    # get the final distance matrix
    mDist_final <- calc_sq_euc_dist_m(mX, centroids)
    
    # current within sums and total sums
    within_SS <- calc_within_ss(mDist_final, clusters)
    tot_within_ss <- sum(within_SS)

    # if the current within variation is less than previous, update
    if (tot_within_ss < curr_tot_within_SS){
      curr_tot_within_SS <- tot_within_ss
      fin_centroids <- centroids
      fin_within_SS <- within_SS
      
    }
  }
  
  result = list(centroids = fin_centroids, tot_within_SS = curr_tot_within_SS, within_SS = fin_within_SS)
  
  return(result)
}

###############################################################################
# D) Compare to package
##############################################################################

K = 3
mX <- as.matrix(cityweather)


result_ours <- K_means(as.matrix(cityweather), K,n_iter = 10, n_random_centroids = 10)
result_ours$centroids
result_ours$within_SS


result_pack <- kmeans(mX, K, iter.max=10, nstart=10)
result_pack$centers
result_pack$withinss

