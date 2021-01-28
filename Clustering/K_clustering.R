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

###############################################################################
# C) K-means
##############################################################################

K_means <- function(mX, K, n_iter = 10, n_random_centroids=10, eps = 1e-6){
  
  # counter for number of random points to pick
  j = 1
  
  # save the results here
  results_centroid <- list()
  
  # total sum of square
  curr_tot_within_SS = Inf
  
  while(j <= n_random_centroids){
    
    # update
    j = j + 1
    
    # pick the initial centroids 
    centroids <- mX[sample(nrow(mX), K), ]
    
    # counter for iterations once random points picked
    i = 1

    while(i <= n_iter){
      
      # update
      i = i + 1
      
      # create distance matrix between centroids and all observations
      mDist <- calc_sq_euc_dist_m(mX, centroids)
      
      # get the clusters for each centroid
      clusters <- apply(mDist, 1, FUN=which.min)
      
      # check cluster per observation
      mX_cluster <- cbind(mX, clusters)
      
      # get new centroids 
      centroids <- aggregate(mX_cluster,by= list(mX_cluster[,'clusters']), mean)[,-1]
    }
    
    # get the final distance matrix
    mDist_final <- calc_sq_euc_dist_m(mX, centroids)
    
    # get the minimum distance per observation, combined with cluster
    mX_dist_c <- cbind(apply(t, 1, FUN=min), clusters)
    
    # current within sums and total sums
    within_SS <- aggregate(mX_dist_c,by= list(mX_dist_c[,'clusters']), sum)
    
    tot_within_ss <- sum(within_SS)

    if (tot_within_ss < curr_tot_within_SS){
      curr_tot_within_SS <- tot_within_ss
    }
  }
  
  result = list(centroids = centroids, tot_within_SS = curr_tot_within_SS, within_SS = within_SS)
  
  return(result)
}

###############################################################################
# D) Compare to package
##############################################################################

K = 3
mX <- as.matrix(cityweather)

result_ours <- K_means(as.matrix(cityweather), K, n_random_centroids = 10)
result_ours$centroids
result_ours$within_SS


result_pack <- kmeans(mX, K, iter.max=10, nstart=100)
result_pack$centers
result_pack$withinss


