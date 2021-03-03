######################
# Author: Floris Holstege
# Purpose: Set of functions to run a hierarchical clustering algorithm, and compare it to hclust()
# 
# Overview:
#           A) Loading the packages and cleaning the data, defining the distance matrix
#           B) Set of functions used for hierarchical clustering algorithm
#           C) Compare our hierarchical clustering algorithm to hclust()
#
#####################

#####################
# A) Loading the packages and cleaning the data, defining the distance matrix
#####################

# empty global environment
rm(list=ls())

# Packes required for subsequent analysis. P_load ensures these will be installed and loaded. 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(pdist,
               tidyverse) 

# load data to compare
load("Data/cityweather.RData")

cityweather <- head(cityweather,20)

# scale the data to be between 0, 1
distObj <- dist(scale(cityweather))
mDist <- as.matrix(distObj)


######################
# B) Set of functions used for hierarchical clustering algorithm
######################


#####################
# calc_method_cluster: calculates a distance between an nxk matrix of n observations and k clusters, based on the method
#
# Arguments: 
#     mDist_labels: nxk matrix of n observations and k clusters, 
#     method: string, one of "single" (minimum), "complete" (maximum) or "average" (avg. distance)
#
# Output:
#     mDist_labels_calced: kxk matrix of the distance between the labels
#     
#
calc_method_cluster <- function(mDist_labels,method){
  
  # if single, get  observation with minimum distance
  if(method == "single"){
    mDist_labels_calced <- mDist_labels %>% summarise_all(min)
  # if complete, get observation with maximum distance  
  }else if (method == "complete"){
    mDist_labels_calced <- mDist_labels %>% summarise_all(max)
  } else if (method == "average"){
  # if average, get average distance of the observations in a cluster
    mDist_labels_calced <- mDist_labels %>% summarise_all(mean)
  }
  # remove labels for cluster and return
  mDist_labels_calced$labels <- NULL
  
  return(mDist_labels_calced)
}

#####################
# update_mDist_clusters: turns an nxn distance matrix to an kxk distance matrix of the distances between clusters
#
# Arguments: 
#     mDist: nxn matrix with distance between points
#     vLabels: vector with labels per observation
#     method: string, one of "single" (minimum), "complete" (maximum) or "average" (avg. distance)
#
# Output:
#     dfDist_labels_clusters: kxk matrix of the distance between the labels
#     
#
update_mDist_clusters <- function(mDist, vLabels, method = "single"){
  
  print("AT START")
  print(dim(mDist))

  # group by clusters
  dfDist_labels <- data.frame(cbind(labels = vLabels,mDist)) %>% group_by(labels)
  
  # calculate the distance between each point and the clusters
  dfDist_labels_Obs_cluster <- calc_method_cluster(dfDist_labels, method)
  
  # transpose this matrix and determine which clusters belong to each point
  dfDist_labels_Obs_cluster_t <-  data.frame(cbind(labels = vLabels,t(dfDist_labels_Obs_cluster))) %>% group_by(labels)
  
  
  dfDist_labels_clusters <- calc_method_cluster(dfDist_labels_Obs_cluster_t, method)
  
  print("AT END")
  print(dim(dfDist_labels_clusters))
  
  colnames(dfDist_labels_clusters) <- unique(vLabels)
  rownames(dfDist_labels_clusters) <- unique(vLabels)
  
  
  return(dfDist_labels_clusters)
  
}
#####################
# Hierarchical_clustering: using a distance matrix, applies hierarchical clustering method
#
# Arguments: 
#     mDist: nxn matrix with distance between points
#     method: string, one of "single" (minimum), "complete" (maximum) or "average" (avg. distance)
#
# Output:
#     dfDist_labels_clusters: kxk matrix of the distance between the labels
#     
#


hierarchical_clustering <- function(mDist, method = "single"){
  
  # turn distance matrix to lower triangular matrix
  mDist[upper.tri(mDist, diag = TRUE)]<- NA
  
  # save parameters - number of observations
  n <- nrow(mDist)
  
  # save here the updated distance matrix - distance between clusters
  mDist_current <- mDist
  
  # save here the allocated clusters per iteration - starting at n, then n_1, n_2, ..
  dfLabels <- data.frame(obs = rownames(mDist), n = 1:nrow(mDist))
  
  
  # iteratively add clusters
  for(i in n:1){
    
    # get parameters - n of clusters
    n_cluster <- n - i + 1
    col_name_labels <- paste0("n_", n_cluster)

    # update the distance matrix
    mDist_current <- update_mDist_clusters(mDist, vLabels = dfLabels[,n_cluster+1], method = method)
    
    # get obs where the distance is the lowest 
    dfWhich_min = which(mDist_current == min(mDist_current, na.rm = TRUE), arr = TRUE)
    dfWhich_min <- matrix(dfWhich_min[!duplicated(dfWhich_min[,1]),], ncol = 2)
    
    print("df to check")
    print(dfWhich_min)
        
    
    # get index of the observations, and the cluster they are assigned to
    index_cluster_obs1 = c(dfWhich_min[,1]) 
    index_cluster_obs2 = c(dfWhich_min[,2])
    
    print("index label input")
    print(index_cluster_obs1)
    
    print("index label replace")
    print(index_cluster_obs2)
    
    cluster_obs1 <- dfLabels[,n_cluster+1][index_cluster_obs1]
    cluster_obs2 <- dfLabels[,n_cluster+1][index_cluster_obs2]
    
    
    print("Label input")
    print(cluster_obs1)
    
    print("Label replaced")
    print(cluster_obs2)
    
    print("DISTANCE")
    print(mDist_current)
    
    # Update labels for this iteration
    dfLabels[,n_cluster+2] <- dfLabels[,n_cluster + 1]
    index_toChange <- which(cluster_obs2 == dfLabels[,n_cluster + 1])
    dfLabels[index_toChange,n_cluster+2] <- cluster_obs1

    
    # set col name
    colnames(dfLabels)[2 + n_cluster] <- col_name_labels
    
  }
  
  return(dfLabels)
  
}

our_result<- hierarchical_clustering(mDist, method = "single")
hclust_result  <- hclust(distObj, method = 'single')


theirs = cutree(hclust_result, k = 18)
dfComparison <- data.frame(ours = our_result$n_2, theirs = theirs)


