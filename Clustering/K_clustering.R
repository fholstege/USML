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
               tidyverse, 
               dplyr, 
               tidyr) 

# load the data
load("Data/cityweather.Rdata")

# set seed
set.seed(123)

################################################################################
# Data prep functions
################################################################################

# Clean the data
data.cleaning <- function(data, dup.cols) {
  
  # Enable variable calling
  attach(data)
  
  # check if there are missing values
  if (any(is.na(data)) == T) {
    na.rows <- length(which(apply(data, 1, function(X) any(is.na(X)))))
    cat('Data set contains', sum(is.na(data)),'missing values, in', 
        na.rows, 'rows.')
    
    # Remove rows with missing values (if less than 1% of data)
    if (na.rows/nrow(data) < 0.01) {
      data <- na.omit(data)
    } else {
      print('\n\nWarning: too many missing values removed.')
    }}
  
  # check for duplicates and remove them
  duplicate_indexes <- which(duplicated(data[dup.cols]),) 
  data <- data[-duplicate_indexes,]
  
  return(data)
}

# Select variables
data.selection <- function(data, select.cols) {
  
  # store selected variables in dataframe
  data <- as.data.frame(data)
  data <- subset(data, select = select.cols)
  
  return(data)
}

# Standardize the data (if not same unit of measurement)
data.scaling <- function(data, method = "stand", select.vars) {
  
  # select subset to be scaled
  vars <- subset(data, select = select.vars)
  
  # scale using chosen method
  if (method == "minmax"){
    vars <- apply(vars, 2, function (x) (x-min(x))/(max(x)-min(x)))
  } else if (method == "stand") {
    vars < scale(vars, center = T, scale = T)
  } else {
    print("Specify: method = c(stand, minmax). Stand used by default.")
  }
  
  # obtain all data incl. scaled vars
  data <- data[,!(names(data) %in% select.vars)]
  data <- cbind(data,vars)
  
  return(data)
}

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
            
### This is how I called the functions
### Required packages: dplyr

# Obtain cleaned data
dup.cols <- c('track_name', 'track_artist')
data <- data.cleaning(data, dup.cols)

# Obtain filtered data
data <- filter(data, track_popularity > 0)
select.cols <- c('danceability', 'energy', 'loudness', 'speechiness',
                 'acousticness', 'instrumentalness', 'liveness', 'valence',
                 'tempo')
data <- data.selection(data, select.cols)

# Adapt variables that need scaling
select.vars <- c('loudness', 'tempo')
data <- data.scaling(data, method = "minmax", select.vars)
