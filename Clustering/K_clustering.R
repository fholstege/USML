#################
# Authors: Floris Holstege, Yuchou Peng, Chao Liang, Eva Mynott, Gabriele Mingoli
# Purpose: Code for Unsupervised Machine Learning Assignment 2 - applying K-means and K-medioids
# Sections: 
#           A) Load packages and data
#           B) Data prep functions
#           C) Helper functions for K-means implementation
#           C) K-means implementation
#           D) Compare to our implementation package 'kmeans'
#           E) Clean spotify data for analysis 
#           F) Find optimal K, Lambda, and random starts 
#           G) Analyse the results
#
#
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
               tidyr,
               readr,
               psych,
               tidyverse,
               reshape2,
               GGally,
               clustMixType,
               fastDummies,
               parallel,
               devtools,
               xtable) 

# add repository to make mclapply() run in parallel (only necessary on windows)
install_github('nathanvan/parallelsugar')
library(parallelsugar)


# load the raw data on spotify songs
spotify_songs <- read_csv('https://raw.githubusercontent.com/rfordatascience/tidytuesday/master/data/2020/2020-01-21/spotify_songs.csv')

# set seed
set.seed(123)

################################################################################
# Data prep functions
################################################################################

####
# data.cleaning: Removes NA's and duplicate rows 
# 
#  Arguments:
#     data: dataframe
#     dup.cols: columns from which to check if there is a duplicate
# 
# Output:
#     data: same dataframe, without NA's and duplicates
####
data.cleaning <- function(data, dup.cols) {
  
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


####
# data.scaling: scales chosen columns
# 
#  Arguments:
#     data: dataframe
#     method: string, one of "minmax" or "stand"
#     select.vars: vars to scale
# 
# Output:
#     data: same dataframe, with scaled columns
####

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
# B) Helper Functions for K_means
###############################################################################

####
# calc_sq_euc_dist_m: calculates squared euclidean distance matrix between X and centroids
# 
#  Arguments:
#     mX: matrix of characteristics
#     centroid: matrix of centroids
# 
# Output:
#     n x K matrix of squared euclidean distances
####

calc_sq_euc_dist_m <- function(mX, centroids){ as.matrix(pdist(mX,centroids))^2}

####
# calc_within_ss: calculates within sum of squares for a given distance matrix
# 
#  Arguments:
#     mDist: n x K matrix of distances
#     clusters: nx 1 vector of assigned clusters 
# 
# Output:
#    k floats, within ss
####

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

####
# K_means: implements K_means algorithm
# 
#  Arguments:
#     mX: matrix of characteristics
#     K: number of clusters
#     n_iter: number of iterations before stopping
#     n_random_centroids: the number of random starts to try
# 
# 
# Important: 
#     Current stopping criterion, apart from max iterations, is if there are no changes in assigned clusters     
# 
# Output:
#     n x K matrix of squared euclidean distances
####

K_means <- function(mX, K, n_iter = 10, n_random_centroids=10){
  
  # parameters
  n <- nrow(mX)
  
  # counter for number of random points to pick
  j = 1
  
  # total sum of square
  curr_tot_within_SS = Inf
  
  # placeholder for saving centroids
  fin_centroids <- matrix(0, K, ncol(mX))
  fin_within_SS <- matrix(0,1,K)
  fin_clusters <- rep(0, n)
  
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
      new_clusters <- apply(mDist, MARGIN=1, FUN=which.min)

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
      fin_clusters <- clusters
      
    }
  }
  
  # save results in a list
  result = list(centroids = fin_centroids, tot_within_SS = curr_tot_within_SS, within_SS = fin_within_SS, clusters = fin_clusters)
  
  return(result)
}

###############################################################################
# D) Compare our implementation to package
##############################################################################

# show for K = 3
K = 3

load("Data/cityweather.Rdata")
mX_city <- as.matrix(cityweather)

# check our results 
result_ours <- K_means(mX_city, K,n_iter = 10, n_random_centroids = 100)
result_ours$centroids
result_ours$within_SS

dfResult <- cbind(result_ours$centroids, result_ours$within_SS)[-c(6,7)]
xtable(dfResult)

# check results package
result_pack <- kmeans(mX_city, K, iter.max=10, nstart=100)

xtable(cbind(result_pack$centers, result_pack$withinss))

help(kmeans)

###############################################################################
# E) Clean Spotify Data
##############################################################################

# Obtain cleaned data, with duplicates removed
dup.cols <- c('track_name', 'track_artist')
spotify_clean <- data.cleaning(spotify_songs, dup.cols)

# select columns
select.cols.con <- c('danceability', 'energy', 'loudness', 'speechiness',
                        'acousticness', 'instrumentalness', 'liveness', 'valence',
                        'tempo')

# create dummies for key variable
vKey <- spotify_clean %>% select('key')
dummy_key <- dummy_cols(vKey)[,-1]
dummy_key_factor <- data.frame(lapply(dummy_key, factor))

# select data for kmeans, and kproto
spotify_selected <- spotify_clean %>% select(all_of( select.cols.con))
spotify_selected_cat <- cbind(spotify_selected, mode = spotify_clean$mode,dummy_key_factor) 

# Adapt variables that need scaling
select.scale.vars <- c('loudness', 'tempo')
spotify_final <- data.scaling(spotify_selected, method = "minmax", select.scale.vars)
spotify_final_cat <- data.scaling(spotify_selected_cat, method = "minmax", select.scale.vars)

# turn to matrix for kmeans
mX_spotify <- as.matrix(spotify_final)

# set mode and key as factors for the kproto
spotify_final_cat$mode <- as.factor(spotify_final_cat$mode)

###############################################################################
# F) Determine optimal K and number of random starts
##############################################################################


######
# find_K_elbow: generate plot with within ss per K 
# 
# Arguments; 
#     mX: matrix X of attributes, numeric continuous
#     K_range; 1:max_k, try these values of K
#     n_random_centroids: number of random starts
#     type: string, one of 'kmeans' or 'kproto'
#     lambda: if type='kproto', then set lambda
#   
#  output: ggplot object with elbowplot
######
find_K_elbow <- function(mX, K_range, n_iter, n_random_centroids, type = "kmeans", lambda = NULL){
  
  # list of predefined length with the within ss
  lResults <-  vector(mode = "list", length = length(K_range))
  result_index <- 1

  # go over all K's and get the results
  for(K in K_range){
    
    # pick kmeans or kproto based on the user specification
    if(type == "kmeans"){
      result_k <- kmeans(mX, K, n_iter, n_random_centroids)$tot.withinss
    }else if (type == "kproto"){
      result_k <- kproto(mX, K, n_iter, n_random_centroids,lambda=lambda)$tot.withinss
    }else{
      stop("Specify one of 'kmeans' or 'kproto' for 'type' variable")
    }
    
    # add to predefined list of results
    lResults[[result_index]] <- result_k
    result_index <- result_index + 1
    
  }
  
  # get list of results, create dataframe with K and total within ss
  ltot.withinss <-  unlist(lResults)
  dfElbowPlot <- data.frame(K = K_range_spotify, tot.withinss = ltot.withinss)
  
  # define the elbow plot
  elbowPlot <- ggplot(data = dfElbowPlot, aes(x = K, y =tot.withinss))+
    geom_line(color = "Red", size = 2) + 
    labs(x = "K", y = "Total Within Cluster Distance")+
    theme_bw() + 
    theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 16),
          axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16)) 
  
  return(elbowPlot)
  
  
} 


######
# check_initialPoint_range: checks distribution of the scores one gets for n starts with kproto 
# 
# Arguments:
#     dfX: dataframe X
#     K: number of clusters
#     iter.max: number of max iterations
#     nstart: number of simulations, standard = 1000
#
# Output:`
#     ltot.withinss: list of length nstart of total within sum of squares for each simulation
######

check_initialPoint_range <- function(dfX, K, lambda, iter.max, nstart=1000){
  
  range <- 1:nstart
  lResults = vector(mode = "list", length = length(nstart))
  
  get_initial_result <- function(i,dfX, K, iter.max, lambda, l){l[[i]] <- kproto(dfX, K, iter.max = iter.max, nstart = 1,lambda = lambda)$tot.withinss}
  
  par_result <- mclapply(range, get_initial_result, dfX=dfX, K=K, iter.max=iter.max, lambda=lambda, l=lResults, mc.cores = detectCores()/2)
  
  ltot.withinss <-  unlist(par_result)
  
  return(ltot.withinss)
  
}


########
# Check how many random starts we should have
########

# get random sample of 5000 observations from datsaet
spotify_final_limitedSample <- spotify_final_cat[sample(nrow(spotify_final_cat), 5000), ]

# get the total within ss for 1000 random starts
sim_tot.withinss <- check_initialPoint_range(spotify_final_limitedSample, K_ideal,lambda, iter.max = 100, nstart = 1000)

# create plot that shows distribution of results for these 1000 random starts
Convergence_Plot <- ggplot() + 
  geom_histogram(aes(x=sim_tot.withinss), stat='bin', colour="black", fill="#3090C7", bins = 100)+
  labs(x = "Total Within Cluster Distance", y = "Count")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16))+
  lims(x=c(1300,1500))

########
# Find Lambda
########


# check the average standard deviation
sigma_col <- apply(mX_spotify, 2, sd)
avg_sigma <- mean(sigma_col)
avg_sigma

# paper says - suitable gamma is between 1/3*sigma and 2/3* sigma
lambda = 0.5 * avg_sigma

########
# Check optimal K
########


# range of values for K - check between 2 and 10
K_range_spotify <- 2:10
elbowPlot_kmeans <- find_K_elbow(mX_spotify, K_range = K_range_spotify, n_iter = 100, n_random_centroids = 10, type='kmeans')
elbowPlot_kproto <- find_K_elbow(spotify_final_cat, K_range = K_range_spotify, n_iter = 100, n_random_centroids = 10, type='kproto', lambda = lambda)
elbowPlot_kmeans + geom_vline(xintercept = 4, linetype="dotted", color = "black", size=1.5)
elbowPlot_kproto + geom_vline(xintercept = 4, linetype="dotted", color = "black", size=1.5)  

###############################################################################
# G) Analysis of results
##############################################################################

# define ideal K based on elbow plot, run the kmeans and kproto results
K_ideal <- 4
result_kmeans <- kmeans(mX_spotify, K_ideal, iter.max = 100, nstart = 50)
result_kproto <- kproto(spotify_final_cat, K_ideal, iter.max = 100, nstart = 50,lambda = lambda, keep.data=TRUE)


# add the dataframes of the results of the cluster analysis and the variables
dfSpotify_clustered <- data.frame(cbind(spotify_final_cat, 
                                  cluster_kmeans = result_kmeans$cluster), 
                                  cluster_kproto = result_kproto$cluster )



##############
# First, lets grab several summary statistics for both kmeans and kproto
##############


# summarise the average values for k-means
dfSummary_kmeans <- dfSpotify_clustered %>% 
  group_by(cluster_kmeans) %>%
  summarise_all(mean) %>%
select(-starts_with("key_"), - mode, -cluster_kproto)
dfSummary_kmeans - dfSummary_kproto[c(1,3,4,2),]

# summarise the average values for the k-prototype
dfSummary_kproto <- dfSpotify_clustered %>% 
  group_by(cluster_kproto) %>%
  summarise_all(mean) %>%
  select(-starts_with("key_"), - mode, -cluster_kmeans)
dfSummary_kproto

# summarise the percentage for the mode variable per cluster
dfSummary_kproto_mode <- dfSpotify_clustered %>% 
  group_by(cluster_kproto, mode) %>%
  summarise(n_mode = n()) %>%
  mutate(perc_mode = (n_mode/sum(n_mode)))
dfSummary_kproto_mode


# summarise the percentage for the key variable(s) per cluster
dfSummary_kproto_key <- dfSpotify_clustered %>% 
  select(starts_with("key_"), cluster_kproto) %>%
  melt(id = 'cluster_kproto') %>%
  group_by(cluster_kproto, variable) %>%
  summarise(n_genre = sum(as.numeric(value)))%>%
  mutate(perc_genre = (n_genre/sum(n_genre))) %>%
  filter(n_genre != 0)
dfSummary_kproto_key



##############
# Second, lets make several plots about the mode and keys (k-proto)
##############


# create a plot of  distribution of mode
mode_cluster_plot <- ggplot(data = dfSummary_kproto_mode, aes(x = cluster_kproto, y = perc_mode, fill = mode))+
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_brewer(palette="Paired", name = "Mode", labels = c("Minor (0)", "Major (1)")) + 
  labs(x = "Cluster", y = "% Of Cluster")+
  scale_y_continuous(labels = function(x) paste0(x*100, "%"))+
  theme_bw()+
  theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16))
mode_cluster_plot


# get percentage of key in all genres
dfKeys <- dfSpotify_clustered %>%
  select(starts_with("key_"))
dfKeys_num <- apply(dfKeys,2, as.numeric)
dfKeys_perc_all <- colSums(dfKeys_num)/nrow(dfKeys_num)

# calculate the diff between cluster key percentage and average
unmelted_dfKeys <- dcast(data = dfSummary_kproto_key,formula = cluster_kproto~variable,fun.aggregate = sum,value.var = "perc_genre")
dfKeys_perc_all_rep <- do.call("rbind", replicate(K_ideal, dfKeys_perc_all, simplify = FALSE))
dfKeys_perc_diff <- melt(cbind(cluster_kproto = unmelted_dfKeys[,1], unmelted_dfKeys[,-1]- dfKeys_perc_all_rep),id="cluster_kproto")

# change the keys in order to make them more understandable
current_keys <- c(as.character((unique(dfKeys_perc_diff$variable))))
Key_names <- c("C", "C/D", "D", "D/E", "E", "F","F/G","G", "G/A", "A", "A/B", "B")
new_keys <- mapvalues(dfKeys_perc_diff$variable, from = current_keys, to = Key_names)
dfKeys_perc_diff$variable <- new_keys

# plot the keys 
key_cluster_plot <- ggplot(data = dfKeys_perc_diff, aes(x = cluster_kproto, y = value, fill = variable))+
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_brewer(palette="Set3", name = "Key") + 
  labs(x = "Cluster", y = "% Of Cluster, Difference from overall % ")+
  scale_y_continuous(labels = function(x) paste0(x*100, "%"))+
  theme_bw()+
  theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16))
key_cluster_plot

##################################
# Third, lets make some graphs about average values and genre for both k-proto and kmeans
##################################

## get mean of all numeric variables in whole dataset
dfSummary_mean_all <- dfSpotify_clustered %>%
  select(-starts_with("key_"), - mode, -cluster_kmeans) %>%
  summarise_all(mean)

# get the diff between mean in clusters and overall mean
dfSummary_mean_all_rep <- do.call("rbind", replicate(K_ideal, dfSummary_mean_all[-10], simplify = FALSE))
dfSummary_kproto_diff <- cbind(dfSummary_kproto[,1], dfSummary_kproto[,-1]-dfSummary_mean_all_rep)
dfSummary_kmeans_diff <- cbind(dfSummary_kmeans[,1], dfSummary_kmeans[,-1]-dfSummary_mean_all_rep)

  
# melt the datasets in order to visualise the average values
melted_summary_kproto <- melt(dfSummary_kproto_diff, id = "cluster_kproto")
melted_summary_kmeans <- melt(dfSummary_kmeans_diff, id = "cluster_kmeans")

# graph of the average values for the numeric values for k-prototypes
Avg_Numeric_kproto <- ggplot(data = melted_summary_kproto, aes(x = cluster_kproto, y = value, fill = variable))+
  geom_bar(stat="identity", position=position_dodge()) + 
  scale_fill_brewer(palette="Set1", name = "Variable") + 
  labs(x = "Cluster", y = "Difference in cluster from Avg. Value")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16))
Avg_Numeric_kproto + geom_vline(xintercept = 1.5) + geom_vline(xintercept = 2.5) + geom_vline(xintercept = 3.5)


# create df to compare both averages
dfCompare_avg_methods <- cbind(melted_summary_kmeans, melted_summary_kproto)

# and for k-means
Avg_Numeric_kmeans <- ggplot(data = melted_summary_kmeans, aes(x = cluster_kmeans, y = value, fill = variable))+
  geom_bar(stat="identity", position=position_dodge()) + 
  scale_fill_brewer(palette="Set1", name = "Variable") + 
  labs(x = "Cluster", y = "Difference in cluster from Avg. Value")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16))
Avg_Numeric_kmeans  + geom_vline(xintercept = 1.5) + geom_vline(xintercept = 2.5) + geom_vline(xintercept = 3.5)


# create a dataframe with the independent variables, the clusters, and further info (track name, genre)
dfSpotify_complete <- cbind(dfSpotify_clustered, 
                            name = spotify_clean$track_name, 
                            artist = spotify_clean$track_artist,
                            popularity =spotify_clean$track_popularity,
                            genre = spotify_clean$playlist_genre )


# check the average popularity per cluster
dfSummary_popularity <- dfSpotify_complete %>% 
                          group_by(cluster_kproto) %>%
                          summarise(mean_popularity =mean(popularity))
dfSummary_popularity

# summary df of the percentage per genre
dfSummary_genres_kproto <- dfSpotify_complete %>%
  group_by(cluster_kproto, genre) %>%
  summarise(n_genre = n()) %>%
  mutate(perc_genre = (n_genre/sum(n_genre)))
dfSummary_genres_kproto

# summary df of the percentage per genre
dfSummary_genres_kmeans <- dfSpotify_complete %>%
  group_by(cluster_kmeans, genre) %>%
  summarise(n_genre = n()) %>%
  mutate(perc_genre = (n_genre/sum(n_genre)))
dfSummary_genres_kmeans

# distribution of genre per cluster - kproto
Genre_plot_kproto <- ggplot(data = dfSummary_genres_kproto[,-3], aes(x = cluster_kproto, y = perc_genre, fill = genre))+
  geom_bar(stat="identity", position=position_dodge()) + 
  scale_fill_brewer(palette="Set2",  name = "Genre") + 
  labs(x = "Cluster", y = "% Of Cluster")+ 
  scale_y_continuous(labels = function(x) paste0(x*100, "%"))+
  theme_bw()+
  theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16))
Genre_plot_kproto + geom_vline(xintercept = 1.5) + geom_vline(xintercept = 2.5) + geom_vline(xintercept = 3.5)

# distribution of genre per cluster - kmeans
Genre_plot_kmeans <- ggplot(data = dfSummary_genres_kmeans[,-3], aes(x = cluster_kmeans, y = perc_genre, fill = genre))+
  geom_bar(stat="identity", position=position_dodge()) + 
  scale_fill_brewer(palette="Set2",  name = "Genre") + 
  labs(x = "Cluster", y = "% Of Cluster")+ 
  scale_y_continuous(labels = function(x) paste0(x*100, "%"))+
  theme_bw()+
  theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16))
Genre_plot_kmeans + geom_vline(xintercept = 1.5) + geom_vline(xintercept = 2.5) + geom_vline(xintercept = 3.5)


# histogram of the continuous variables
dfX_spotify_melted <- melt(data.frame(mX_spotify))

##################################
# Check distribution per variable
##################################

Hist_Var_Plot <- ggplot(data = dfX_spotify_melted, aes(x = value))+
  geom_histogram(stat='bin', colour="black", fill="#3090C7")+
  facet_wrap(~toupper(variable), scales = "free")+
  labs(x = "Value", y = "Count")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16))
Hist_Var_Plot
