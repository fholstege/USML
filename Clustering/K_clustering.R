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
#           F) Implement K-means and K-proto algorithms - find optimal K and analyse clusters
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
               plyr) 
# load the data
spotify_songs <- read_csv('https://raw.githubusercontent.com/rfordatascience/tidytuesday/master/data/2020/2020-01-21/spotify_songs.csv')

# set seed
set.seed(123)

################################################################################
# Data prep functions
################################################################################

# Clean the data
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
  
  result = list(centroids = fin_centroids, tot_within_SS = curr_tot_within_SS, within_SS = fin_within_SS, clusters = fin_clusters)
  
  return(result)
}

###############################################################################
# D) Compare to package
##############################################################################

# show for K = 3
K = 3

load("Data/cityweather.Rdata")
mX_city <- as.matrix(cityweather)

# check our results 
result_ours <- K_means(mX_city, K,n_iter = 10, n_random_centroids = 10)
result_ours$centroids
result_ours$within_SS

# check results package
result_pack <- kmeans(mX_city, K, iter.max=10, nstart=10)
result_pack$centers
result_pack$withinss

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

vKey <- data.selection(spotify_clean , 'key')
dummy_key <- dummy_cols(vKey)[,-1]
dummy_key_factor <- data.frame(lapply(dummy_key, factor))

spotify_selected <- data.selection(spotify_clean, select.cols.con)
spotify_selected_cat <- cbind(spotify_selected, mode = spotify_clean$mode,dummy_key_factor) 


# Adapt variables that need scaling
select.scale.vars <- c('loudness', 'tempo')
spotify_final <- data.scaling(spotify_selected, method = "minmax", select.scale.vars)
spotify_final_cat <- data.scaling(spotify_selected_cat, method = "minmax", select.scale.vars)

# turn to matrix
mX_spotify <- as.matrix(spotify_final)

# set mode and key as factors for the kproto
spotify_final_cat$mode <- as.factor(spotify_final_cat$mode)




###############################################################################
# E) Determine optimal K
##############################################################################


######
# find_K_elbow: generate plot with within ss per K 
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
    labs(x = "K", y = "Total Within Sum of Squared Distances")+
    theme_bw() + 
    theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 16),
          axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16)) 
  
  return(elbowPlot)
  
  
} 

##########################
# Find ideal Lambda and K
#########################


# check the average standard deviation
sigma_col <- apply(mX_spotify, 2, sd)
avg_sigma <- mean(sigma_col)

# paper says - suitable gamma is between 1/3*sigma and 2/3* sigma
lambda = 0.5 * avg_sigma

# range of values for K - check between 2 and 10
K_range_spotify <- 2:10
elbowPlot_kmeans <- find_K_elbow(mX_spotify, K_range = K_range_spotify, n_iter = 100, n_random_centroids = 10, type='kmeans')
elbowPlot_kproto <- find_K_elbow(spotify_final_cat, K_range = K_range_spotify, n_iter = 100, n_random_centroids = 10, type='kproto', lambda = lambda)
elbowPlot_kproto 
elbowPlot_kproto + geom_vline(xintercept = 4, linetype="dotted", color = "black", size=1.5) + 



###############################################################################
# E) Analyse with optimal K
##############################################################################

# ideal K 
K_ideal <- 4
result_kmeans <- kmeans(mX_spotify, K_ideal, iter.max = 100, nstart = 10)
result_kproto <- kproto(spotify_final_cat, K_ideal, iter.max = 100, nstart = 10,lambda = lambda)





# add the dataframes of the results of the cluster analysis and the variables
dfSpotify_clustered <- data.frame(cbind(spotify_final_cat, 
                                  cluster_kmeans = result_kmeans$cluster), 
                                  cluster_kproto = result_kproto$cluster )


# summarise the average values for k-means
dfSummary_kmeans <- dfSpotify_clustered %>% 
  group_by(cluster_kmeans) %>%
  summarise_all(mean) %>%
select(-starts_with("key_"), - mode, -cluster_kproto)
dfSummary_kmeans

# summarise the average values for the k-prototype
dfSummary_kproto <- dfSpotify_clustered %>% 
  group_by(cluster_kproto) %>%
  summarise_all(mean) %>%
  select(-starts_with("key_"), - mode, -cluster_kmeans)
dfSummary_kproto

dfSummary_mean_all <- dfSpotify_clustered %>%
select(-starts_with("key_"), - mode, -cluster_kmeans) %>%
summarise_all(mean)

dfSummary_mean_all_rep <- do.call("rbind", replicate(K_ideal, dfSummary_mean_all[-10], simplify = FALSE))
dfSummary_kproto_diff <- cbind(dfSummary_kproto[,1], dfSummary_kproto[,-1]-dfSummary_mean_all_rep)


# summarise the percentage for the mode variable per cluster
dfSummary_kproto_mode <- dfSpotify_clustered %>% 
  group_by(cluster_kproto, mode) %>%
  summarise(n_mode = n()) %>%
  mutate(perc_mode = (n_mode/sum(n_mode)))
dfSummary_kproto_mode

# plot the distribution of mode
mode_cluster_plot <- ggplot(data = dfSummary_kproto_mode, aes(x = cluster_kproto, y = perc_mode, fill = mode))+
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_brewer(palette="Paired", name = "Mode", labels = c("Minor (0)", "Major (1)")) + 
  labs(x = "Cluster", y = "% Of Cluster")+
  scale_y_continuous(labels = function(x) paste0(x*100, "%"))+
  theme_bw()+
  theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16))
mode_cluster_plot

# summarise the percentage for the key variable(s) per cluster
dfSummary_kproto_key <- dfSpotify_clustered %>% 
  select(starts_with("key_"), cluster_kproto) %>%
  melt(id = 'cluster_kproto') %>%
  group_by(cluster_kproto, variable) %>%
  summarise(n_genre = sum(as.numeric(value)))%>%
  mutate(perc_genre = (n_genre/sum(n_genre))) %>%
  filter(n_genre != 0)
dfSummary_kproto_key

# get percentage of key in all genres
dfKeys <- dfSpotify_clustered %>%
  select(starts_with("key_"))
dfKeys_num <- apply(dfKeys,2, as.numeric)
dfKeys_perc_all <- colSums(dfKeys_num)/nrow(dfKeys_num)

unmelted_dfKeys <- dcast(data = dfSummary_kproto_key,formula = cluster_kproto~variable,fun.aggregate = sum,value.var = "perc_genre")
unmelted_dfKeys

dfKeys_perc_all_rep <- do.call("rbind", replicate(K_ideal, dfKeys_perc_all, simplify = FALSE))
dfKeys_perc_diff <- melt(cbind(cluster_kproto = unmelted_dfKeys[,1], unmelted_dfKeys[,-1]- dfKeys_perc_all_rep),id="cluster_kproto")
dfKeys_perc_diff

current_keys <- c(as.character((unique(dfKeys_perc_diff$variable))))
Key_names <- c("C", "C/D", "D", "D/E", "E", "F","F/G","G", "G/A", "A", "A/B", "B")

new_keys <- mapvalues(dfKeys_perc_diff$variable, from = current_keys, to = Key_names)
dfKeys_perc_diff$variable <- new_keys

dfKeys_perc_diff
key_cluster_plot <- ggplot(data = dfKeys_perc_diff, aes(x = cluster_kproto, y = value, fill = variable))+
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_brewer(palette="Set3", name = "Key") + 
  labs(x = "Cluster", y = "% Of Cluster, Difference from overall % ")+
  scale_y_continuous(labels = function(x) paste0(x*100, "%"))+
  theme_bw()+
  theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16))
key_cluster_plot
  
# melt the datasets in order to visualise the average values
melted_summary_kproto <- melt(dfSummary_kproto_diff, id = "cluster_kproto")


# graph of the average values for the numeric values for k-prototypes
Avg_Numeric_kproto <- ggplot(data = melted_summary_kproto, aes(x = cluster_kproto, y = value, fill = variable))+
  geom_bar(stat="identity", position=position_dodge()) + 
  scale_fill_brewer(palette="Set1", name = "Variable") + 
  labs(x = "Cluster", y = "Difference in cluster from Avg. Value")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16))
Avg_Numeric_kproto

# create a dataframe with the independent variables, the clusters, and further info (track name, genre)
dfSpotify_complete <- cbind(dfSpotify_clustered, 
                            name = spotify_clean$track_name, 
                            artist = spotify_clean$track_artist,
                            popularity =spotify_clean$track_popularity,
                            genre = spotify_clean$playlist_genre )

# check the average popularity per cluster
dfSummary_popularity <- dfSpotify_complete %>% 
  group_by(cluster_kproto) %>%
  summarise(avg_popularity = mean(popularity))
dfSummary_popularity


# summary df of the percentage per genre
dfSummary_genres <- dfSpotify_complete %>%
  group_by(cluster_kproto, genre) %>%
  summarise(n_genre = n()) %>%
  mutate(perc_genre = (n_genre/sum(n_genre)))
dfSummary_genres

# distribution of genre per cluster
Genre_plot <- ggplot(data = dfSummary_genres[,-3], aes(x = cluster_kproto, y = perc_genre, fill = genre))+
  geom_bar(stat="identity", position=position_dodge()) + 
  scale_fill_brewer(palette="Set2",  name = "Genre") + 
  labs(x = "Cluster", y = "% Of Cluster")+ 
  scale_y_continuous(labels = function(x) paste0(x*100, "%"))+
  theme_bw()+
  theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16))
Genre_plot

# make ggpairs object
ggpairs(dfSpotify_clustered, mapping = aes(colour = factor(dfSpotify_clustered$cluster_kproto)))


dists_kproto <- result_kproto$dists
index_min_dist_1 <- which.min(result_kproto$dists[,1])
index_min_dist_2 <- which.min(result_kproto$dists[,2])
index_min_dist_3 <- which.min(result_kproto$dists[,3])
index_min_dist_4 <- which.min(result_kproto$dists[,4])

dfSpotify_complete[c(23982),]
