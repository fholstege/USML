################################################################################
# A) Load packages & clean data
################################################################################

# empty global environment
rm(list=ls())

# Packes required for subsequent analysis. P_load ensures these will be installed and loaded. 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(contextual,
               tidyverse,
               rlist,
               reshape2,
               readr,
               data.table,
               rlist,
               parallel,
               flexclust,
               class,
               Rfast,
               devtools) 

# add repository to make mclapply() run in parallel (only necessary on windows)
install_github('nathanvan/parallelsugar')
library(parallelsugar)


# get clean data with arms, rewards, days, and user features
dfCleanData <- read.csv("FinalProject/dfUserData_readyForAnalysis.csv")
dfCleanData[,1]<-NULL

# get dataframe of the first 5 days for the customer segmentation
dfSample_Kmeans <- dfCleanData %>% filter(day <= 5) %>% sample_n(10000000)
dfTestSet <-  dfCleanData %>% filter(day > 5)


# Function to compute total within-cluster sum of square 
wss <- function(k) {
  result <- kmeans(dfSample_Kmeans %>% select(-day, -arm, -reward), k)$tot.withinss
  return(result)
}

# Compute and plot wss for k = 1 to k = 15
k.values <- 1:5
tot.withinss_results <- map_dbl(k.values, wss)

# plot 
plot(k.values, tot.withinss_results,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")


# create model to determine in which cluster observations belong 
model_kmeans <- kcca(dfSample_Kmeans %>% select(-day, -arm, -reward), k=3, family = kccaFamily("kmeans"))

# assign to which cluster the customers belong in the testset
clusters_TestSet <- predict(model_kmeans, dfTestSet %>% select(-day, -arm, -reward))

dfTestSet$cluster <- clusters_TestSet
dfTestSet$user_vist_id <- 1:10000000

getwd()
write.csv(dfTestSet, "FinalProject/dfTestSet.csv")







