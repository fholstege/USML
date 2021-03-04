######################
#
#
#
######################

# Packes required for subsequent analysis. P_load ensures these will be installed and loaded. 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(readxl,
               quantmod,
               dendextend,
               lubridate,
               reshape2,
               RColorBrewer,
               PerformanceAnalytics,
               rlist,
               tidyverse) 

# load the prices for all stocks in the NASDAQ index
coltypes <- c("date", rep("numeric", 102))
dfStockPrices <- read_xlsx("Data/stockprices_nasdaq.xlsx", sheet = 1, col_names = TRUE, col_types = coltypes, na = "", skip = 0)

# load the info for all the stocks in the NASDAQ index
dfStockInfo <- read_xlsx("Data/nasdaq_current_index_info.xlsx", sheet = 1)
colnames(dfStockInfo)[3:4]<- c("GISC_Sector", "GISC_SubSector")



#######################################
# A) Clean the data for subsequent analysis
#######################################

# calculate arithmetic returns for each stock
dfStockReturns <- data.frame(apply(dfStockPrices[,-1], 2, function(prices){
  exp(diff(log(prices))) - 1
}) )

# get for each stock the date at which the data is at first available
first_available_date_index <- apply(dfStockPrices,2,function(x){min(which(!is.na(x)))})
first_available_dates <- dfStockPrices$Date[first_available_date_index]

# select those stocks for which there is data available after 2000
index_available_after_2000 <- first_available_dates < "2000-01-01"

# add the date column, and define a 
dfStockReturns <- cbind(dfStockPrices$Date[-1], dfStockReturns)
colnames(dfStockReturns)[1] <- "Date"

# define dataframe with only stocks available after 2000, get returns of the last 5 years
dfStockReturns_after2000 <- dfStockReturns[,index_available_after_2000]


#######################################
# B) Functions and analysis to illustrate the (present) clusters
#######################################


##
calc_clusters <- function(index_date, years_for_correlation, dfReturns, method = "average"){
  
  # get the returns from the index date until X years back 
  min_date_backhistory <- index_date - (years_for_correlation  * 365)
  dfStockReturns_backHistory <- dfReturns %>% filter(Date >= min_date_backhistory & Date <= index_date)

  
  # calc correlations, and turn correlations to euclidean distances
  dfStockReturns_correlationMatrix <- cor(dfStockReturns_backHistory %>% select(-Date))
  
  mDist_stockCorrelations <- sqrt(2 * (1- dfStockReturns_correlationMatrix))
  
  # turn the matrix to a distance object for hierarchical clustering
  diag(mDist_stockCorrelations) <- 0
  mDist_stockCorrelations[upper.tri(mDist_stockCorrelations)] <- NA
  mDist_stockCorrelations_distObj <- as.dist(mDist_stockCorrelations, diag = TRUE)
  

  # get results of hierarchical clustering - using average distance
  cluster_results <- hclust(mDist_stockCorrelations_distObj, method = method)
  
  return(cluster_results)
  
}


# get the returns from the index date until X years back 
cluster_results_past5Years <- calc_clusters(as.Date("2020-12-31"), 5, dfStockReturns_after2000)

## cut the tree into  clusters, allocate each stock to a cluster
clusters_stocks <- cutree(cluster_results_past5Years, k=n_clusters)

dfStockClusters <- data.frame(Ticker = names(clusters_stocks),cluster = clusters_stocks)
dfStockInfo_clusters <- merge(dfStockInfo, dfStockClusters, on = "Ticker")


# check which stocks are in each cluster
dfStockClusters_sum <- dfStockInfo_clusters %>% 
    group_by(cluster, GISC_Sector) %>%
    summarise(count = n())

# make plot to illustrate the difference in sector division across clusters
ggplot(data = dfStockClusters_sum, aes(x = as.factor(cluster), y = count, fill = GISC_Sector))+
  geom_bar(stat = 'identity') +
  theme_classic() + 
  labs(fill = "GISC Sector", x = "Cluster", y = "Count") + 
  theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16))


# create parameters dendogram visualisation
cluster_results_dendogram <- as.dendrogram(cluster_results_past5Years)
nodePar <- list(lab.cex = 0.6, pch = c(NA, 19), 
                cex = 0.7, col = "black")
colors_dendogram <-brewer.pal(n_clusters, "Set1")


# based on plot, cut in number of clusters, set to company names
n_clusters <- sum(cluster_results_past5Years$height>1.1)
idx <- sapply(cluster_results_past5Years$labels, function(x) {
  which(dfStockInfo$Ticker == x)
})

cluster_results_past5Years$labels <- dfStockInfo[idx,]$Company
cluster_results_past5Years$labels <-sapply(strsplit(cluster_results_past5Years$labels," "), `[`, 1)


# create dendogram visualisation 
cluster_results_dendogram %>% set("branches_k_color", 
                                  value = colors_dendogram, k = n_clusters, branches_lwd = 3) %>% 
  plot(xlab = "Avg. Distance Between Clusters",
       nodePar = nodePar, horiz = TRUE) %>% 
  abline(v = 1.1, lty = 2) 





#######################################
# C) Backtest diversification benefits from clusters
#######################################

calc_portfolio_vol <- function(w0, mCovar){
  
  n <- nrow(mCovar)
  vWeights <- rep(w0/n, n)
  
  portfolio_var <- t(vWeights) %*% mCovar %*% vWeights
  
  
  return(sqrt(portfolio_var))
  
}


calc_portfolioStats_period_clusters <- function(dfReturns, dfClusterInfo, w0 = 100){
  
  dfReturns_melted <- melt(dfReturns, id.vars = 'Date')
  colnames(dfReturns_melted)[2] <- "Ticker"
  
  dfReturns_withClusters  <- merge(dfReturns_melted, dfClusterInfo, on = "Ticker")
  
  dfClusterReturns <- dfReturns_withClusters %>%
    group_by(Date, cluster) %>%
    summarise(avg_returns = mean(value))
  
  
  # calculate the portfolio variance, accounting for correlation
  dfClusterReturns_transposed <- dcast(dfClusterReturns, Date~cluster, value = 'avg_returns')
  mCovar_clusters <- cov(dfClusterReturns_transposed %>% select(-Date))
  portfolio_var_clusters <- calc_portfolio_vol(w0, mCovar_clusters)
  
  # calculate the average volatility of the assets
  dfClusterVol <- dfClusterReturns %>%
    group_by(cluster) %>%
    summarise(vol = sd(avg_returns))
  Avg_vol_clusters <- mean(dfClusterVol$vol)*w0
  

  #calculate the average correlation among the assets
  mCor <- cor(dfClusterReturns_transposed %>% select(-Date))
  mCor[upper.tri(mCor, diag = TRUE)]<- NA
  avg_cor <- mean(mCor, na.rm = TRUE)
  
  list_results <- list(portfolio_var = portfolio_var_clusters, avg_vol = Avg_vol_clusters,avg_cor =avg_cor, all_returns =dfClusterReturns_transposed )
  
  return(list_results)
}

backtest_clustering <- function(index_date, years_for_recalc,years_for_correlation,dfReturns_all ,cutoff_distance_clusters=1,sectors=FALSE,dfSectorInfo=NA,min_date = "2000-01-01", threshold_n_clusters = 0){
  
  # check number of years inbetween to get number of periods that can be used for backtest
  timeDiff <- (as.Date(index_date) - as.Date(min_date))
  n_days_diff <- as.integer(timeDiff)
  n_periods_backtest <- round(n_days_diff/(365*years_for_recalc) ,0)
  
  # account for number of years of backhistory that need to be available for the backtest
  vPotential_Index_dates <- as.Date(index_date) - seq(1:n_periods_backtest) * (365 * years_for_recalc)
  vIndex_dates <- head(vPotential_Index_dates, n=-years_for_correlation/years_for_recalc)
  
  lResults <- lapply(as.character(vIndex_dates), function(date){
    
    # get clusters for the previous x periods
    clusterObj <- calc_clusters(index_date=as.Date(date), years_for_correlation=years_for_correlation, dfReturns_all)
    
    # cut the tree into ten clusters, allocate each stock to a cluster
    if(!sectors){
      
      # based on plot, cut in number of clusters
      n_clusters <- sum(clusterObj$height>cutoff_distance_clusters)
      
      clusters_stocks <- cutree(clusterObj, k=n_clusters)
      dfStockClusters <- data.frame(Ticker = names(clusters_stocks),cluster = clusters_stocks)

      # get min of companies in cluster
      meet_threshold <- c(which(table(clusters_stocks)>threshold_n_clusters))
      dfStockClusters <- dfStockClusters %>% filter(cluster %in% meet_threshold)
      
      # if the boolean sectors = TRUE, then use the sectors as clusters
    }else{
      dfStockClusters <- dfSectorInfo
    }
   
    # using the returns from the subsequent period, determine what happened to the clusters
    dfReturns_nextPeriod <- dfReturns_all %>% filter(Date >= as.Date(date) & Date <= as.Date(date) + (years_for_recalc*365))
    portfolio_stats_clusters <- calc_portfolioStats_period_clusters(dfReturns_nextPeriod, dfStockClusters)
    
    return(portfolio_stats_clusters)
    
  })
  
  portfolio_vars <- unlist(lapply(lResults, function(x){return(x$portfolio_var)}))
  avg_vol <- unlist(lapply(lResults, function(x){return(x$avg_vol)}))
  avg_cor <- unlist(lapply(lResults, function(x){return(x$avg_cor)}))
  
  dfDiversification_results <- data.frame(index_date =vIndex_dates,  portfolio_var = portfolio_vars, avg_vol_portfolioItem = avg_vol)
  
  avg_diversification_benefit <- mean(-((portfolio_vars - avg_vol)/avg_vol))
  
  lReturn <- list(full_results = lResults,
                  avg_diversification_benefit = avg_diversification_benefit,
                  avg_vol = avg_vol,
                  avg_cor = avg_cor,
                  dfDiversification = dfDiversification_results, 
                  rolling_correlation = years_for_correlation, 
                  duration_before_switch = years_for_recalc)
    
    
    
    
  return(lReturn )
}


index_date = "2020-12-31"
min_date = "2000-01-01"

colnames(dfStockInfo)[3] <- "cluster"


backtest_result_clusters <- backtest_clustering(index_date, 
                                                years_for_recalc =1, 
                                                years_for_correlation = 3, 
                                                dfReturns_all = dfStockReturns_after2000, 
                                                cutoff_distance_clusters=1.1, 
                                                threshold_n_clusters = 1) 
backtest_result_sectors <- backtest_clustering(index_date, years_for_recalc =1, years_for_correlation = 5, dfReturns_all = dfStockReturns_after2000,sectors=TRUE,dfSectorInfo=dfStockInfo) 


years_for_correlation <- c(1,2,3,4,5)
cutoff_distance_clusters <- c(0.95,1,1.05, 1.1)
years_for_recalc = c(1)

param_grid <- expand.grid(years_for_correlation, cutoff_distance_clusters, years_for_recalc)
colnames(param_grid) <- c("years_for_corr", "cutoff_distance_clusters", "years_for_recalc")

backtest_param <- apply(param_grid, 1, function(param){
  
  print(param)
  
  result_param <- backtest_clustering(index_date, 
                                      years_for_recalc =param[3],
                                      years_for_correlation =param[1], 
                                      dfReturns_all = dfStockReturns_after2000, 
                                      cutoff_distance_clusters=param[2], 
                                      threshold_n_clusters = 1) 
  
  avg_diversification_benefit <- result_param$avg_diversification_benefit
  avg_cor <- result_param$avg_cor
  
  row_result <- list(years_for_correlation= param[1], years_for_recalc =param[3], cutoff_distance_clusters=param[2],avg_div_benefit = avg_diversification_benefit, avg_cor = mean(avg_cor, na.rm = TRUE))

  return(row_result)
})

backtest_param_results <- t(matrix(unlist(backtest_param), nrow = 5))
dfBacktest_param_results <- data.frame(backtest_param_results)
colnames(dfBacktest_param_results) <- c("years_for_correlation", "years_for_recalc", "cutoff_distance_clusters", "avg_diversification_benefit", "avg_correlation_clusters")
dfBacktest_param_results


dfDiversification_clusters <- backtest_result_clusters$dfDiversification
dfDiversification_clusters$Diversification_benefit <- -(dfDiversification_clusters$portfolio_var - dfDiversification_clusters$avg_vol_portfolioItem)/dfDiversification_clusters$avg_vol_portfolioItem

dfDiversification_sectors <- backtest_result_sectors$dfDiversification
dfDiversification_sectors$Diversification_benefit <- -(dfDiversification_sectors$portfolio_var - dfDiversification_sectors$avg_vol_portfolioItem)/dfDiversification_sectors$avg_vol_portfolioItem


dfCompareDiversification_benefit <- data.frame(date = year(dfDiversification_clusters$index_date), 
                                               clustering = dfDiversification_clusters$Diversification_benefit,
                                               sectors = dfDiversification_sectors$Diversification_benefit)

dfCompareDiversification_benefit_melted <- melt(dfCompareDiversification_benefit, id.vars = 'date')
dfCompareDiversification_benefit_melted %>% 
  group_by(variable)%>%
  summarise(avg_reduction_vol = mean(value))

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

ggplot(data = dfCompareDiversification_benefit_melted, aes(x = date, y = value, fill = firstup(as.character(variable)))) + 
  geom_bar(stat = 'identity', position = position_dodge())+
  theme_bw() + 
  scale_fill_manual(values = c("#00ccff", "#ff0066"))+
  labs(x = "Year", y = "Reduction In Volatility From Diversification", fill = "Method")+
  scale_y_continuous(labels = function(x) paste0(x*100, "%"))+
       theme(axis.text.x = element_text(size = 13), axis.title.x = element_text(size = 14),
                             axis.text.y = element_text(size = 13), axis.title.y = element_text(size = 14))

mean(backtest_result_clusters$avg_cor, na.rm= TRUE)
mean(backtest_result_sectors$avg_cor)
