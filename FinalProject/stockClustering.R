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
               rlist) 

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


# based on plot, cut in number of clusters
n_clusters <- 10

# create parameters dendogram visualisation
cluster_results_dendogram <- as.dendrogram(cluster_results_past5Years)
nodePar <- list(lab.cex = 0.6, pch = c(NA, 19), 
                cex = 0.7, col = "black")
colors_dendogram <-brewer.pal(n_clusters, "Spectral")

# create dendogram visualisation 
cluster_results_dendogram %>% set("branches_k_color", 
    value = colors_dendogram, k = n_clusters) %>% 
  plot(xlab = "Avg. Distance Between Clusters",
       nodePar = nodePar, horiz = TRUE) +
  abline(v = 1.084)


## cut the tree into ten clusters, allocate each stock to a cluster
clusters_stocks <- cutree(cluster_results, k=n_clusters)
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
  labs(fill = "GISC Sector", x = "Cluster", y = "Count")


#######################################
# C) Backtest diversification benefits from clusters
#######################################

calc_portfolio_var <- function(w0, mCovar){
  
  n <- nrow(mCovar)
  vWeights <- rep(w0/n, n)
  
  portfolio_var <- t(vWeights) %*% mCovar %*% vWeights
  return(portfolio_var)
  
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
  portfolio_var_clusters <- calc_portfolio_var(w0, mCovar_clusters)
  
  # calculate the average volatility of the assets
  dfClusterVol <- dfClusterReturns %>%
    group_by(cluster) %>%
    summarise(vol = sd(avg_returns))
  Avg_vol_clusters <- mean(dfClusterVol$vol)*w0
  
  # calculate the geometric cumulative returns
  vReturns <-  apply(dfClusterReturns_transposed[,-1], 2,Return.cumulative)
  
  list_results <- list(portfolio_var = portfolio_var_clusters, avg_vol = Avg_vol_clusters, cumulative_returns =vReturns , all_returns =dfClusterReturns_transposed )
  
  return(list_results)
}

backtest_clustering <- function(index_date, years_for_recalc,years_for_correlation,dfReturns_all ,n_clusters,sectors=FALSE,dfSectorInfo=NA,min_date = "2000-01-01"){
  
  # check number of years inbetween to get number of periods that can be used for backtest
  timeDiff <- (as.Date(index_date) - as.Date(min_date))
  n_days_diff <- as.integer(timeDiff)
  n_periods_backtest <- round(n_days_diff/(365*years_for_recalc) ,0)
  
  # account for number of years of backhistory that need to be available for the backtest
  vPotential_Index_dates <- as.Date(index_date) - seq(1:n_periods_backtest) * (365 * 1)
  vIndex_dates <- head(vPotential_Index_dates, n=-years_for_correlation)
  print(vIndex_dates)

  lResults <- lapply(as.character(vIndex_dates), function(date){
    
    # get clusters for the previous x periods
    clusterObj <- calc_clusters(index_date=as.Date(date), years_for_correlation=years_for_correlation, dfReturns_all)
    
    # cut the tree into ten clusters, allocate each stock to a cluster
    if(!sectors){
      clusters_stocks <- cutree(clusterObj, k=n_clusters)
      dfStockClusters <- data.frame(Ticker = names(clusters_stocks),cluster = clusters_stocks)
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
  dfCumulativeReturns <- cbind(index_date = vIndex_dates, list.rbind(lapply(lResults, function(x){return(x$cumulative_returns)})))
  
  
  dfDiversification_results <- data.frame(index_date =vIndex_dates,  portfolio_var = portfolio_vars, avg_vol_portfolioItem = avg_vol)
  
  lReturn <- list(full_results = lResults, 
                  dfDiversification = dfDiversification_results, 
                  dfCumulativeReturns = dfCumulativeReturns,
                  rolling_correlation = years_for_correlation, 
                  duration_before_switch = years_for_recalc,
                  n_clusters = n_clusters)
    
    
    
  return(lReturn )
}


index_date = "2021-12-31"
min_date = "2000-01-01"

colnames(dfStockInfo)[3] <- "cluster"

backtest_result_clusters <- backtest_clustering(index_date, years_for_recalc =1, years_for_correlation = 3, dfReturns_all = dfStockReturns_after2000, n_clusters = 10) 
backtest_result_sectors <- backtest_clustering(index_date, years_for_recalc =1, years_for_correlation = 3, dfReturns_all = dfStockReturns_after2000, n_clusters = 7,sectors=TRUE,dfSectorInfo=dfStockInfo) 


dfDiversification_clusters <- backtest_result_clusters$dfDiversification
dfDiversification_clusters$Diversification_benefit <- -(dfDiversification_clusters$portfolio_var - dfDiversification_clusters$avg_vol_portfolioItem)/dfDiversification_clusters$avg_vol_portfolioItem

dfDiversification_sectors <- backtest_result_sectors$dfDiversification
dfDiversification_sectors$Diversification_benefit <- -(dfDiversification_sectors$portfolio_var - dfDiversification_sectors$avg_vol_portfolioItem)/dfDiversification_sectors$avg_vol_portfolioItem


dfCompareDiversification_benefit <- data.frame(date = strftime(dfDiversification_clusters$index_date[c(-1, -2)], format = "%Y"), 
                                               clustering = dfDiversification_clusters$Diversification_benefit[c(-1, -2)],
                                               sectors = dfDiversification_sectors$Diversification_benefit[c(-1, -2)])


dfCompareDiversification_benefit_melted <- melt(dfCompareDiversification_benefit, id.vars = 'date')


ggplot(data = dfCompareDiversification_benefit_melted, aes(x = as.factor(date), y = value, fill = variable))+
  geom_bar(stat="identity", color="black", position=position_dodge())+
  labs(x = "Date")
