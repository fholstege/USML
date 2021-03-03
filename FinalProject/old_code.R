

sim_agents_offlineReplay <- function(df,formula, n_sim, size_sim){
  
  bandit <- OfflineBootstrappedReplayBandit$new(formula = formula, data = df, replacement = TRUE, randomize = FALSE)
  
  agents <- list(Agent$new(LinUCBDisjointOptimizedPolicy$new(alpha = 0.25), bandit, "Lin UCB - Alpha=0.25")) 
  
  simulator <- Simulator$new(agents, horizon= size_sim, do_parallel = TRUE, simulations = n_sim)
  
  history <- simulator$run()
  
  return(history)
}

sim_agents_offlineReplay_allDays <- function(dfAllDays, formula, n_sim, size_sim){
  
  
  lDayData <- split(dfAllDays, f = dfAllDays$day)
  
  lResults_day <- lapply(lDayData, function(dfDay){
    
    dfDay$day <- NULL
    
    
    result_sim <- sim_agents_offlineReplay(dfDay, formula = formula, n_sim = n_sim, size_sim = size_sim)
    
    return(result_sim)
  })
  
}






# set up formula 
f <- formula("reward ~ arm | user1 + user2 + user3 + user4 + user5")


dfSmallSim <- dfCleanData %>%  group_by(day) %>% sample_n(1000) %>% filter(day ==1 )


lResults <- sim_agents_offlineReplay_allDays(dfSmallSim, f, n_sim = 1, size_sim = 200)
hist1 <- lResults[[1]]
hist1$data$t


# And the cumulative reward rate, which equals the Click Through Rate)
plot(lResults[[1]], type = "cumulative", regret = FALSE, rate = TRUE)



