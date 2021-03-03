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

# get clean data with arms, rewards, days, and user features
dfCleanData <- read.csv("FinalProject/dfUserData_readyForAnalysis.csv")
dfCleanData[,1]<-NULL




################################################################################
# B) Functions for running contextual bandit in package
################################################################################


##########################
# sim_agents_offlineReplay: simulates a list of agent for a given dataframe
#

sim_agents_offlineReplay <- function(df,formula, n_sim, size_sim, n_arms, agent_type,param,...){
  
  # get the arms in df, and sample a random n amount
  arms_in_df <- unique(df$arm)
  selected_arms <- sample(arms_in_df, n_arms)
  
  # leave an df with only these arms
  df_selected <- df %>% filter(arm %in% selected_arms)
  
  # notify user if smaller amount of data than wanted because of the arm restriction
  if(nrow(df_selected) < size_sim){
    print("The data does not have enough observations - reducing the size of the simulation")
    size_sim <- nrow(df_selected)
  }
  
  # define the bandit
  bandit <- OfflineBootstrappedReplayBandit$new(formula = formula, data = df_selected, replacement = TRUE, randomize = FALSE)
  print(param)
  # 
  if(agent_type == "LinUCB"){
    agents <- list(Agent$new(LinUCBDisjointOptimizedPolicy$new(alpha = param), bandit,"LinUCB"))
  }else if (agent_type == "LinTS"){
    agents <- list(Agent$new(ContextualLinTSPolicy$new(v = param), bandit, "LinTS"))
  }
  
  
  simulator <- Simulator$new(agents, horizon= size_sim, do_parallel = TRUE, simulations = n_sim, progress_file = TRUE)
  
  history <- simulator$run()
  
  return(history)
}

##########################
# sim_agents_offlineReplay_allDays: simulates a list of agent for a set of days, equally spreading the sample across days
#
sim_agents_offlineReplay_allDays <- function(dfAllDays, formula, n_sim, size_sim, n_arms, agent_type,param){
  
  lDayData <- split(dfAllDays, f = dfAllDays$day)
  
  lResults_day <- lapply(lDayData, function(dfDay){
    
    # indicate to user whats happening
    print(paste0("Currently simulation for day: " ,unique(dfDay$day)))
    
    # remove day variable
    dfDay$day <- NULL
    
    result_sim <- sim_agents_offlineReplay(dfDay, formula = formula, n_sim = n_sim, size_sim = size_sim, n_arms=n_arms, agent_type = agent_type,param = param)
    
    return(result_sim)
  })
  
  return(lResults_day)
  
}

# define this dataframe as the one used in simulations, with enough data points per day
dfSims <- dfCleanData %>% group_by(day) %>% sample_n(10000)


formula_sim =  formula("reward ~ arm | user1 + user2 + user3 + user4 + user5")
lresult_linUCB_alpha_025 <- sim_agents_offlineReplay_allDays(dfAllDays = dfSims %>% filter(day %in% c(1,2)),formula=formula_sim, n_sim = 10, size_sim = 100, n_arms = 10, agent_type = "LinUCB", param = 0.25)
lresult_linUCB_alpha_05 <- sim_agents_offlineReplay_allDays(dfAllDays = dfSims %>% filter(day %in% c(1,2)),formula=formula_sim, n_sim = 10, size_sim = 10000, n_arms = 10, agent_type = "LinUCB", param = 0.5)
lresult_linUCB_alpha_1 <- sim_agents_offlineReplay_allDays(dfAllDays = dfSims %>% filter(day %in% c(1,2)),formula=formula_sim, n_sim = 10, size_sim = 10000, n_arms = 10, agent_type = "LinUCB", param = 1)
lresult_linTS <- sim_agents_offlineReplay_allDays(dfAllDays = dfSims,formula=formula_sim, n_sim = 10, size_sim = 10000, n_arms = 10, agent_type = "LinTS")

agent_names <- c("Lin UCB - Alpha=0.25", "Lin UCB - Alpha=0.5", "Lin UCB - Alpha=1")
lresult_linUCB_alpha_025[[2]]


analyse_results_sim <- function(df_sim_agent, max_obs){
  
  lDayResults <- lapply(lresult_linUCB_alpha_025, function(dfDay){
    
    dfDay_select <- dfDay %>% filter(t <= max_obs) %>% select(t, cum_reward, sim)
    dfDay_select_transposed <- dcast(data = dfDay_select,formula = t~sim,fun.aggregate = sum,value.var = "cum_reward")
    dfDay_select_transposed$t <- NULL
    
  })
  
  dfAllResults <- list.cbind(lDayResults)
  
  n_sims <- ncol(dfAllResults)
  
  SD_CumReward <- apply(dfAllResults, 1, sd)
  avg_CumReward <- apply(dfAllResults, 1, mean)
  lower_bound <- pmax(avg_CumReward - (1.96 * SD_CumReward/sqrt(n_col)),0)
  upper_bound <- avg_CumReward + (1.96 * SD_CumReward/sqrt(n_col))
  
  dfSimResults <- data.frame(t = 1:max_obs,  avg_CumReward = avg_CumReward, sd_CumReward = SD_CumReward, lower_bound = lower_bound, upper_bound = upper_bound)
  
  list_results <- list(dfSimResults =dfSimResults, dfAllResults = cbind(t = 1:max_obs ,dfAllResults) )
  return(list_results)
}



ggplot(data = dfAnalyse_sim, aes(x = t, y = avg_CumReward/t)) + 
  geom_line()

plot(t, type = "cumulative", regret = FALSE, rate = TRUE,
     legend_position = "topleft")




################################################################################
# B) Functions for non-contextual bandit
################################################################################


###
# policy_UCB: picks the arm according to the upper confidence bound criterion 
# 
# Arguments: 
#      dfResults_arms: a dataframe with the following columns
#           Arms: contains the names of the arms used in the simulation                      
#           sample_size: how many times, at that point, the arm has been pulled
#           succes_rate: the % that arm has yielded a good result
#
# Output:
#   chosen arm: index of the arm chosen
#
###
policy_UCB <- function(dfResults_arms, fC){
  
  # get the parameters for UCB - est. the success rate, log of number of steps, and n of times arms puled
  Est_val_arms <- dfResults_arms$succes_rate
  log_t <- log(sum(dfResults_arms$sample_size))
  N_t_arms <- dfResults_arms$sample_size
  
  # calculate the UCB per arm, pick arm with highest
  UCB_score <- Est_val_arms + (fC * (log_t/N_t_arms)^(1/2))
  chosen_arm <- which.max(UCB_score)
  
  return(chosen_arm)
}
###
# policy_TS: picks the arm according to Thompson sampling
# 
# Arguments: 
#      dfResults_arms: a dataframe with the following columns
#           Arms: contains the names of the arms used in the simulation                      
#           sample_size: how many times, at that point, the arm has been pulled
#           succes_rate: the % that arm has yielded a good result
#      n_arms: number of arms
#      eps: % of random arms
#
# Output:
#   chosen arm: index of the arm chosen
#
###

policy_TS <- function(dfResult, n_arms,index_arms ){
  
  # you need to draw from as many distributions as the number of arms. 
  # From each distribution, we draw a 100. Each arm has a different alpha, beta
  n <- rep(100, n_arms)
  alpha <- dfResult$alpha
  beta <- dfResult$beta
  
  # draw observations from distributions using lapply, get the mean of each distribution (sample = 100)
  lObs_result <- lapply(n, rbeta, alpha, beta)
  mean_distribution <- unlist(lapply(lObs_result, mean))
  
  # pick arm with highest sample value
  chosen_arm <- which.max(mean_distribution)
  
  return(chosen_arm)
}


################################################################################
# C) define function to run simulations
################################################################################


###
# simulation_policy: simulates performance of a policy with parameters for both the policy and simulation
# 
#   Arguments: 
#     dfBandit: n x 2 dataframe, with colnames "Arm", "Result" required.
#     policy_type; string, one of "greedy", 
#     exploration: the % of data that first needs to be explored before we start the policy (so for example, time taken before determining best option for greedy)
#     ...: remaining arguments that are necessary for the policy 
#
#   Output: list with following
#     overall; df, success rate per arm
#     after expore: df, results after exploration phase
#     before explore: df, results before exploration phase
#     total: df, results combined (before and after exploration phase)
#
#   
###
sim_policy <- function(dfBandit, policy_type, exploration,size_experiment,only_total=TRUE,...){
  
  ## Part 1: set parameters, number of observations in random exploration and after
  n_obs <- nrow(dfBandit)
  
  if(size_experiment > n_obs){
    print(paste0("The indicated size of the experiment is bigger than the data provided - shrinking to size ", n_obs))
    size_experiment <- n_obs
  }
  n_explore <- round(size_experiment* exploration,0)
  n_exploit <- size_experiment - n_explore
  
  # get random of size
  dfBandit <- dfBandit[sample(1:n_obs, size_experiment), ]
  
  ## part 2: exploration phase 
  # Set the results for the exploration phase by finding a certain number of random observations
  index_exploration <- sample(1:size_experiment,n_explore)
  
  dfResult_explore <- dfBandit[index_exploration,]
  dfBandit_postExplore <- dfBandit[-index_exploration,]
  
  # in exploration phase, find the most succesful arm
  dfResult <- dfResult_explore %>% 
    group_by(Arm)%>%
    summarise(succes_size = sum(result),
              sample_size = n(),
              succes_rate = succes_size/sample_size
    )
  # if policy type is thompson sampling, add alpha and beta param to the dataframe\
  if(policy_type == "TS"){
    start_alpha = 2
    start_beta = 2
    
    dfResult$alpha <- start_alpha + dfResult$succes_size 
    dfResult$beta <- start_beta + (dfResult$sample_size - dfResult$succes_size)
  }
  # get index of best arm
  index_best_arm <- which.max(dfResult$succes_rate)
  
  # get the arms and number of arms 
  Arms <- unique(dfResult_explore$Arm)
  
  # get list of dataframes (one per arm & its results) for remainder
  ldf_arms <- split(dfBandit_postExplore , f = dfResult_explore$Arm )
  index_arms <- 1:length(ldf_arms)
  
  # create vector to save remainder of results
  dfResult_exploit <-  data.frame(matrix(NA, nrow = n_exploit, ncol = 2))
  colnames(dfResult_exploit) <- c("Arm", "result")
  
  ## part 3: appply approach, constantly updating
  i = 1
  while(i <= n_exploit){
    
    if(policy_type == "greedy"){
      
      # pick an arm according to the 'greedy' policy
      chosen_arm <- policy_greedy(index_arms=index_arms, index_best_arm=index_best_arm,...)
      
    }else if (policy_type == "UCB"){
      
      # pick an arm according to the 'UCB' policy
      chosen_arm <- policy_UCB(dfResult, ...)
    } else if (policy_type == "TS"){
      
      # pick an arm according to the thompson sampling policy
      chosen_arm <- policy_TS(dfResult=dfResult, n_arms= length(index_arms), index_arms=index_arms,...)
    }
    
    # get dataframe with remaining results for chosen arm, and pick a random instance
    chosen_arm_df <- ldf_arms[[chosen_arm]]
    
    # obtain the result, remove this from the overall remaining data
    index_result <- sample(1:length(chosen_arm_df$result), 1)
    result <- chosen_arm_df[index_result,]$result
    
    if(length(result) == 0){
      print(i)
      print("You have run out of observations from a chosen arm - decrease the size of the experiment, or increase the size of the dataset")
      index_arms <- index_arms[-chosen_arm]
      break
    }
    
    # get vector of results (arm, result)
    vResult <- c(dfResult[chosen_arm,]$Arm, result)
    dfResult_exploit[i,] <- vResult
    
    # update the overall results data
    dfResult[chosen_arm,]$succes_size <- dfResult[chosen_arm,]$succes_size +  result
    dfResult[chosen_arm,]$sample_size <- dfResult[chosen_arm,]$sample_size + 1
    dfResult[chosen_arm,]$succes_rate <- dfResult[chosen_arm,]$succes_size/dfResult[chosen_arm,]$sample_size
    
    # if policy type is thompson sampling, add alpha and beta param to the dataframe\
    if(policy_type == "TS"){
      dfResult$alpha <- start_alpha + dfResult$succes_size 
      dfResult$beta <- start_beta + (dfResult$sample_size - dfResult$succes_size)
      
    }else if (policy_type == "greedy"){
      # update the best arm
      index_best_arm <- which.max(dfResult$succes_rate)
    }
    # onto the next
    i <- i + 1
  }
  
  # if true, only return the df with all results
  if(only_total){
    results <- list(total = rbind(dfResult_explore,dfResult_exploit))
  }else{
    # save data in list
    results <- list(Overall = dfResult, 
                    after_explore = dfResult_exploit, 
                    before_explore = dfResult_explore, 
                    total = rbind(dfResult_explore,dfResult_exploit),
                    total_succes_rate = sum(dfResult$succes_size)/sum(dfResult$sample_size))
  }
  return(results)
  
}

###
# sim_experiment: In parallel, simulates performance of a policy for an certain number of sims
# 
#   Arguments: 
#     dfBandit: n x 2 dataframe, with colnames "Arm", "Result" required. If day_split = TRUE, then needs to be n x 3 with colnames "Arm", "Result", "day"
#     n_sim: number of simulations
#     n_per_sim: number of observations per sim
#     policy_type: string, one of "greedy", "UCB", "TS
#     exploration: % of the data that should be spent on random exploratoin
#     day_split: boolean, if true, split per day
#     ...: remaining arguments that are necessary for the policy 
#
#   Output: list with following
#     overall; df, success rate per arm
#     after expore: df, results after exploration phase
#     before explore: df, results before exploration phase
#     total: df, results combined (before and after exploration phase)
#
###
sim_experiment <- function(dfBandit, n_sim, n_per_sim, policy_type, exploration,day_split=TRUE,...){
  
  # if true, split per day
  if(day_split){
    # list of dataframes split out per day
    ldf_perDay <- split(dfBandit[,1:2] , f = dfBandit$day )
    
    # repeat enough to get enough observations for sim
    list_dfBandit_exp <- rep(ldf_perDay,n_sim/length(ldf_perDay))
    
  }else{
    # get list of the dataframes - same dataframe, n_sim times
    list_dfBandit_exp <- rep(list(dfBandit),n_sim)
    
  }
  
  
  # get result n_sim times 
  result_per_random_sample <- mclapply(list_dfBandit_exp, 
                                       sim_policy,
                                       policy_type = policy_type, 
                                       exploration = exploration, 
                                       size_experiment=n_per_sim, 
                                       mc.cores = detectCores()-1,...)
  
  return(result_per_random_sample)
}
###
# sim_params: simulates performance of a policy for certain parameters
#
#   Arguments: 
#     dfBandit: n x 2 dataframe, with colnames "Arm", "Result" required. If day_split = TRUE, then needs to be n x 3 with colnames "Arm", "Result", "day"
#     n_sim: number of simulations
#     n_per_sim: number of observations per sim
#     policy_type: string, one of "greedy", "UCB", "TS
#     dfParams: dataframe with parameters, suitable for the policy. 
#           greedy: exploration, eps
#           UCB: exploration, fC
#           TS: exploration
###

sim_params <- function(dfBandit, n_sim, n_per_sim, policy_type, dfParams){
  
  # if policy is greedy
  if(policy_type == "greedy"){
    results_params <- apply(dfParams, MARGIN=1,function(row_param){
      print(paste0("Running simulations for the ", policy_type, " policy, with parameters exploration: ", row_param[1], " and epsilon: ", row_param[2]))
      
      # run the sim with exploration and epsilon param
      exploration_sim <- row_param[1]
      eps_sim <- row_param[2]
      result_for_param <- sim_experiment(dfBandit, n_sim, n_per_sim, policy_type, exploration = exploration_sim, eps = eps_sim)
      
      # get the aggregate result
      aggregate_result <- calc_rewards_sims(result_for_param)
      aggregate_result$param <- paste0("exploration=", exploration_sim, ", epsilon=", eps_sim)
      
      # save parameters, all individual results, and return
      params <- row_param
      results <- list(aggregate = aggregate_result, details = result_for_param, params = params)
      
      return(results)
    })
  }else if (policy_type == "UCB"){
    results_params <- apply(dfParams, MARGIN=1,function(row_param){
      print(paste0("Running simulations for the ", policy_type, " policy, with parameters exploration: ", row_param[1], " and C: ", row_param[2]))
      
      # run the sim with exploration and C param
      exploration_sim <- row_param[1]
      C_sim <- row_param[2]
      result_for_param <- sim_experiment(dfBandit, n_sim, n_per_sim, policy_type, exploration = exploration_sim, fC = C_sim)
      
      # get the aggregate result
      aggregate_result <- calc_rewards_sims(result_for_param)
      aggregate_result$param <- paste0("exploration=", exploration_sim, ", C=", C_sim)
      
      # save parameters, all individual results, and return
      params <- row_param
      results <- list(aggregate = aggregate_result, details = result_for_param, params = params)
      
      return(results)
    })
  }else if (policy_type == "TS"){
    results_params <- apply(dfParams, MARGIN=1,function(row_param){
      print(paste0("Running simulations for the ", policy_type, " policy, with parameters exploration: ", row_param[1]))
      
      # run the sim with exploration and C param
      exploration_sim <- row_param[1]
      result_for_param <- sim_experiment(dfBandit, n_sim, n_per_sim, policy_type, exploration = exploration_sim)
      
      # get the aggregate result
      aggregate_result <- calc_rewards_sims(result_for_param)
      aggregate_result$param <- paste0("exploration=", exploration_sim)
      
      # # save parameters, all individual results, and return
      params <- row_param
      results <- list(aggregate = aggregate_result, details = result_for_param, params = params)
      
      return(results)
    })
  }
  return(results_params)
}


