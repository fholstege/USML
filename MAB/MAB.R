#################
# Authors: Floris Holstege, Yuchou Peng, Chao Liang, Eva Mynott, Gabriele Mingoli
# Purpose: Code to clean the data file for yahoo data - MAB assignment 
# Sections: 
#           A) Load packages
#           B) Define policy functions
#           C) Define functions to run simulations 
#           D) Define functions to analyse the results (plots) 
#           E) Run simulations to get results per method 
#           F) Compare to package
#################


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
               utils,
               readr,
               purrr,
               dplyr,
               furrr,
               Rfast,
               devtools,
               parallel,
               TTR) 
# make sure runs in parallel
future::plan(multiprocess)

# add repository to make mclapply() run in parallel (only necessary on windows)
install_github('nathanvan/parallelsugar')
library(parallelsugar)


# start of each file name with yahoo data
dfArticles <- read.csv("Data/all_data_n=10000000.csv")[,2:3]
colnames(dfArticles)<- c("Arm", "result")
n_per_day <- 1000000

# create day column
day <- c()

for(i in 1:10){
  
  day <- c(day, rep(i,n_per_day))
}

dfArticles$day <- day

# create small df to test
df_small_sim <- dfArticles[1:1000,]
dfSims <- dfArticles[sample(1:1000000, 1000000),]

# set seed to ensure reproducability
set.seed(123)

################################################################################
# B) define policies
################################################################################


###
# policy_greedy: picks an arm, basedon the greedy algorithm for multi-armed bandits
# 
#   Arguments: 
#     eps: percentage of times a random arm needs to be picked
#     index_best_arm; index of the current, best arm - picked in 1-eps\% of cases
#
#   Output:
#     chosen_arm; index of the arm chosen
###
policy_greedy <- function(index_arms, index_best_arm, eps){
  
  if (runif(1,min=0,max=1) < eps){
    chosen_arm <- sample(index_arms,1)
  }else{
    chosen_arm <- index_best_arm
  }
  
  return(chosen_arm)
}

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
  
  Est_val_arms <- dfResults_arms$succes_rate
  log_t <- log(sum(dfResults_arms$sample_size))
  
  N_t_arms <- dfResults_arms$sample_size

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
  
  #print("Sampled")
  #print(lObs_result)

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
  
  ## Part 1: set parameters 
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
    
    # remove said result in the data, and allocate new dataframe without that result back to list of dataframes
    #ldf_arms[[chosen_arm]]<- chosen_arm_df[-index_result,]

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

sim_experiment <- function(dfBandit, n_sim, n_per_sim, policy_type, exploration,day_split=TRUE,...){

  if(day_split){
    
    ldf_perDay <- split(dfBandit[,1:2] , f = dfBandit$day )
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

sim_params <- function(dfBandit, n_sim, n_per_sim, policy_type, dfParams){


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
      
      # save details and parameters
      details_result <- result_for_param
      params <- row_param
      
      print(nrow(aggregate_result))
      
      
      results <- list(aggregate = aggregate_result, details = details_result, params = params)
      
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
      
      # save details and parameters
      details_result <- result_for_param
      params <- row_param
      
      results <- list(aggregate = aggregate_result, details = details_result, params = params)
      
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
      
      # save details and parameters
      details_result <- result_for_param
      params <- row_param
      

      results <- list(aggregate = aggregate_result, details = details_result, params = params)
      
      return(results)
    
    })
  }
  
  return(results_params)
  

}



  
  

################################################################################
# D) Define functions to analyse the results  
################################################################################

create_arm_dummy <- function(arm_name, df){
  return(ifelse(df$Arm == arm_name, 1,0))
}

calc_chosenArm_experiment <- function(df){
  
  # get per time step the chosen arm 
  arms <- unique(df$Arm)
  lResult_dummies <- lapply(arms, create_arm_dummy, df=df)
  count_chosen <- list.cbind(lResult_dummies)
  
  # calculate the % chosen per arm over time
  cumsum_chosen <- apply(count_chosen, MARGIN=2, cumsum)
  index_arms <- 1:length(df$result)
  dfPerc_chosen <- data.frame(index_arms,cumsum_chosen/index_arms)
  
  # change colnames
  colnames(dfPerc_chosen) <- c("index", arms)
  
  
  return(dfPerc_chosen)
  
}

create_armChoice_plot <- function(dfPerc_chosen){

  # add column names, melt dataframe for ggplot viz
  dfPerc_chosen_melted <- melt(dfPerc_chosen, id.vars = "index")

  # define ggplot visualization
  armChoice_plot <- ggplot(data=dfPerc_chosen_melted %>% filter(index >10), aes(x=index, y = value, fill = variable))+
                    geom_area()+
                    theme_bw()+
                    labs(x = "Time Steps", y = "% Of Arm chosen", col = "Article ID")

  
  
  return(armChoice_plot)
}


calc_rewards_sims <- function(result_sim){
  

  # get list with all the results
  lReward_sims <- lapply(result_sim, function(x){return(x$total$result)})

  # get avg. total reward
  dfTot_reward_sims <- list.cbind(lReward_sims)
  dfTot_cumReward <- apply(dfTot_reward_sims, 2, cumsum)
  dfAvg_cumReward <- apply(dfTot_cumReward, 1, mean)
  
  n <- ncol(dfTot_cumReward)
  dfSD <- apply(dfTot_cumReward, 1, sd)

  lower_bound <- pmax(dfAvg_cumReward - (1.96 * dfSD/sqrt(n)),0)
  upper_bound <- dfAvg_cumReward + (1.96 * dfSD/sqrt(n))
  
  
  # return these in dataframe
  dfSim_plot <- data.frame(index = 1:length(dfAvg_cumReward), avg_cumReward = dfAvg_cumReward, lower = lower_bound, upper = upper_bound)
  
  return(dfSim_plot)
  
  
}

################################################################################
# E) Run simulations to get results
################################################################################

# first, define parameters for simulation
n_sims = 100
n_per_sim = 10000

# 27 combinations
vEps <- c(0.2,0.1,0.05)
vC <- c(0.3,0.2,0.1)
vExploration <- c(0.1,0.2,0.3)

nrow(dfSims)

# create df for sim_params functions - greedy and TS
dfParams_eps <- expand.grid(vExploration, vEps)
colnames(dfParams_eps) <- c("exploration", "eps")

# create df for sim_params functions - UCB
dfParam_UCB <-  expand.grid(vExploration, vC)
colnames(dfParam_UCB) <- c("exploration", "C")

# create df for sim_params functions - TS
dfParam_TS <- expand.grid(vExploration)
colnames(dfParam_TS) <- c("exploration")

param_results_greedy <- sim_params(dfSims, n_sims, n_per_sim, "greedy", dfParams_eps)
param_results_C <- sim_params(dfSims, n_sims, n_per_sim, "UCB", dfParam_UCB)
param_results_TS <- sim_params(dfSims, n_sims, n_per_sim, "TS", dfParam_TS)


dfAggregate_results_eps <- rbind(param_results_greedy[[1]]$aggregate, param_results_greedy[[4]]$aggregate, param_results_greedy[[7]]$aggregate)
dfAggregate_results_eps_explore <- rbind(param_results_greedy[[1]]$aggregate, param_results_greedy[[2]]$aggregate, param_results_greedy[[3]]$aggregate)
dfAggregate_results_C <- rbind(param_results_C[[1]]$aggregate, param_results_C[[4]]$aggregate, param_results_C[[7]]$aggregate)
dfAggregate_results_C_explore <- rbind(param_results_C[[1]]$aggregate, param_results_C[[2]]$aggregate, param_results_C[[3]]$aggregate)
dfAggregate_results_TS <- rbind(param_results_TS[[1]]$aggregate, param_results_TS[[2]]$aggregate, param_results_TS[[3]]$aggregate)

ggplot(data = dfAggregate_results_eps, aes(x =index, y = avg_cumReward, col=param)) +
  geom_line()+
  geom_ribbon(aes(ymax = upper, ymin = lower, fill = param), alpha=0.1, colour = NA)+
  theme_bw()+
  labs(x = "Time Steps", y = "Avg. Total Clicks", title="Greedy")

ggplot(data = dfAggregate_results_C, aes(x =index, y = avg_cumReward, col=param)) +
  geom_line()+
  geom_ribbon(aes(ymax = upper, ymin = lower, fill = param), alpha=0.1, colour = NA)+
  theme_bw()+
  labs(x = "Time Steps", y = "Avg. Total Clicks", title="UCB,  exploration = 10%")



ggplot(data = dfAggregate_results_TS, aes(x =index, y = avg_cumReward, col=param)) +
  geom_line()+
  geom_ribbon(aes(ymax = upper, ymin = lower, fill = param), alpha=0.1, colour = NA)+
  theme_bw()+
  labs(x = "Time Steps", y = "Avg. Total Clicks", title="TS,  exploration = 10%")

################################################################################
# F) Compare to package
################################################################################

ts <- sim_policy(df_small_sim[,1:2], "TS", exploration=0.3, size_experiment=1000, only_total = FALSE)

ts_perc_chosen <- calc_chosenArm_experiment(ts$total)

create_armChoice_plot(ts_perc_chosen)

# simulate with package - key difference is that they seem to know beforehand what the top choice is
horizon <- 5000
simulations <- 1 
conversionProbabilities <- dfSummary$succes_rate
bandit <- BasicBernoulliBandit$new(weights = conversionProbabilities)
policy <- EpsilonGreedyPolicy$new(epsilon = 0.10) 
agent <- Agent$new(policy, bandit) 
historyEG <- Simulator$new(agent, horizon, simulations)$run() 
plot(historyEG, type = "arms",      legend_title = 'Epsilon Greedy',      legend_position = "topright",      smooth = TRUE) 



