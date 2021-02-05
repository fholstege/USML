#################
# Authors: Floris Holstege, Yuchou Peng, Chao Liang, Eva Mynott, Gabriele Mingoli
# Purpose: Code to clean the data file for yahoo data - MAB assignment 
# Sections: 
#           A) Load packages
#           B) Define policy functions
#           C) Define functions to run simulations 
#           D) Define functions to analyse the results (plots) 
#           E) Compare to package
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
               Rfast) 
# make sure runs in parallel
future::plan(multiprocess)


# start of each file name with yahoo data
start_fileName <- "Data/YahooOpenData/ydata-fp-td-clicks-v1_0.200905"

#specify how many to grab
n_grab = 1000
n_grab_in_batch <- 100


# function for cleaning data
clean_yahooData <- function(df){
  
  df_necessary <- df[,2:3]
  colnames(df_necessary) <- c("Arm", "result")
  
  check_if_error <- ifelse(df_necessary$result== 1 | df_necessary$result == 0, TRUE, FALSE)
  
  df_cleaned <- df_necessary[check_if_error,]
  
  return(df_cleaned)
}

# for all the files
list_file_ends <- c("01.gz", "02.gz", "03.gz", "04.gz", "05.gz", "06.gz", "07.gz", "08.gz", "09.gz", "10.gz")
list_file_names <- paste0(rep(start_fileName, 10),list_file_ends)

#function to grab random sample
get_random_sample_fromFile <- function(fileName, n_grab, n_grab_in_batch){
  
  # get how many per batch 
  n_batch = n_grab/n_grab_in_batch
  
  # where to start at (random), ordered
  start_at  <- floor(runif(n_batch, min = 1, max = (n_grab - n_grab_in_batch) ))
  start_at  <- start_at[order(start_at)]
  
  dfRandom <- future_map_dfr(start_at, ~clean_yahooData(read.table(fileName, nrow= n_grab_in_batch, fill = TRUE, sep = " ",skip = .x))) 
  
  return(dfRandom)
  
}

list_random_samples <- lapply(list_file_names, get_random_sample_fromFile, n_grab = n_grab, n_grab_in_batch = n_grab_in_batch)
dfBandit_sim <- list.rbind(list_random_samples)
colnames(dfBandit_sim) <- c("Arm", "result")

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
policy_greedy <- function(eps, index_arms, index_best_arm){
  
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


policy_TS <- function(dfResult, n_arms){
  
  n <- rep(100, n_arms)
  
  alpha <- dfResult$alpha
  beta <- dfResult$beta
  
  lObs_result <- lapply(n, rbeta, alpha, beta)
  mean_distribution <- unlist(lapply(lResult, mean))
  
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
sim_policy <- function(dfBandit, policy_type, exploration,size_experiment=10000,only_total=TRUE,...){
  
  
  
  ## Part 1: set parameters 
  n_obs <- size_experiment
  n_explore <- round(n_obs* exploration,0)
  n_exploit <- n_obs - n_explore
  
  # get random of size
  dfBandit <- dfBandit[sample(1:n_obs, size_experiment), ]

  ## part 2: exploration phase 
  # Set the results for the exploration phase by finding a certain number of random observations
  index_exploration <- sample(1:n_obs,n_explore)
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
  Arms <- unique(dfResult$Arm)
  index_arms <- 1:length(Arms)
  
  # get list of dataframes (one per arm & its results) for remainder
  ldf_arms <- split(dfBandit_postExplore , f = dfBandit_postExplore$Arm )
  
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
      chosen_arm <- policy_TS(dfResult, n_arms)
    }
    
    # get dataframe with remaining results for chosen arm, and pick a random instance
    chosen_arm_df <- ldf_arms[[chosen_arm]]
    

    # obtain the result, remove this from the overall remaining data
    index_result <- sample(1:length(chosen_arm_df$result), 1)
    result <- chosen_arm_df[index_result,]$result
    
    if(length(result) == 0){
      print(i)
      cat("You have run out of observations from a chosen arm - decrease the size of the experiment, or increase the size of the dataset")
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
    
    # update the best arm
    index_best_arm <- which.max(dfResult$succes_rate)
    
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


sim_experiment <- function(dfBandit, n_sim, policy_type, exploration,...){
  

  # get n of observations, and get how many you want in sim
  n_obs <- length(dfBandit$result)
  n_in_sim <- n_obs/n_sim
  
  # get list of the dataframes - same dataframe, n_sim times
  list_dfBandit_exp <- rep(list(dfBandit),n_sim)
    
  # get result 10x
  result_per_random_sample <- lapply(list_dfBandit_exp, sim_policy, policy_type = policy_type, exploration = exploration, size_experiment=n_in_sim,...)
  
  return(result_per_random_sample)
}
  
  
  
  

################################################################################
# D) Define functions to analyse the results (plots) 
################################################################################

create_arm_dummy <- function(arm_name, df){
  return(ifelse(df$Arm == arm_name, 1,0))
}



create_armChoice_plot <- function(dfResult){
  
  # get per time step the chosen arm 
  arms <- unique(dfResult$Arm)
  lResult_dummies <- lapply(arms, create_arm_dummy, df=dfResult)
  count_chosen <- list.cbind(lResult_dummies)
  
  # calculate the % chosen per arm over time
  cumsum_chosen <- apply(count_chosen, MARGIN=2, cumsum)
  index_arms <- 1:length(dfResult$result)
  dfPerc_chosen <- data.frame(index_arms,cumsum_chosen/index_arms)
  
  # add column names, melt dataframe for ggplot viz
  colnames(dfPerc_chosen) <- c("index", arms)
  dfPerc_chosen_melted <- melt(dfPerc_chosen, id.vars = "index")

  # define ggplot visualization
  armChoice_plot <- ggplot(data=dfPerc_chosen_melted %>% filter(index >10), aes(x=index, y = value, col = variable))+
                    geom_line()+
                    theme_bw()+
                    labs(x = "Time Steps", y = "% Of Arm chosen", col = "Article ID")

  
  
  return(armChoice_plot)
}



n_sims = 100
result_sims <- sim_experiment(dfBandit_sim, n_sims,"greedy", 0.3, eps=0.1)

lReward_sims <- lapply(result_sims, function(x){return(x$total$result)})
dfTot_reward_sims <- list.cbind(lReward_sims)
dfTot_cumReward <- apply(dfTot_reward_sims, 2, cumsum)
dfAvg_cumReward <- apply(dfTot_cumReward, 1, mean)

ind_5perc <- n_sims/20
ind_95perc <- n_sims - (n_sims/20)


lower_bound <- apply(dfTot_cumReward, 1, nth, k=ind_5perc)
upper_bound <- apply(dfTot_cumReward, 1, nth, k=ind_95perc)


dfSim_plot <- data.frame(index = 1:1000, avg_cumReward = dfAvg_cumReward, lower = lower_bound, upper = upper_bound)


ggplot(data = dfSim_plot, aes(x =index , y = avg_cumReward)) +
  geom_line()+
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey70", alpha = 0.2) +
  theme_bw()+
  labs(x = "Time Steps", y = "Avg. Total Clicks")



#check greedy with epsilon of 0.1, 0.05, 0.001 and exploration of 20% of the data
result_sim_eps_01 <- sim_policy(dfBandit_sim, "greedy", 0.2, eps=0.1)
result_sim_eps_005 <- sim_policy(dfBandit_sim, "greedy", 0.2, eps=0.05)
result_sim_eps_001 <- sim_policy(dfBandit_sim, "greedy", 0.2, eps=0.01)

#check greedy with epsilon of 0.1, 0.05, 0.001 and exploration of 20% of the data
result_sim_c_01 <- sim_policy(dfBandit_sim, "UCB", 0.2, fC=0.1)
result_sim_c_02 <- sim_policy(dfBandit_sim, "UCB", 0.2, fC=0.2)
result_sim_c_03 <- sim_policy(dfBandit_sim, "UCB", 0.2, fC=0.3)

result_TS<- sim_policy(dfBandit_sim, "TS", 0.2)


length(unique(dfYahoo_day1$Arm))
unique(dfYahoo_day1$Arm)

# visualize the results for epsilon difference 
dfResult_sims_eps <- data.frame(1:length(result_sim_eps_01$total$result), 
                                  cumsum(result_sim_eps_01$total$result), 
                                  cumsum(result_sim_eps_005$total$result),
                                  cumsum(result_sim_eps_001$total$result))
colnames(dfResult_sims_eps) <- c("index", "total_reward_eps_01", "total_reward_eps_005", "total_reward_eps_001")
dfResult_eps_melted <- melt(dfResult_sims_eps, id.vars = "index")

ggplot(data = dfResult_eps_melted, aes(x = index, y = value, col = variable)) +
  geom_line()+
  theme_bw()+
  labs(x = "Time Steps", y = "Total Clicks", col = "Parameter")+
  scale_color_discrete(name = "Epsilon = ", labels = c("10%", "5%", "1%"))



# visualize the results for c difference 
dfResult_sims_c <- data.frame(1:length(result_sim_c_01$total$result), 
                                cumsum(result_sim_c_01$total$result), 
                                cumsum(result_sim_c_02$total$result),
                                cumsum(result_sim_c_03$total$result))
colnames(dfResult_sims_c) <- c("index", "total_reward_c_01", "total_reward_c_02", "total_reward_c_03")
dfResult_c_melted <- melt(dfResult_sims_c, id.vars = "index")

ggplot(data = dfResult_c_melted, aes(x = index, y = value, col = variable)) +
  geom_line()+
  theme_bw()+
  labs(x = "Time Steps", y = "Total Clicks", col = "Parameter")+
  scale_color_discrete(name = "C = ", labels = c("10%", "20%", "30%"))


# check for thompson sampling 
dfResult_sims_ts <- data.frame(1:length(result_TS$total$result),
                               cumsum(result_TS$total$result))

colnames(dfResult_sims_ts) <- c("index", "total_reward")
dfResult_ts_melted <- melt(dfResult_sims_ts, id.vars="index")
colnames(dfResult_ts_melted)
ggplot(data = dfResult_ts_melted, aes(x = index, y = value, col = variable)) +
  geom_line()+
  theme_bw()+
  labs(x = "Time Steps", y = "Total Clicks", col = "Parameter")

create_armChoice_plot(result_TS$total)


# get the actual probabilities 
dfSummary <- dfBandit_sim %>% 
  group_by(Arm)%>%
  summarise(succes_size = sum(result),
            sample_size = n(),
            succes_rate = succes_size/sample_size)
dfSummary$succes_rate



# simulate with package - key difference is that they seem to know beforehand what the top choice is
horizon <- 5000
simulations <- 1 
conversionProbabilities <- dfSummary$succes_rate
bandit <- BasicBernoulliBandit$new(weights = conversionProbabilities)
policy <- EpsilonGreedyPolicy$new(epsilon = 0.10) 
agent <- Agent$new(policy, bandit) 
historyEG <- Simulator$new(agent, horizon, simulations)$run() 
plot(historyEG, type = "arms",      legend_title = 'Epsilon Greedy',      legend_position = "topright",      smooth = TRUE) 


n_arms <- 5
n <- rep(5, n_arms)

alpha <- c(1,2,3, 4, 5)
beta <- c(5,4,3,2,1)

lObs <- lapply(n, rbeta, alpha, beta)
mean_distribution <- unlist(lapply(lResult, mean))
  
vectorize_rbeta <- function()
mean(rbeta(100, 7,97))


