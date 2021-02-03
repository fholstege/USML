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
               reshape2) 

# start of each file name with yahoo data
start_fileName <- "Data/YahooOpenData/ydata-fp-td-clicks-v1_0.2009050"

# for now, just start with the first day
dfYahoo <- read.table(paste0(start_fileName, "1.gz"), nrow=10000,header = FALSE, fill = TRUE, sep = " ")[,1:3]

# only grab arm and result data, call columns
dfBandit_sim <- dfYahoo[,2:3]
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
sim_policy <- function(dfBandit, policy_type, exploration,...){
  
  ## Part 1: set parameters 
  n_obs <- length(dfBandit$result)
  n_explore <- round(n_obs* exploration,0)
  n_exploit <- n_obs - n_explore
  
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
    }
    
    # get dataframe with remaining results for chosen arm, and pick a random instance
    chosen_arm_df <- ldf_arms[[chosen_arm]]

    # obtain the result, remove this from the overall remaining data
    index_result <- sample(1:length(chosen_arm_df$result), 1)
    result <- chosen_arm_df[index_result,]$result

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
  
  # save data in list
  results <- list(Overall = dfResult, 
                  after_explore = dfResult_exploit, 
                  before_explore = dfResult_explore, 
                  total = rbind(dfResult_explore,dfResult_exploit),
                  total_succes_rate = sum(dfResult$succes_size)/sum(dfResult$sample_size))

  return(results)
  
}


sim_experiment <- function(dfBandit, n_sim, policy_type, exploration,...){
  
  
  
  
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




 


#check greedy with epsilon of 0.1, 0.05, 0.001 and exploration of 20% of the data
result_sim_eps_01 <- sim_policy(dfBandit_sim, "greedy", 0.2, eps=0.1)
result_sim_eps_005 <- sim_policy(dfBandit_sim, "greedy", 0.2, eps=0.05)
result_sim_eps_001 <- sim_policy(dfBandit_sim, "greedy", 0.2, eps=0.01)

#check greedy with epsilon of 0.1, 0.05, 0.001 and exploration of 20% of the data
result_sim_c_01 <- sim_policy(dfBandit_sim, "UCB", 0.2, fC=0.1)
result_sim_c_02 <- sim_policy(dfBandit_sim, "UCB", 0.2, fC=0.2)
result_sim_c_03 <- sim_policy(dfBandit_sim, "UCB", 0.2, fC=0.3)

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
