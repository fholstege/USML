# data directory
data_dir = "C:/Users/flori/OneDrive/Documents/Tinbergen/Courses/USML/Data/YahooOpenData/"



read_YahooData <- function(data_dir, rsample.size = 1000000){
  
  # set temporary directory for data
  setwd(data_dir)
  
  # Obtain files from folder
  files <- grep(".*gz", list.files(), value=T)
  
  print(paste0("Reading and cleaning the following files: ", files))
  
  # Initialize final data matrix
  df <- as.data.frame(matrix(0, nrow=length(files)*rsample.size, ncol=9))
  
  # Loop over files in folder
  result_loading <- lapply(files, function(file){
    
    print(file)
    
    # load data file
    rawdata <- read.table(file, header = F, sep = " ",fill=TRUE, blank.lines.skip = FALSE, colClasses = c("NULL", "character", "character",rep("character", 7),rep("NULL", 140))) 
    
    colnames(rawdata)[1:2] <- c("Arm", "Reward")
    
    # filter out wrong observations
    data <- rawdata[which(rawdata$Reward == 0 | rawdata$Reward == 1),]

    if(nrow(data) < rsample.size){
      print("Requested sample size larger than available data - taking all data")
      rsample.size <- nrow(data)
    }
    
    # take a random sample of size rsample.size
    r.ind <- sample(nrow(data), size = rsample.size)
    data <- data[r.ind,]
    
    print(data)
    
    return(data)
    
   
    
  })
  
  return(result_loading)
}


result_loading <- read_YahooData(data_dir, rsample.size = 2000000)


library(rlist)

i = 1
setwd("FinalProject/Data/")


write.csv(data.frame(result_loading[1]), "userData_day1.csv")
write.csv(data.frame(result_loading[2]), "userData_day2.csv")
write.csv(data.frame(result_loading[3]), "userData_day3.csv")
write.csv(data.frame(result_loading[4]), "userData_day4.csv")
write.csv(data.frame(result_loading[5]), "userData_day5.csv")
write.csv(data.frame(result_loading[6]), "userData_day6.csv")
write.csv(data.frame(result_loading[7]), "userData_day7.csv")
write.csv(data.frame(result_loading[8]), "userData_day8.csv")
write.csv(data.frame(result_loading[9]), "userData_day9.csv")
write.csv(data.frame(result_loading[10]), "userData_day10.csv")


# add repository to make mclapply() run in parallel (only necessary on windows)
install_github('nathanvan/parallelsugar')
library(parallelsugar)

# get all the raw data 
dfRawData <- read.csv("Data/dfUserData_raw.csv")

# get the features on users
dfUserFeatures <- dfRawData[,6:10]

# clean the features on users
cleaned_user_features <- apply(dfUserFeatures, 2, function(x){
  
  cleaned_col <- as.numeric(gsub(".*:","",x))
  return(cleaned_col)
})

# create dataframe with clean user features
dfClean_user_features <- data.frame(unlist(cleaned_user_features))

colnames(dfRawData)[3:4]

# combine clean user features with arm and reward data
dfComplete <- cbind(dfRawData[,3:4],dfClean_user_features)
colnames(dfComplete) <- c("arm", "reward","user1", "user2", "user3", "user4", "user5")

write.csv(dfComplete,"Data/dfUserData_clean.csv")



# load the clean data
dfUserData <- read.csv("Data/dfUserData_clean.csv")
n_per_day <- 2000000

# create day column to be able to split per day
day <- c()
for(i in 1:10){
  
  day <- c(day, rep(i,n_per_day))
}
colnames(dfUserData)[1] <- "day"
dfUserData$day <- day

# ensure arms start with 1 to ensure it works in simulation
i <- 1
for(a in unique(dfComplete$arm)){
  dfComplete[dfComplete$arm == a, "arm"] <- i
  i <- i + 1
}

# set arm and reward to numerics
dfUserData$reward <- as.numeric(dfComplete$reward)
dfUserData$arm <- as.numeric(dfComplete$arm)


# Agent$new(LinUCBDisjointPolicy$new(alpha = 0.5), bandit, "Lin UCB - Alpha=0.5"),
# Agent$new(LinUCBDisjointPolicy$new(alpha = 1), bandit, "Lin UCB - Alpha=1"),
# Agent$new(ContextualLinTSPolicy$new(v = 0.2), bandit, "Lin TS - v=0.2"),
# Agent$new(ContextualLinTSPolicy$new(v = 0.1), bandit, "Lin TS - v=0.1"), 
# Agent$new(UCB1Policy$new(), bandit, "UCB1 - Alpha=0.25"),
# Agent$new(ThompsonSamplingPolicy$new(), bandit, "Thompson"))




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

