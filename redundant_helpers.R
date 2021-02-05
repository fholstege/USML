
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