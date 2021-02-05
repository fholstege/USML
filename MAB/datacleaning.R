
data.cleaning <- function(funcdir, datadir, nrow = NULL, rsample.size = NULL) {
  '
  Inputs:
    funcdir         directory where functions are saved
    datadir         directory where data can be found
    nrow            import fixed number of rows from each file (max on my PC: n=250,000)
    rsample.size    random sample of data per file
    
  Purpose:
    Import data from multiple files and store as csv
  
  Returns:
    data        dataframe with combined data from files
  '

  # set temporary directory for data
  setwd(datadir)
  
  # Obtain files from folder
  files <- grep(".*gz", list.files(), value=T)
  
  # Initialize final data matrix
  df <- as.data.frame(matrix(0, nrow=length(files)*rsample.size, ncol=2))

  # Loop over files in folder
  i <- 0
    for (file in files) {
      i <- i + 1
      print(i)
      
      # load data file
      rawdata <- read.table(file, header = F, sep = " ", colClasses = c("NULL", "character", "character", rep("NULL", 147)),
                            fill=TRUE, blank.lines.skip = FALSE) 
      
      # set new column names
      colnames(rawdata) <- c("Arm", "Reward")
      
      # filter out wrong observations
      data <- rawdata[which(rawdata$Reward == 0 | rawdata$Reward == 1),]
      
      # take a random sample of size rsample.size
      r.ind <- sample(nrow(data), size = rsample.size)
      data <- data[r.ind,]
      
      # store data in pre-created matrix
      df[((i-1)*rsample.size + 1):(i*rsample.size), ] <- as.matrix(data)
    }
  
  # create dataframe
  df <- as.data.frame(df)
  
  # go back to initial directory 
  setwd(funcdir)
  
  # store data in gz file
  write.csv(df, paste("all_data.csv", sep=''))
  
  return(df)
  
}
