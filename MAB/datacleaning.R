
data.cleaning <- function(funcdir, datadir, nrow = NULL, rsample = NULL) {
  '
  Inputs:
    funcdir     directory where functions are saved
    datadir     directory where data can be found
    nrow        import fixed number of rows from each file (max on my PC: n=250,000)
    rsample     random sampling does not work yet
    
  Purpose:
    Import data from multiple files and store as csv
  
  Returns:
    data        dataframe with combined data from files
  '

  # set temporary directory for data
  setwd(datadir)
  
  # Obtain files from folder
  files <- grep(".*gz", list.files(), value=T)
  
  # Determine file sample size
  if (is.null(nrow) & is.null(rsample)) {
    n <- 10000
  } else if (is.null(rsample)) {
    n <- nrow
  } else {
    n <- rsample
  }
  
  # Initialize final data matrix
  data <- as.data.frame(matrix(0, nrow=length(files)*n, ncol=3))
  
  # Loop over files in folder
  tryCatch(
    expr = 
      {i <- 0
      for (file in files) {
        i <- i + 1
        
        # load data file
        rawdata <- read.table(file, header = F, sep = " ", fill = T, nrows = n)[,1:3] 
    
        # store data in pre-created matrix
        data[((i-1)*n + 1):(i*n), ] <- as.matrix(rawdata)
      }
    },
    # if an error occurs
    error = function(err) {
        setwd(funcdir)
        print('Data importing went wrong. Try decreasing your sample size.')
      }
    )
  
  # create dataframe
  data <- as.data.frame(data)
  
  # go back to initial directory 
  setwd(funcdir)
  
  # store data in gz file
  write.csv(data, paste("data_n=",n,".csv", sep=''))
  
  return(data)
  
}
