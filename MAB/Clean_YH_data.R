#################
# Authors: Floris Holstege, Yuchou Peng, Chao Liang, Eva Mynott, Gabriele Mingoli
# Purpose: Code to clean the data file for yahoo data - MAB assignment 
# Sections: 
#           A) Load packages
#           B) Load Data
#           C) Clean data
#           D) Define policies
#           E) Define agents 
#           F) Define simulation function
#           G) Check results vs "contextual" packages
#################


################################################################################
# A) Load packages
################################################################################

# empty global environment
rm(list=ls())

# Packes required for subsequent analysis. P_load ensures these will be installed and loaded. 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(contextual) 

help(read.table)
n_files <- 1:10
start_fileName <- "Data/YahooOpenData/ydata-fp-td-clicks-v1_0.2009050"


dfYahoo <- read.table(gzfile("Data/YahooOpenData/ydata-fp-td-clicks-v1_0.20090501.gz"), nrow=100000,header = FALSE, fill = TRUE, sep = " ")
