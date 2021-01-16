
# Packes required for subsequent analysis. P_load ensures these will be installed and loaded. 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(PMA, 
               ggplot2,
               factoextra,
               xtable)

# load the data
load("../Data/FIFA2017_NL.RData")

# load helper functions
source("../Helpers/utilities.R")
source("PCA.R")

# set seed to ensure reproducability
set.seed(123)

# create dataframe with the player attributes
df_fifa_attributes <- fifa[,1:32]

# put attributes (without player name, club etc. ) in matrix, and scale these
mX_attributes <- scale(as.matrix(df_fifa_attributes[,4:ncol(df_fifa_attributes)]))

### Compare our PCA function to the prcomp() function in R
# We get the exact same results as prcomp 
result_prcomp <- prcomp(mX_attributes)
result_PCA <- svd_components(mX_attributes)

### create an scree plot
screePlot <- create_screePlot(result_PCA)

## show how variance improves with each component added
varPlot <- create_varPlot(result_PCA)

## create table of the loadings 
create_LoadingTable(result_PCA, K=5)



### Compare our sparse PCA function to the PMD() function in R, for rank-1
# the results appear very, very similar - when sumabsu, sumabsv = 1, then difference of 1-e16. No difference when other values are picked.
# our sparse PCA
result_sparse_svd <- sparse_svd(mX_attributes, c1=1.5, c2=1.5)
result_sparse_svd$v
result_sparse_svd$u
result_sparse_svd$sigma

# sparce PCA from PMA package
result_pmd <- PMD(mX_attributes, sumabsu = 1.5, sumabsv = 1.5, center = FALSE)
result_pmd$v
result_pmd$v.init
result_pmd$u
result_pmd$d

### Compare our sparse PCA function to the PMD() function with a higher rank 
result_pmd_rank <- PMD(mX_attributes, sumabsu = 1.5, sumabsv = 1.5, K=5)
result_pmd_rank$v
result_pmd_rank$u

result_sparse_PCA <- sparse_PCA(mX_attributes, c1= 1.5,c2=1.5,K=5)
result_sparse_PCA$v
result_sparse_PCA$u

abs(result_sparse_PCA$v) - abs(result_pmd_rank$v)

