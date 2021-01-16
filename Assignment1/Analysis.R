
# Packes required for subsequent analysis. P_load ensures these will be installed and loaded. 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(PMA, 
               ggplot2)

# load the data
load("../Data/cityweather.RData")

# load helper functions
source("../Helpers/utilities.R")
source("PCA.R")

# set seed to ensure reproducability
set.seed(123)


### Compare our PCA function to the prcomp() function in R
# We get the exact same results as prcomp 
result_prcomp <- prcomp(scale(as.matrix(cityweather)))
result_ours <- svd_components(scale(as.matrix(cityweather)))

### create an scree plot
screePlot <- create_screePlot(as.matrix(scale(cityweather)))
screePlot


### Compare our sparse PCA function to the PMD() function in R
result_ours <- sparce_PCA(as.matrix(cityweather), c1=1,c2=1)
result_ours$v
result_ours$u
result_ours$sigma

result_pmd <- PMD(as.matrix(cityweather), sumabsu = 1, sumabsv = 1, center = FALSE)
result_pmd$v
result_pmd$u
result_pmd$d

# the results appear quite similar, but not exactly - not sure if due to convergence, our initial u


result_pmc <- SPC(as.matrix(cityweather),sumabs = 1, center = FALSE, K=1)
result_pmc$v
result_pmc$u
result_pmc$d

SPC
