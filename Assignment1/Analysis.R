
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

### Compare our sparse PCA function to the PMD() function in R
# our sparse PCA
result_ours <- sparce_PCA(as.matrix(cityweather), c1=1.5,c2=1.5)
result_ours$v
result_ours$u
result_ours$sigma

# sparce PCA from PMA package
result_pmd <- PMD(as.matrix(cityweather), sumabsu = 1.5, sumabsv = 1.5, center = FALSE)
result_pmd$v
result_pmd$u
result_pmd$d

# the results appear very, very similar - when sumabsu, sumabsv = 1, then difference of 1-e16. No difference when other values are picked. 





