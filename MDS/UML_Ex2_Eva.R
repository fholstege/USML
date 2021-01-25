################################################################################
# COURSE: Unsupervised Machine Learning
#
# Exercise 2
#
# Team 4
#
# DATE: 2021-01-23
################################################################################

rm(list=ls())

################################################################################
# Initialize local settings
################################################################################

# Specify working directory
setwd('C:/Users/evamy/Documents/1. University/MPhil Tinbergen/Year 2/Block 3/UML/Exercise 2/')

# load data
load("basket.RData")

# load dependencies
if(!require('smacof'))install.packages('smacof', quiet=T)
if(!require('ggplot2'))install.packages('ggplot2', quiet=T)
if(!require('rlist'))install.packages('rlist', quiet=T)

set.seed(42)

################################################################################
# Data prep
################################################################################

# labels 
labels <- colnames(basket)

# Create dissimilarity matrix
create.dissim <- function(x){
  
  # scale similarities (turn into correlations)
  x <- cor(x)

  # obtain dissimilarities from scaled similarities
  x <- 1-x 
  
  # make sure diagonal is exactly 0 (not rounded)
  diag(x) <- 0
  
  return(x)
}

# Ordinal transformation of dissimilarities
ord.transform <- function(x) {
  
  # initialize new matrix of ranks
  ord.mat <- matrix(0L, nrow=nrow(x), ncol=ncol(x))
  
  # replace upper and lower matrix with rank values
  ord.mat[lower.tri(ord.mat, diag = FALSE)] <- rank(x[lower.tri(x, diag = FALSE)])
  ord.mat[upper.tri(ord.mat, diag = FALSE)] <- rank(x[upper.tri(x, diag = FALSE)])
  
  # ensure values are integers
  ord.mat <- mapply(ord.mat, FUN=as.integer)
  ord.mat <- matrix(data=ord.mat, ncol=ncol(x), nrow=nrow(x))
  
  return(ord.mat)
}

# obtain dissimilarities
dissim <- create.dissim(as.matrix(basket))

# obtain ordinal dissimilarities (disparities)
disparity <- ord.transform(dissim)


################################################################################
# Predefined functions
################################################################################

# Calculate Euclidean distance
eu.dist <- function(X) {

  # Predefine XX'
  X.tX <- X%*%t(X)
  
  # Calculate Euclidean distance
  c <- diag(X.tX)
  prod <- outer(c,c,"+") - 2*X.tX
  ED <- sqrt(prod)
  
  return(ED)
}

# Obtain B matrix
get.B <- function(x, dissim, ED) {
  
  # Define matrix F
  F <- dissim/ED
  diag(F) <- 0
  
  # Obtain matrix B from F
  F1 <- rowSums(F)
  B <- diag(F1)-F
  
  return(B)
}

# Calculate (normalized) stress
get.stress <- function(dissim, X, tX.X, B, ED) {
  
  n <- nrow(X)

  # Stress can be divided into three elements
  cons <- sum(upper.tri(dissim, diag = FALSE)*(dissim^2))
  eta <- n*sum(diag(tX.X))
  rho <- sum(diag(t(X) %*% B %*% X))
  
  # Calculate raw stress (squared) from the 3 elements
  raw.stress.sq <- cons + eta - 2*rho
  
  # Obtain normalized stress (squared) from raw stress
  norm.stress.sq <- raw.stress.sq/cons
  
  # Take the square root to obtain final stress
  norm.stress <- sqrt(norm.stress.sq)
  
  return(norm.stress)
}

################################################################################
# Smacof algorithm
################################################################################

SMACOF <- function(dissim, config = NULL, eps = 1e-06) {
  
  # Number of variables
  n <- nrow(dissim)
  
  # Use random or initial Torgerson configuration
  if (is.null(config) == T) {
    X <- matrix(data = rnorm(n*2,), nrow = n, ncol = 2)
  } else {
    X <- config
  }
  
  # Set Z to be equal to the initial configuration
  Z <- list(X)
  
  # Get initial euclidean distance
  DZ <- eu.dist(X)
  
  # Get initial B
  BZ <- get.B(X, dissim, DZ)
  
  # Pre-calculate stress value
  stress <- c(0,get.stress(dissim, X, t(X)%*%X, BZ, DZ))
  
  # Initialize counter
  k <- 1
  
  # Start loop
  while (k == 1 | stress[k]-stress[k+1] > eps){
    k <- k+1
    
    # Update X
    X <- (1/n)*BZ%*%Z[[k-1]]
    Z <- list.append(Z,X)
    
    # Obtain new distances
    DZ <- eu.dist(Z[[k]])
    
    # Calculate updated B
    BZ <- get.B(Z[[k]], dissim, DZ)
    
    # Compute updated stress
    stress <- c(stress, get.stress(dissim, Z[[k]], t(Z[[k]])%*%Z[[k]], BZ, DZ))
    
    # Check whether stress goes down
    if (stress[k+1]-stress[k] > 0) {
      cat('\n\nSomething goes wrong! Stress is increasing by ', stress[k+1]-stress[k])
    }
    
  } # End loop
   
  # Output a list containing the results
  output <- list(conf = data.frame(Z[[k]]), stress = stress[k], niter = k)
  
  return(output)
}

################################################################################
# Compare implementation with package
################################################################################

### Ordinary MDS

# obtain initial configeration using classical torgerson
initial <- mds(dissim, ndim = 2, init='torgerson', itmax=1)
initial.conf <- initial$conf

# own implementation
own.result <- SMACOF(dissim, config = initial.conf, eps=1e-06) 
own.result$stress
own.result$conf
own.result$niter

# package results
pack.result <- mds(dissim, init=initial.conf, eps=1e-06)
pack.result$stress
pack.result$conf
pack.result$niter

### Random MDS

# own implementation
own.random.result <- SMACOF(dissim, eps=1e-06)
own.random.result$stress

# package results
pack.random.result <- mds(dissim, init='random', eps=1e-06)
pack.random.result$stress
  
### Ordinal MDS

# package results
pack.ord.result <- mds(dissim, itmax = 1000, init=initial.conf, type = 'ordinal', eps=1e-06)
pack.ord.result$stress
pack.ord.result$conf
pack.ord.result$niter

### Plots

# plot of own results
own.config <- data.frame(own.result$conf)
rownames(own.config) <- labels
own.plot <- ggplot(data = own.config, aes(x = D1, y = D2))+
  geom_point(size = 1.5) + 
  ggtitle('Multidimensional scaling: product co-purchases (own implementation)') + 
  ylab("") + xlab("") + theme_bw() + 
  coord_cartesian(ylim = c(-1, 1), xlim = c(-1, 1)) + 
  geom_text(aes(label=labels),hjust=0.3, vjust=1.3, size=4)
print(own.plot)
ggsave(own.plot,filename='Plots/own_plot.png',width=8, height=8)

# plot of package results
pack.config <- data.frame(pack.result$conf)
rownames(pack.config) <- labels
pack.plot <- ggplot(data = pack.config, aes(x = D1, y = D2))+
  geom_point(size = 1.5) + ggtitle('Multidimensional scaling: product co-purchases (package)') + 
  ylab("") + xlab("") + theme_bw() + 
  coord_cartesian(ylim = c(-1, 1), xlim = c(-1, 1)) + 
  geom_text(aes(label=labels),hjust=0.3, vjust=1.3, size=4)
print(pack.plot)
ggsave(pack.plot,filename='Plots/package_plot.png',width=8, height=8)

# comparison plot of own implementation and package results


# plot of ordinal results
pack.ord.config <- data.frame(pack.ord.result$conf)
rownames(pack.ord.config) <- labels
pack.plot <- ggplot(data = pack.ord.config, aes(x = D1, y = D2))+
  geom_point(size = 1.5) + ggtitle('Multidimensional scaling (ordinal transformation): product co-purchases') + 
  ylab("") + xlab("") + theme_bw() + 
  coord_cartesian(ylim = c(-1, 1), xlim = c(-1, 1)) + 
  geom_text(aes(label=labels),hjust=0.3, vjust=1.3, size=4)
print(pack.plot)
ggsave(pack.plot,filename='Plots/package_ord_plot.png',width=8, height=8)
