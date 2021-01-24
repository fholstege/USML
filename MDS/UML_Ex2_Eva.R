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

# Create dissimilarity matrix
dissim <- 1-basket/diag(as.matrix(basket))

# labels 
labels <- colnames(basket)

################################################################################
# Predefined functions
################################################################################

# Calculate Euclidean distance
dist <- function(X) {
   
  n <- nrow(X)
  ED <- matrix(, n, n)
  X.tX <- X%*%t(X)
  c <- diag(X.tX)
  prod <- outer(c,c,"+") - 2*X.tX
  for(i in 1:n){
    for(j in 1:n){
      ED[i,j] = sqrt(prod[i,j])
    }
  }
  return(ED)
}

# Obtain lower triangular sum
lower.tri.sum <- function(x) {
  return(sum(lower.tri(x, diag = FALSE)))
}

# Obtain B matrix
get.B <- function(x) {
  n <- nrow(x)
  F <- matrix(0L, nrow = n, ncol = n)
  ED <- dist(x)
  for(i in 1:n){
    for(j in 1:n){
      if (ED[i,j]!=0){
        F[i,j] <- dissim[i,j]/ED[i,j]
      }
      else {
        F[i,j] <- 0
      }
      }
  }
  F1 <- rowSums(F)
  B <- diag(F1)-F
  return(B)
}

# Get squared elements of a matrix
squared.mat <- function(x) {
  sq.mat <- matrix(, nrow = nrow(x), ncol = ncol(x))
  for(i in 1:nrow(x)){
    for (j in 1:ncol(x)){
      sq.mat[i,j] <- x[i,j]^2
    }
  }
  return(sq.mat)
}

# Calculate (normalized) stress value
get.stress <- function(dissim, X, tX.X, B) {
  n <- nrow(X)
  cons <- lower.tri.sum(squared.mat(dissim))
  eta <- n*sum(diag(tX.X))
  rho <- sum(diag(t(X) %*% B %*% X))
  raw.stress <- cons + eta - 2*rho
  norm.stress <- raw.stress/cons
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
  print(X)
  
  # Set Z to be equal to the initial configuration
  Z <- list(X)
  
  # Get initial euclidean distance
  DZ <- dist(X)
  
  # Get initial B
  BZ <- get.B(X)
  
  # Pre-calculate stress value
  stress <- c(0)
  stress <- c(stress,get.stress(dissim, X, t(X)%*%X, BZ))
  
  # Initialize counter
  k <- 1
  
  # Start loop
  while (k == 1 | stress[k]-stress[k+1] > eps){
    k <- k+1
    
    # Update X
    X <- (1/n)*BZ%*%Z[[k-1]]
    Z <- list.append(Z,X)
    
    # Obtain new distances
    DZ <- dist(Z[[k]])
    
    # Calculate updated B
    BZ <- get.B(Z[[k]])
    
    # Compute updated stress
    stress <- c(stress, get.stress(dissim, Z[[k]], t(Z[[k]])%*%Z[[k]], BZ))
  
    # Check whether stress goes down
    if (stress[k+1]-stress[k] > 0) {
      cat('\n\nSomething goes wrong! Stress is increasing by ', stress[k+1]-stress[k])
    }
    
  } # End loop
   
  cat('\n\nIterations:', k) # Total number of iterations
  
  # Output a list containing the results
  cat('\n\nFinal stress:', stress[k])
  cat('\n\nFinal configuration:\n')
  print(data.frame(Z[[k]]))
  output <- list(conf = data.frame(Z[[k]]), stress = stress[k], niter = k)
  
  return(output)
}


################################################################################
# Compare implementation with package
################################################################################

### Ordinary MDS

# obtain initial configeration using classical torgerson
initial <- mds(dissim, init='torgerson', itmax=1)
initial.conf <- initial$conf

# own implementation
own.result <- SMACOF(dissim, config = initial.conf, eps=1e-06) 
own.result$stress
own.result$conf
own.result$niter

# package results
pack.result <- mds(dissim, itmax = 1000, init=initial.conf, eps=1e-06)
pack.result$stress
pack.result$conf
pack.result$niter

# plot of package results
pack.config <- data.frame(pack.result$conf)
rownames(pack.config) <- labels
pack.plot <- ggplot(data = pack.config, aes(x = D1, y = D2))+
  geom_point(size = 1.5) + ggtitle('MDS co-purchases') + 
  ylab("") + xlab("") + theme_bw() + 
  coord_cartesian(ylim = c(-1, 1), xlim = c(-1, 1)) + 
  geom_text(aes(label=labels),hjust=0.3, vjust=1.3, size=4)
print(pack.plot)
ggsave(pack.plot,filename='Plots/package_plot.png',width=8, height=8)

### Ordinal MDS

# own function

# package results
pack.ord.result <- mds(df, itmax = 1000, init='torgerson', type = 'ordinal', eps=1e-06)
pack.ord.result$stress
pack.ord.result$conf
pack.ord.result$niter

# plot of package results
pack.ord.config <- data.frame(pack.ord.result$conf)
rownames(pack.ord.config) <- labels
pack.plot <- ggplot(data = pack.ord.config, aes(x = D1, y = D2))+
  geom_point(size = 1.5) + ggtitle('MDS co-purchases (ordinal transformation)') + 
  ylab("") + xlab("") + theme_bw() + 
  coord_cartesian(ylim = c(-1, 1), xlim = c(-1, 1)) + 
  geom_text(aes(label=labels),hjust=0.3, vjust=1.3, size=4)
print(pack.plot)
ggsave(pack.plot,filename='Plots/package_ord_plot.png',width=8, height=8)

