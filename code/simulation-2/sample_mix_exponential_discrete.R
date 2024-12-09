#-------------------------------------------------------------------
# Generate data from a mixture of exponential-discrete distributions
#-------------------------------------------------------------------
mixture_exp_distr <- function(to, by, lambda1, lambda1_2, lambda2, lambda2_2, n){
  # Sequence generation
  t = seq(0, to, by)
  ### Exponential mixture
  # Non-normalized data function
  f <- function(t){
    lambda1*exp(-lambda1_2*t) - lambda2*exp(-lambda2_2*t)
  }
  # normalization constant
  cons_norm <- integrate(f,0,to)$value
  # PDF
  f_norm <- function(t){
    (lambda1*exp(-lambda1_2*t) - lambda2*exp(-lambda2_2*t))/cons_norm
  }
  # CDF
  CDF <- function(t){
    integrate(f_norm, 0,t)$value
  }
  CDF <- Vectorize(CDF)
  # Generate a random sample from an Uniform distribution
  u <- runif(n, min=0, max=1)
  # Inverse CDF function
  CDF_inv <- function(u){
    uniroot(function(t){CDF(t)-u},
            interval=c(0,to))$root
  }
  CDF_inv <- Vectorize(CDF_inv)
  return(CDF_inv(u))
}
#---------------------------------------------
# 3 CLUSTERS
#---------------------------------------------
generate_mix_exp_disc_k3 <- function(size, con_var, cat_var, bin_var, mix_weights) {
  
  # libraries
  library(dplyr)
  library(Rlab)
  
  # initialize matrix
  con_data <- matrix(0, ncol = con_var, nrow = size)
  cat_data <- matrix(0, ncol = cat_var, nrow = size)
  bin_data <- matrix(0, ncol = bin_var, nrow = size)
  true_id <- integer(size)
  
  # Define the categories and probabilities
  prob1 <- c(0.25, 0.2, 0.15, 0.1, 0.1, 0.05, 0.03, 0.02, 0.02, 0.02, 0.02, 0.02, 0.01, 0.01)
  prob2 <- c(0.01, 0.01, 0.02, 0.02, 0.02, 0.02, 0.02, 0.03, 0.05, 0.1, 0.1, 0.15, 0.2, 0.25)
  prob3 <- c(0.01, 0.02, 0.02, 0.05, 0.1, 0.15, 0.25, 0.15, 0.1, 0.1, 0.02, 0.01, 0.01, 0.01)
  
  for (s in 1:size) {
    
    pr <- mix_weights
    # pr <- rep(1/3, times = 3)
    b <- which.max(rmultinom(1, size = 1, prob = pr))
    
    for (con in 1:con_var) {
      # print(paste("con:", con))
      
      for (cat in 1:cat_var) {
        # print(paste("cat:", cat))
        
        for (bin in 1:bin_var) {
          # print(paste("bin:", bin))
          
          if (b == 1) {
            # Continuous variable: mixture of AH exponential distributions
            x <- mixture_exp_distr(to = 4000, by = 0.1, lambda1 = 3.5, lambda1_2 = 0.05,
                                     lambda2 = 20, lambda2_2 = 1, n = 1)
            
            con_data[s, con] <- x
            
            # Categorical variable
            r <- rmultinom(1, 1, prob1)
            cat_data[s, cat] <- which(r==1)
            
            # Binary variable
            bin_data[s, bin] <- rbinom(1, 1, 0.8)
            
            # True label
            t_id <- b
          } else if (b == 2) {
            x <- mixture_exp_distr(to = 4000, by = 0.1, lambda1 = 1.3, lambda1_2 = 0.05,
                                     lambda2 = 20, lambda2_2 = 1, n = 1)
            con_data[s, con] <- x
            
            r <- rmultinom(1, 1, prob2)
            cat_data[s, cat] <- which(r==1)
            
            bin_data[s, bin] <- rbinom(1, 1, 0.5)
            t_id <- b
          } else {
            x <- mixture_exp_distr(to = 4000, by = 0.1, lambda1 = 1.1, lambda1_2 = 0.05,
                                     lambda2 = 20, lambda2_2 = 1, n = 1)
            con_data[s, con] <- x
            
            r <- rmultinom(1, 1, prob3)
            cat_data[s, cat] <- which(r==1)
            
            bin_data[s, bin] <- rbinom(1, 1, 0.2)
            t_id <- b
          }
        }
        
      }
      
    }
    true_id[s] <- t_id
  }
  con_names <- sprintf("con%d",1:con_var)
  colnames(con_data) <- con_names
  cat_names <- sprintf("cat%d",1:cat_var)
  colnames(cat_data) <- cat_names
  mode(cat_data) <- "character"
  bin_names <- sprintf("bin%d",1:bin_var)
  colnames(bin_data) <- bin_names
  mode(bin_data) <- "character"
  # data frame
  df <- data.frame(con_data, as.data.frame(cat_data, stringsAsFactors = TRUE),
                   as.data.frame(bin_data, stringsAsFactors = TRUE),
                   true_id = as.factor(true_id))
}

#---------------------------------------------
# 2 CLUSTERS
#---------------------------------------------
generate_mix_exp_disc_k2 <- function(size, mix_weigts, con_var, cat_var, bin_var) {
  
  # libraries
  library(dplyr)
  library(Rlab)
  
  # initialize matrix
  con_data <- matrix(0, ncol = con_var, nrow = size)
  cat_data <- matrix(0, ncol = cat_var, nrow = size)
  bin_data <- matrix(0, ncol = bin_var, nrow = size)
  true_id <- integer(size)
  
  # Define the categories and probabilities
  prob1 <- c(9524,316,252,548,326,428,316,192,232,1134,1888,1542,982,1250)
  prob2 <- c(512,129,519,876,352,230,777,304,53,1063,85,1312,216,295)
  
  for (s in 1:size) {
    b <- rbern(1, prob = mix_weigts)
    
    for (con in 1:con_var) {
      # print(paste("con:", con))
      
      for (cat in 1:cat_var) {
        # print(paste("cat:", cat))
        
        for (bin in 1:bin_var) {
          # print(paste("bin:", bin))
          
          if (b == 1) {
            # Continuous variable: mixture of AH exponential distributions
            x <- mixture_exp_distr(to = 4000, by = 0.1, lambda1 = 1.3, lambda1_2 = 0.05,
                                     lambda2 = 20, lambda2_2 = 1, n = 1)
            
            con_data[s, con] <- x
            
            # Categorical variable
            r <- rmultinom(1, 1, prob1)
            cat_data[s, cat] <- which(r==1)
            
            # Binary variable
            bin_data[s, bin] <- rbinom(1, 1, 0.64)
            
            # True label
            t_id <- b
          } else {
            x <- mixture_exp_distr(to = 8000, by = 0.1, lambda1 = 1.1, lambda1_2 = 0.05,
                                     lambda2 = 20, lambda2_2 = 1, n = 1)
            con_data[s, con] <- x
            
            r <- rmultinom(1, 1, prob2)
            cat_data[s, cat] <- which(r==1)
            
            bin_data[s, bin] <- rbinom(1, 1, 0.3)
            t_id <- b
          }
        }
        
      }
      
    }
    true_id[s] <- t_id
  }
  con_names <- sprintf("con%d",1:con_var)
  colnames(con_data) <- con_names
  cat_names <- sprintf("cat%d",1:cat_var)
  colnames(cat_data) <- cat_names
  mode(cat_data) <- "character"
  bin_names <- sprintf("bin%d",1:bin_var)
  colnames(bin_data) <- bin_names
  mode(bin_data) <- "character"
  # data frame
  df <- data.frame(con_data, as.data.frame(cat_data, stringsAsFactors = TRUE),
                   as.data.frame(bin_data, stringsAsFactors = TRUE),
                   true_id = as.factor(true_id))
}
