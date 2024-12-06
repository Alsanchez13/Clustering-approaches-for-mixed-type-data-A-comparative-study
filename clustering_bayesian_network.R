library(bnlearn)
library(magrittr)
library(dplyr)

#-------------------------------------------
# learning Bayesian Netoworks for clustreing
#-------------------------------------------

split_dataframe <- function(df, k) {
  
  prob_vector <- rep(1/k, k)
  
  labels <- sample(x = seq_along(prob_vector), size = nrow(df), replace = TRUE, prob = prob_vector)
  
  # Add a new column with the labels to the dataframe
  df$class <- as.factor(labels)
  
  return(df)
}


clust_bn <- function(df, K, score_bn, multi_start) {
  
  # START CLUSTERING BN
  multi_labels <- matrix(NA, nrow = nrow(df), ncol = multi_start)
  multi_score <- numeric(multi_start)
  epsilon <- 10e-8
  
  for (multi in 1:multi_start) {
    # cat("Multi:", multi, "\n")
    
    stoch_score <- numeric(200)
    stoch_labels <- matrix(NA, nrow = nrow(df), ncol = 200)
    scores <- c()
    
    for (st in 1:200) {
      # cat("stochastic:", st, "\n")
      
      df <- split_dataframe(df, k=K)
      
      # E-STEP
      pp_x <- matrix(NA, nrow = nrow(df), ncol = K)
      
      for (e in 1:K) {
        # restricted df
        df_k <- df[(df$class==as.character(e)),-ncol(df)]
        # learn network structure
        learn_structure <- bnlearn::tabu(df_k, score = score_bn)
        # fit the model
        fitted <- bnlearn::bn.fit(learn_structure, df_k, replace.unidentifiable=TRUE)
        # structure plot 
        # bnlearn::graphviz.plot(fitted)
        # compute posterior probabilities
        post_prob_x <- exp(logLik(fitted, df[,-ncol(df)], by.sample = TRUE))
        # store posterior probabilities
        pp_x[, e] <- post_prob_x
        
      }
      
      # S-STEP
      #  get new Partition P1
      pp_x <- pp_x/rowSums(pp_x)
      label <- character(nrow(df))
      
      for (s in 1:nrow(df)) {
        label[s] <- sample.int(K, 1, prob = pp_x[s, ])
      }
      
      # labels <- cbind(labels, label)
      df$class <- as.factor(label)
      
      # M-STEP
      loglikelihood <- 0
      for (m in 1:K) {
        # restricted df
        df_k <- df[(df$class==as.character(m)),-ncol(df)]
        # learn network structure
        learn_structure <- bnlearn::tabu(df_k, score = score_bn)
        # fit the model
        fitted <- bnlearn::bn.fit(learn_structure, df_k, replace.unidentifiable=TRUE)
        # structure plot 
        # bnlearn::graphviz.plot(fitted)
        # compute posterior probabilities
        ll <- sum(logLik(fitted, df_k) + nrow(df_k)*log(nrow(df_k)/nrow(df)))
        # store posterior probabilities
        loglikelihood <- loglikelihood + ll
        
      }
      stoch_score[st] <- loglikelihood
      stoch_labels[,st] <- label
      
    }
    best_score <- max(stoch_score)
    scores <- c(best_score)
    best_label <- as.factor(stoch_labels[,which.max(stoch_score)])
    df$class <- best_label
    
    # CEM starting from the best stochastic partition
    for (i in 1:800) {
      # cat("cem:", i, "\n")
      
      # E-STEP
      pp_x <- matrix(NA, nrow = nrow(df), ncol = K)
      
      for (e in 1:K) {
        # restricted df
        df_k <- df[(df$class==as.character(e)),-ncol(df)]
        # learn network structure
        learn_structure <- bnlearn::tabu(df_k, score = score_bn)
        # fit the model
        fitted <- bnlearn::bn.fit(learn_structure, df_k, replace.unidentifiable=TRUE)
        # structure plot 
        # bnlearn::graphviz.plot(fitted)
        # compute posterior probabilities
        post_prob_x <- exp(logLik(fitted, df[,-ncol(df)], by.sample = TRUE))
        # store posterior probabilities046678
        pp_x[, e] <- post_prob_x
        
      }
      
      # C-STEP
      #  get new Partition P1
      label <- as.factor(apply(pp_x, 1, which.max))
      df$class <- label
      
      # M-STEP
      loglikelihood <- 0
      for (m in 1:K) {
        # restricted df
        df_k <- df[(df$class==as.character(m)),-ncol(df)]
        # learn network structure
        learn_structure <- bnlearn::tabu(df_k, score = score_bn)
        # fit the model
        fitted <- bnlearn::bn.fit(learn_structure, df_k, replace.unidentifiable=TRUE)
        # structure plot 
        # bnlearn::graphviz.plot(fitted)
        # compute posterior probabilities
        ll <- sum(logLik(fitted, df_k) + nrow(df_k)*log(nrow(df_k)/nrow(df)))
        # store posterior probabilities
        loglikelihood <- loglikelihood + ll
        
      }
      scores <- c(scores, loglikelihood)
      # Check if significant change in likelihood
      if (abs(scores[i] - loglikelihood) < epsilon) break
      
    }
    multi_labels[,multi] <- df$class
    multi_score[multi] <- max(scores)
    
  }
  df$class <- as.factor(multi_labels[,which.max(multi_score)])
  
  result <- list(
    class = df$class,
    loglikelihood = multi_score
    
  )
  return(result)

}
