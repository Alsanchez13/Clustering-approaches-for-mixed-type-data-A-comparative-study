library(bnlearn)

# Generate random probability vectors
generate_prob_vector <- function(K, n, seed = NULL) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  result <- c()
  
  for (i in 1:K) {
    # Generate random numbers
    # Set seed for reproducibility if provided
    prob_vector <- sample(seq(0.05, 1, by = 0.4), size = n, replace = FALSE)
    
    # Normalize to sum up to 1
    prob_vector <- prob_vector / sum(prob_vector)
    
    result <- cbind(result, prob_vector)
    
  }
  results_vector <- as.vector(result)
  return(results_vector)
}


# sample data from an Bayesian Network with SBN structure
sim_model_sbn_bn <- function(K, seed = NULL) {
  
  # Set seed for reproducibility if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # network specification.
  dag = bnlearn::model2network("[C][X6|C][X1|C][X2|X1:C][X4|X2:C][X5|X2:C][X3|X1:C]")

  # Build BN
  L <- c("a", "b", "c")
  
  bn =  bnlearn::custom.fit(dag,
                            list(
    C = matrix(rep(1/K, K), ncol = K, dimnames = list(NULL, c(sprintf("%d",1:K)))),
    X1 = matrix(generate_prob_vector(K=K, n=3, seed = seed), ncol = K,
                dimnames = list(X1 = L, C = c(sprintf("%d",1:K)))),
    X2 = array(generate_prob_vector(K=3*K, n=3, seed = seed), dim = c(3, 3, K),
               dimnames = list(X2 = L, X1 = L, C = c(sprintf("%d",1:K)))),
    X3 = list(coef = matrix(seq(0.5, 1.5*K, by = 0.5), ncol = 3*K,
                            dimnames = list("(Intercept)", NULL)),
              sd = sample(seq(0.5, 1.5, by = 0.1), size = 3*K, replace = TRUE)),
    X4 = list(coef = matrix(seq(2*K, (2*K + 3*K)-1, by = 1), ncol = 3*K,
                            dimnames = list(X4 ="(Intercept)", NULL)),
              sd = sample(seq(0.5, 1.5, by = 0.1), size = 3*K, replace = TRUE)),
    X5 = array(generate_prob_vector(K=3*K, n=2, seed = seed), dim = c(2, 3, K),
               dimnames = list(X5 = c("a", "b"), X2 = L, C = c(sprintf("%d",1:K)))),
    X6 = list(coef = matrix(seq((2*K + 3*K)-1, (2*K + 3*K)-2 + K, by = 1), ncol = K,
                            dimnames = list(X6 = "(Intercept)", NULL)),
              sd = sample(seq(0.5, 1.5, by = 0.1), size = K, replace = TRUE))
  )
  )
  
  output <- list(
    dag = dag,
    model = bn
  )
  
  return(output)
}
