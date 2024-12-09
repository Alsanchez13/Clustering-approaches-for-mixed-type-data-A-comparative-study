library(bnlearn)
library(dplyr)

generate_single_prob_vector <- function(C, n, seed = NULL) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  result <- c()
  
  for (i in 1:C) {
    # Generate random numbers
    prob_vector <- sample(seq(0.05, 1, by = 0.4), size = n, replace = FALSE)
    
    # Normalize to sum up to 1
    prob_vector <- prob_vector / sum(prob_vector)
    
    result <- cbind(result, prob_vector)
    
  }
  results_vector <- as.vector(result)
  return(results_vector)
}

# Generate sample data from a mixture of BN with the same DAG structure
sim_model_mixture_bn <- function(K, seed = NULL) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # network specification.
  dag = bnlearn::model2network("[X6][X1][X2|X1][X4|X2][X5|X2][X3|X1]")
  
  # Build BN
  L <- c("a", "b", "c")
  
  # Generate mixture coefficients
  pr <- rep(1/K, times = K)

  # initialize empty lists to store results
  list_bn <- list()
  
  for (i in 1:K) {
    
    bn =  bnlearn::custom.fit(dag,
                              list(
      X1 = array(generate_single_prob_vector(C=1, n=3, seed = seed), dim=c(1, 3),
                  dimnames = list(NULL, X1 = L)),
      X2 = array(generate_single_prob_vector(C=3, n=3, seed = seed), dim = c(3, 3),
                 dimnames = list(X2 = L, X1 = L)),
      X3 = list(coef = matrix(seq(0+2*i, 2+2*i, by = 1), ncol = 3,
                              dimnames = list("(Intercept)", NULL)),
                sd = sample(seq(0.5, 1.5, by = 0.1), size = 3, replace = TRUE)),
      X4 = list(coef = matrix(seq(5+2*i, 7+2*i, by = 1), ncol = 3,
                              dimnames = list(X4 ="(Intercept)", NULL)),
                sd = sample(seq(0.5, 1.5, by = 0.1), size = 3, replace = TRUE)),
      X5 = array(generate_single_prob_vector(C=3, n=2, seed = seed), dim = c(2, 3),
                 dimnames = list(X5 = c("a", "b"), X2 = L)),
      X6 = list(coef = c("(Intercept)" = seq(9+2*i, 9+2*i, by = 1)),
                sd = sample(seq(0.5, 1.5, by = 0.1), size = 1, replace = TRUE))
    ))
    
    list_bn[[i]] <- bn

  }
  
  output <- list(
    dag = dag,
    model = list_bn,
    mix.coeff = pr
  )
  
  return(output)
}


sample_data_sim_model_mixture_bn <- function(model, size) {
  
  # Generate labels according to mixing coefficients probabilities
  b <- rmultinom(size, size = 1, prob = model$mix.coeff)
  
  # Get labels
  labels <- apply(b, 2, function(x) which(x == 1, arr.ind = TRUE))
  
  # initialize empty lists to store results
  list_data <- list()
  
  for (i in 1:length(model$mix.coeff)) {
    
    data <-  bnlearn::rbn(model$model[[i]], sum(labels == i))
    data$label <- as.factor(i)
    list_data[[i]] <- data
  }
  
  final_data <- dplyr::bind_rows(list_data)
  return(final_data)
}
