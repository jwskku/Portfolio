# 1/13 수정

# Load necessary libraries


# Marginal Feature Screening (FSM) with parallelization
marginal_screening <- function(X, Y, d) {
    
  #parallelization setting
  library(energy)
  library(foreach)
  library(doParallel)
  cores <- detectCores()/2  
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  
  n <- ncol(X)
  correlations <- foreach(j = 1:n, .combine = 'c', .packages = 'energy') %dopar% {
    dcor(X[, j], Y)
  }
  m_order <- order(correlations, decreasing = TRUE)
  selected <- m_order[1:d]

  stopCluster(cl)
  
  list(selected = selected, correlations = correlations,m_order = m_order)
}

# Two-Stage Sufficient Variable Selection (SVS2) with parallelization
two_stage_svs_parallel <- function(X, Y, d1, d2, S = 2) {
  
  n <- length(Y)
  p <- ncol(X)
  
  # Step 1: Marginal screening
  marginal <- marginal_screening(X, Y, d1)
  m_order <- marginal$m_order
  marginal_indices <- marginal$selected
  
  # Step 2: Check if Y is continuous or discrete
  if (is.numeric(Y) && length(unique(Y)) > 5) {
    # Continuous data - slicing
    slices <- cut(Y, breaks = quantile(Y, probs = seq(0, 1, length.out = S + 1)),
                  include.lowest = TRUE)
  } else {
    # Discrete data - use unique values as slices
    slices <- as.factor(Y)
  }
  
  cores <- detectCores()/2  
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  # Parallel Conditional independence test
  u_double_star <- foreach(i = 1:p, .combine = 'c', .packages = c('energy')) %dopar% {
    sum(sapply(levels(slices), function(s) {
      slice_indices <- which(slices == s)
      dcor(X[slice_indices, i, drop = FALSE], X[slice_indices, -i, drop = FALSE])
    }))
  }
  
  stopCluster(cl)
  
  c_order <- order(u_double_star, decreasing = TRUE)
  condition_indices <- setdiff(c_order, marginal_indices)[1:d2]
  
  # Final selected set
  final_selected <- union(marginal_indices, condition_indices)
  list(selected = final_selected, marginal_indices = marginal_indices,
       condition_indices = condition_indices , m_order = m_order , c_order = c_order)
  
}


