rfe_pseudo_samples_ptest_break <- function(X, Y, c_star, cost = 1, gamma = 1, graph = TRUE) {
  library(e1071)        # For SVM
  library(matrixStats)  # For MAD
  library(scales)       # For normalization
  library(doParallel)   # For parallel processing
  library(foreach)      # For parallel loops
  library(ggplot2)      # For plotting
  
  # Normalize the data
  X_norm <- as.data.frame(scale(X, center = TRUE, scale = TRUE))
  p_star <- ncol(X_norm)  # Initial number of features
  p <- p_star
  Ranks <- vector(mode = "character", length = p_star)  # To store removal order
  svmlist <- list()  # To store SVM models with associated features
  selected <- c()    # To record selected (ineffective) variables
  
  # Set up parallel backend
  cl <- makeCluster(detectCores() / 2)
  registerDoParallel(cl)
  
  while (p >= 2) {
    # Train SVM with RBF kernel on the current set of features
    gamma_i <- gamma / ncol(X_norm)
    SVM_current <- svm(X_norm, Y, kernel = "radial", scale = FALSE, cost = cost, gamma = gamma_i)
    
    # Store the current model and corresponding feature set
    svmlist[[p]] <- list(model = SVM_current, features = colnames(X_norm))
    
    # Compute MAD per variable for current features
    rank_criteria <- foreach(i = 1:p, .combine = c, .packages = c('e1071', 'matrixStats')) %dopar% {
      n <- nrow(X_norm)
      X_i <- X_norm[, i]
      min_val <- min(X_i)
      max_val <- max(X_i)
      quantiles <- seq(min_val, max_val, length.out = c_star)
      
      pseudo_i <- matrix(0, nrow = c_star, ncol = p)
      pseudo_i[, i] <- quantiles
      
      predictions <- predict(SVM_current, pseudo_i, decision.values = TRUE)
      decision_vals <- attr(predictions, "decision.values")
      
      mad(decision_vals)
    }
    
    # Order variables by increasing MAD (lower MAD implies a candidate for removal)
    order_idx <- order(rank_criteria)
    
    candidate_removed <- FALSE
    # Iterate over candidates (ordered by MAD) and apply the permutation test
    for (idx in order_idx) {
      candidate <- colnames(X_norm)[idx]
      # Use the current model and data to test the candidate variable
      result <- per_test_candidate(rep = 2000, model = SVM_current, X = X_norm, Y = Y, candidate = candidate, graph = graph, cost = cost, gamma = gamma)
      
      if (result == FALSE) {
        Ranks[p] <- candidate
        X_norm <- X_norm[, setdiff(colnames(X_norm), candidate), drop = FALSE]
        p <- p - 1
        candidate_removed <- TRUE
        break  # Break out of the candidate loop and retrain with the updated feature set
      }
    }
    
    # If no candidate removed, stop the iteration
    if (!candidate_removed) {
      selected <- colnames(X_norm)
#      print('Remaining variables expected to have effect!')
      break
    }
  }

  stopCluster(cl)
  
  # Return both the final ranking and list of SVM models with feature sets
  return(list(rank = Ranks, svm_models = svmlist, selected_variables = selected))
}



# New imputation test function that tests the current candidate variable directly
per_test_candidate <- function(rep, model, X, Y, candidate, graph = TRUE, cost = cost, gamma = gamma) { 
  
  library(foreach)
  gamma_i <- gamma / ncol(X[,setdiff(colnames(X), candidate)])
  # Calculate baseline error from the current model
  m_svm <- svm(X[,setdiff(colnames(X), candidate)], Y, kernel = "radial", scale = FALSE, cost = cost, gamma = gamma_i)
  m_error <- mean(m_svm$fitted != Y)  

  # Permutation test: permute candidate column and compute misclassification errors
  error_list <- foreach(j = 1:rep, .combine = c, .packages = c('e1071')) %dopar% {
    X_perm <- X
    X_perm[, candidate] <- sample(X_perm[, candidate])
    pred_perm <- predict(model, X_perm)
    mean(pred_perm != Y)
  }
  
  can_error <- mean(model$fitted != Y)
  
  # Calculate differences between permuted errors and the baseline error
  diff_list <- m_error - error_list
  diff_point <- m_error - can_error
  
  # Use the 99th percentile of the differences as the threshold
  quantile_val <- quantile(diff_list, probs = 1)
  
  # Optionally create and print a histogram of the error differences
  if (graph == TRUE) {
    pl <- ggplot(data.frame(diff_list = diff_list), aes(x = diff_list)) +
      geom_histogram(fill = "lightblue", color = "black", bins = 30) +
      geom_vline(aes(xintercept = diff_point, color = "Point"), linewidth = 1) +
      geom_vline(aes(xintercept = quantile_val, color = "Quantile"), linetype = "solid", linewidth = 0.5) +
      geom_vline(aes(xintercept = 0, color = "Quantile"), linetype = "solid", linewidth = 0.5) +
      scale_color_manual(name = "Line Type", values = c("Point" = "red", "Quantile" = "black"), 
                         labels = c("Point (Observed)", "Quantile (100%)")) +
      labs(title = paste0("Histogram of error differences for candidate: ", candidate),
           x = paste0("Error difference for ", candidate),
           y = "Count") +
      theme_minimal() +
      theme(legend.position = "bottom")
    print(pl)
  }
  
  # Decision rule: if the average increase in error upon permutation is sufficiently high,
  # then we conclude the candidate is ineffective (and return TRUE)
  if (diff_point > quantile_val & diff_point > 0) {
    return(TRUE)
  } else { 
    return(FALSE)
  }
}
