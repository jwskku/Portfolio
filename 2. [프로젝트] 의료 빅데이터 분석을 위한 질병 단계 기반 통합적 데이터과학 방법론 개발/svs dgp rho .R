
# Load required libraries
library(MASS)  # For multivariate normal data generation
library(energy)  # For distance correlation calculation
library(ggplot2)
library(gridExtra)
library(reshape2)
library(e1071)
library(SIS)
library(patchwork)
library(progress)
library(GGally)
library(plotly)
library(mt)
library(matrixStats) 
library(Matrix)
library(plotrix)
library(sparseSVM)
library(doParallel)   
library(foreach)      
library(kernlab)
library(caret)



set.seed(1)

# Simulation parameters
# n <- 200
# p <- 1000
# rho <- 0.7  # Different correlation values
# beta <- c(5,5,5, -15 * sqrt(rho), rep(0, p - 4))  # Beta vector for active predictors

# Generate covariance matrix
generate_cov_matrix <- function(p, rho) {
  Sigma <- matrix(rho, nrow = p, ncol = p)
  diag(Sigma) <- 1
  Sigma[4, ] <- Sigma[, 4] <- sqrt(rho)
  Sigma[4, 4] <- 1
  return(Sigma)
}
# 
# generate_data1 <- function(n, p, rho) {
#   Sigma <- generate_cov_matrix(p, rho)
#   X <- mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
#   colnames(X) <- paste0("x", 1:p)
#   eta <- X %*% beta
#   p_i <- exp(eta) / (1+ exp(eta))
#   Y <- rbinom(n = n, size = 1, prob = p_i)
#   Y <- as.factor(Y)
#   return(list(X = X, Y = Y))
# }

# data <- generate_data1(n = n,p = p,rho = rho)
# 
# summary(data$Y)
# ggpairs(data$X, columns = c("x1", "x2", "x3","x4"), aes(color = as.factor(data$Y))) +
#   theme_minimal()

########### setting 2
# generate_cov_matrix <- function(p, rho) {
#   Sigma <- matrix(rho, nrow = p, ncol = p)
#   Sigma[4, ] <- Sigma[, 4] <- sqrt(rho)
#   Sigma[4, 1] <- Sigma[1, 4] <- 0
#   diag(Sigma) <- 1
#   
#   ev <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
#   epsilon <- 1e-6
#   if(min(ev) <= 0){
#     shift <- abs(min(ev)) + epsilon
#     Sigma_pd <- Sigma + shift * diag(p)
#         cat("Added diagonal shift:", shift, "\n")
#   } else {
#     Sigma_pd <- Sigma
#   }
#   return(Sigma_pd)
# }


# set.seed(1)
# n <- 300
# p <- 1000
# rho <- 0.7  
# b1 <- 2.5
# intercept <- 0.5
# beta1 <- 7
# beta2 <- 3
# beta3 <- 6
# gam <- b1*exp(-1*(b1^2/2))*beta3  
# beta4 <- -sqrt(rho) * gam 
# beta_v <- c(beta1, beta2, beta3, beta4, rep(0, p - 4))


generate_data2 <- function(n, p, rho) {
  Sigma <- generate_cov_matrix(p = p,rho = rho)
  X <- mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  colnames(X) <- paste0("x", 1:p)
  
  # Define the intercept
  intercept <- intercept
  b1 <- b1
  
  # Construct the linear predictor
  eta <- intercept + 
    beta_v[1] * X[,1] + 
    beta_v[2] * (X[,2]^2 - 1) + 
    beta_v[3] * sin(b1 * X[,3]) + 
    beta_v[4] * X[,4]
  
  p_i <- exp(eta) / (1 + exp(eta))
  
#  hist(p_i)  # Check the distribution of pi
  
  Y <- rbinom(n = n, size = 1, prob = p_i)
  Y <- as.factor(Y)
  
  return(list(X = X, Y = Y))
}

# data <- generate_data2(n = n,p = p,rho = rho)
# summary(data$Y)






############# setting 3
generate_cov_matrix2 <- function(p, rho) {
  Sigma <- matrix(rho, nrow = p, ncol = p)
  Sigma[4, ] <- Sigma[, 4] <- sqrt(rho)
  Sigma[7, ] <- Sigma[, 7] <- sqrt(rho)
  #Sigma[4, 5] <- Sigma[5, 4] <- Sigma[4, 6] <- Sigma[6, 4] <- Sigma[4, 7] <- Sigma[7, 4] <- 0
  Sigma[4, 6] <- Sigma[6, 4] <- Sigma[4, 7] <- Sigma[7, 4] <- 0
  #Sigma[1, 7] <- Sigma[7, 1] <- Sigma[7, 3] <- Sigma[3, 7] <- Sigma[7, 4] <- Sigma[4, 7] <- 0
  Sigma[1, 7] <- Sigma[7, 1] <- Sigma[7, 4] <- Sigma[4, 7] <- 0
  diag(Sigma) <- 1
  
  ev <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
  epsilon <- 1e-6
  if(min(ev) <= 0){
    shift <- abs(min(ev)) + epsilon
    Sigma_pd <- Sigma + shift * diag(p)
#    cat("Added diagonal shift:", shift, "\n")
  } else {
    Sigma_pd <- Sigma
  }
  return(Sigma_pd)
}

# n <- 300
# p <- 1000
# rho <- 0.5  
# Sigma_pd <- generate_cov_matrix2(p,rho)
# 
# b1 <- 1
# b2 <- 2.3
# b3 <- 0  # Adjust alpha to tune the curvature of x3's effect
# 
# intercept <- -1
# beta1 <- 3
# beta2 <- -4
# beta3 <- -10
# gam <- b1*exp(-1*(b1^2/2))*beta3+beta1 
# beta4 <- -sqrt(rho) * gam * (1/Sigma_pd[1,1])
# beta5 <- 2.5
# beta6 <- -2
# beta7 <- - ( sqrt(rho) * ( beta5 + beta6 * ( 3 * Sigma_pd[6,6] - b2^2 ) ) ) / Sigma_pd[7,7]
# beta_v <- c(beta1, beta2, beta3, beta4, beta5, beta6, beta7,rep(0, p - 7))


generate_data3 <- function(n, p, rho,Sigma_pd) {
  X <- mvrnorm(n, mu = rep(0, p), Sigma = Sigma_pd)
  colnames(X) <- paste0("x", 1:p)
  
  # Define the intercept
  intercept <- intercept
  b1 <- b1
  b2 <- b2
  
  # Construct the linear predictor
  eta <- intercept + 
    beta_v[1] * X[,1] + 
    beta_v[2] * (X[,2]^2) + 
    beta_v[3] * (sin(b1 * X[,3])) +   # Modified non-linear term for x3
    beta_v[4] * X[,4] +
    beta_v[5] * X[,5] +
    beta_v[6] * sin(b2 * X[,6]) +   # Modified non-linear term for x3
    #beta_v[6] * (X[,6]*(X[,6]-b2)*(X[,6]+b2)) +
    beta_v[7] * X[,7] 
  
  p_i <- exp(eta) / (1 + exp(eta))
  
#  hist(p_i)  # Check the distribution of pi
  
  Y <- rbinom(n = n, size = 1, prob = p_i)
  Y <- as.factor(Y)
  
  return(list(X = X, Y = Y))
}

# data <- generate_data3(n = n,p = p,rho = rho,Sigma_pd = Sigma_pd)
# summary(data$Y)





############## setting 4
# set.seed(1)
# n <- 300
# p <- 1000
# rho <- 0.5  
# b1 <- 1.7
# intercept <- 0
# beta1 <- 1
# beta2 <- 1
# beta3 <- 1
# beta4 <- -b1 * sqrt(rho) * (beta1 + beta2) * exp(-1*(beta1^2 + beta2^2 + 2*rho*beta1*beta2)*1/2) - beta3*sqrt(rho)
# beta_v <- c(beta1, beta2, beta3, beta4, rep(0, p - 4))



generate_data4 <- function(n, p, rho) {
  Sigma <- generate_cov_matrix(p, rho)
  X <- mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  colnames(X) <- paste0("x", 1:p)
  
  # Define the intercept
  intercept <- intercept
  b1 <- b1
  
  # Construct the linear predictor
  eta <- intercept + b1*sin( 
    beta_v[1] * X[,1] + 
    beta_v[2] * X[,2]) + 
    beta_v[3] * X[,3] + 
    beta_v[4] * X[,4]
  
  p_i <- exp(eta) / (1 + exp(eta))
  
  #  hist(p_i)  # Check the distribution of pi
  
  Y <- rbinom(n = n, size = 1, prob = p_i)
  Y <- as.factor(Y)
  
  return(list(X = X, Y = Y))
}

# data <- generate_data4(n = n,p = p,rho = rho)
# summary(data$Y)


######### report function
simreport <- function(results, rfe_features, isis_features, sparse_features, raw_features) {
  rpn <- length(results)                       
  
  ## helper to map "xN" → "N" so RFE & ISIS can share a column ----------------
  canonise <- function(v) sub("^x", "", v)
  
  ## unified list of feature labels actually used by either method -----------
  canon_feats <- sort(unique(c(canonise(rfe_features),
                               canonise(isis_features),
                               canonise(raw_features),
                               canonise(sparse_features))))
  
  ## ------- core metric calculator ------------------------------------------
  metric_block <- function(sel_name, truth_vec) {
    truth_canon <- canonise(truth_vec)
    
    ## hit matrix: rows = truth_canon, cols = replicates ---------------------
    hit_mat <- vapply(results, function(res) {
      sel <- canonise(res[[sel_name]])
      truth_canon %in% sel
    }, logical(length(truth_canon)))
    rownames(hit_mat) <- truth_canon
    
    ## per-replicate aggregates ---------------------------------------------
    tp_counts <- colSums(hit_mat)
    fp_counts <- sapply(results, function(res)
      length(setdiff(canonise(res[[sel_name]]), truth_canon)))
    fp_rates  <- sapply(results, function(res) {
      sel <- canonise(res[[sel_name]])
      if (length(sel) == 0) 0 else length(setdiff(sel, truth_canon)) / length(sel)
    })
    tot_sel   <- sapply(results, function(res) length(res[[sel_name]]))
    
    ## per-feature TPR -------------------------------------------------------
    tpr_vec <- rowMeans(hit_mat)
    names(tpr_vec) <- truth_canon
    
    ## fill feature columns; initialise only the *union* labels -------------
    feat_cols <- setNames(rep(NA_real_, length(canon_feats)), canon_feats)
    feat_cols[names(tpr_vec)] <- tpr_vec         # NA where feature not “true”
    
    c(feat_cols,
      tp_count      = mean(tp_counts),
      tp_count_var  = std.error(tp_counts),
      fp_count      = mean(fp_counts),
      fp_count_var  = std.error(fp_counts),
      fp_rate       = mean(fp_rates),
      fp_rate_var   = std.error(fp_rates),
      tot_sel       = mean(tot_sel),
      tot_sel_var   = std.error(tot_sel))
  }
  ## -------------------------------------------------------------------------
  
  rfe_row    <- metric_block("svs_rfe_select", rfe_features)
  isis_row   <- metric_block("isis_select",    isis_features)
  sparse_row <- metric_block("sparse_select",  sparse_features)
  raw_row <- metric_block("raw_rfe_select",  raw_features)
  
  out <- as.data.frame(rbind(rfe = rfe_row, isis = isis_row, raw = raw_row, sparse = sparse_row),
                       check.names = FALSE)   # keep column labels verbatim
  return(out)
}
