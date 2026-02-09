# ============= UTILITY FUNCTIONS # =============
# Helper functions for diagnostics, F-statistics, and visualization

# Suppress CMD check notes for ggplot2 NSE
if(getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "Beta", "CI_Lower_Trad", "CI_Upper_Trad", "X_axis",
    "effect", "effect1_low", "effect1_up", "effect2_low", "effect2_up",
    "effect_low", "effect_up", "time",
    "true_shape", "true_shape1", "true_shape2"
  ))
}

#' Calculate F-statistics and Q-statistic for instrument strength (internal)
#' 
#' @param J Number of genetic instruments
#' @param K Number of exposures
#' @param PC Vector of indices indicating which columns in datafull correspond to the principal components
#' @param datafull Data frame containing instruments (first J columns) and principal components (subsequent columns) [G, X]
#' @param Y Optional outcome vector; if provided, Q-statistic for overidentification is calculated)
#' @return Matrix with columns: PC (component index), RR (R-squared), FF (F-statistic), cFF (conditional F-statistic). If Y is provided, additional columns: Qvalue (Hansen's J overidentification test statistic), df (degrees of freedom for Q-test), pvalue (p-value for Q-test from chi-squared distribution).
#' @export
IS <- function(J, K, PC, datafull, Y = NULL) {
  X <- data.matrix(datafull)
  FF <- rep(NA, K)
  RR <- rep(NA, K)
  cFF <- rep(NA, K)
  
  for (k in 1:K) {
    t <- PC[k]
    fitX <- stats::lm(X[, J + t] ~ matrix(as.numeric(X[, 1:J]), dim(X)[1], J))
    FF[k] <- as.numeric(summary(fitX)$fstatistic[1])
    RR[k] <- as.numeric(summary(fitX)$r.squared[1])
  }
  
  for (k in 1:K) {
    XXm <- c()
    for (kk in (1:K)[-k]) {
      t <- PC[kk]
      fitX <- stats::lm(X[, J + t] ~ matrix(as.numeric(X[, 1:J]), dim(X)[1], J))
      Xhat <- X[, J + t] - summary(fitX)$residuals
      XXm <- cbind(XXm, as.numeric(Xhat))
    }
    t <- PC[k]
    fitX <- stats::lm(X[, J + t] ~ XXm)
    resi <- as.numeric(summary(fitX)$residuals)
    fit <- stats::lm(resi ~ matrix(as.numeric(X[, 1:J]), dim(X)[1], J))
    cFF[k] <- as.numeric(summary(fit)$fstatistic[1]) * summary(fit)$fstatistic[2] / 
      (summary(fit)$fstatistic[2] - (K - 1))
  }
  
  # Calculate Q-statistic (Hansen's J-test for overidentification)
  # Only if outcome Y is provided
  if (!is.null(Y)) {
    n <- nrow(X)
    Z <- X[, 1:J]
    XX <- X[, (J + PC)]
    
    # Center variables
    Y_centered <- Y - mean(Y)
    Z_centered <- scale(Z, center = TRUE, scale = FALSE)
    XX_centered <- scale(XX, center = TRUE, scale = FALSE)
    
    # Two-stage least squares
    # First stage: predict X from Z
    XX_hat <- Z_centered %*% solve(t(Z_centered) %*% Z_centered) %*% t(Z_centered) %*% XX_centered
    
    # Second stage: Y ~ X_hat
    beta_2sls <- solve(t(XX_hat) %*% XX_hat) %*% t(XX_hat) %*% Y_centered
    
    # Residuals
    residuals <- Y_centered - XX_centered %*% beta_2sls
    
    # Hansen's J statistic (overidentification test)
    # Q = n * g'(theta) * W * g(theta) where g = Z'e/n
    g_theta <- t(Z_centered) %*% residuals / n
    W <- solve(t(Z_centered) %*% Z_centered / n)
    Q_stat <- as.numeric(n * t(g_theta) %*% W %*% g_theta)
    
    # Degrees of freedom: J - K (instruments minus endogenous variables)
    Q_df <- J - K
    
    # P-value from chi-squared distribution
    Q_pval <- stats::pchisq(Q_stat, df = Q_df, lower.tail = FALSE)
    
    # Return with Q-statistics
    return(cbind(PC, RR, FF, cFF, 
                 Qvalue = rep(Q_stat, K), 
                 df = rep(Q_df, K), 
                 pvalue = rep(Q_pval, K)))
  }
  
  return(cbind(PC, RR, FF, cFF))
}


#' Get true shape values for simulation
#' 
#' @param workGrid Grid of time points for evaluation
#' @param XYmodel Model code ('0'-'9') specifying the true effect shape
#' @return Vector of true effect values at workGrid time points
#' @keywords internal
get_true_shape_values <- function(workGrid, XYmodel) {
  model_functions <- list(
    '0' = function(t) { 0 * (t < Inf) },
    '1' = function(t) { 0.1 * (t < Inf) },
    '2' = function(t) { 0.02 * t },
    '3' = function(t) { 0.5 - 0.02 * t },
    '4' = function(t) { 0.1 * (t < 20) },
    '5' = function(t) { 0.1 * (t > 30) },
    '6' = function(t) { 0.05 * (-t + 20) * (t < 20) },
    '7' = function(t) { 0.05 * (t - 30) * (t > 30) },
    '8' = function(t) { 0.002 * t^2 - 0.11 * t + 0.5 },
    '9' = function(t) { -0.00002 * t^3 + 0.004 * t^2 - 0.2 * t + 1 }
  )
  
  if (!XYmodel %in% names(model_functions)) {
    stop("Invalid XYmodel provided")
  }
  
  fun <- model_functions[[XYmodel]]
  return(fun(workGrid))
}

#' Get true effect function for simulation
#' 
#' @param model_code Model code ('0'-'9') specifying the effect shape function
#' @return Function that takes time as input and returns effect value
#' @keywords internal
get_true_effect_function <- function(model_code) {
  effect_functions <- list(
    "0" = function(t) 0 * (t < Inf),
    "1" = function(t) 0.1 * (t < Inf),
    "2" = function(t) 0.02 * t,
    "3" = function(t) 0.5 - 0.02 * t,
    "4" = function(t) 0.1 * (t < 20),
    "5" = function(t) 0.1 * (t > 30),
    "6" = function(t) 0.05 * (-t + 20) * (t < 20),
    "7" = function(t) 0.05 * (t - 30) * (t > 30),
    "8" = function(t) 0.002 * t^2 - 0.11 * t + 0.5,
    "9" = function(t) -0.00002 * t^3 + 0.004 * t^2 - 0.2 * t + 1
  )
  
  model_char <- as.character(model_code)
  if (!model_char %in% names(effect_functions)) {
    stop("Unknown model code: ", model_char)
  }
  
  return(effect_functions[[model_char]])
}
