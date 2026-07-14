# ============= GMM ESTIMATION FUNCTIONS # =============
# Generalized Method of Moments estimation for continuous and binary outcomes

#' GMM estimation for continuous outcome 
#'
#' @param X Matrix of exposure principal components (N x K)
#' @param Y Outcome vector (length N)
#' @param Z Genetic instrument matrix (N x J)
#' @param beta0 Initial values for beta (default NA, uses zero initialization)
#' @return List with gmm_est, gmm_se, variance_matrix, gmm_pval, Q_stat, Q_pval
#' @examples
#' set.seed(1)
#' n <- 200; J <- 5; K <- 2
#' Z <- matrix(rbinom(n * J, 2, 0.3), n, J)
#' X <- Z[, 1:K] + matrix(rnorm(n * K, sd = 0.5), n, K)
#' Y <- as.numeric(X %*% c(1, -0.5) + rnorm(n))
#' fit <- gmm_lm_onesample(X, Y, Z)
#' fit$gmm_est
#' @export
gmm_lm_onesample <- function(X, Y, Z, beta0 = NA) {
  n <- nrow(Z)
  J <- ncol(Z)
  K <- ncol(X)
  
  # Mean-center variables  
  Y <- Y - mean(Y)
  X0 <- function(k) {X[, k] - mean(X[, k])}
  X <- sapply(1:K, X0)
  Z0 <- function(k) {Z[, k] - mean(Z[, k])}
  Z <- sapply(1:J, Z0)
  
  # GMM functions
  g <- function(bet) {
    as.vector(t(Z) %*% (Y - (X %*% bet))) / n
  }
  
  Om <- function(bet) {
    (t(Z) %*% Z / n) * as.numeric(t(Y - (X %*% bet)) %*% (Y - (X %*% bet)) / n)
  }
  
  G <- -(t(Z) %*% X) / n
  
  Q <- function(bet) {
    as.numeric(t(g(bet)) %*% solve(Om(bet)) %*% g(bet))
  }
  
  # Preliminary estimate
  Q.gg <- function(bet) {
    as.numeric(t(g(bet)) %*% g(bet))
  }
  bet.gg <- nlminb(rep(0, K), objective = Q.gg)$par
  
  # GMM estimate
  DQ <- function(bet) {
    2 * as.matrix(t(G) %*% solve(Om(bet)) %*% g(bet))
  }
  
  gmm <- nlminb(bet.gg, objective = Q, gradient = DQ)$par
  
  var.gmm <- as.matrix(solve(t(G) %*% solve(Om(gmm)) %*% G))
  
  # Q statistics
  Qstat <- Q(gmm)
  Q.pval <- stats::pchisq(n * Qstat, df = J - K, lower.tail = FALSE)
  
  gmm.pval <- 2 * stats::pnorm(-abs(gmm / sqrt(diag(var.gmm) / n)))
  
  res.list <- list(
    "gmm_est" = gmm,
    "gmm_se" = sqrt(diag(var.gmm) / n),
    "variance_matrix" = var.gmm / n,
    "gmm_pval" = 2 * stats::pnorm(-abs(gmm / sqrt(diag(var.gmm) / n))),
    "Q_stat" = Qstat,
    "Q_pval" = Q.pval
  )
  
  return(res.list)
}

# ============= CONTROL FUNCTION ESTIMATION # =============
# 2SRI-IV estimation (with LASSO implementation)

#' Control function for logit model
#' 
#' @param X Matrix of exposure principal components (N x K)
#' @param Y Binary outcome vector (0/1, length N)
#' @param Z Genetic instrument matrix (N x J)
#' @param alpha Elastic net mixing parameter (1=lasso, 0=ridge)
#' @param nfolds Number of cross-validation folds for lambda selection
#' @param standardize Standardize variables before fitting
#' @param use_lasso Use LASSO regularization in first stage. If FALSE, uses OLS.
#' @return List with gmm_est, gmm_se, variance_matrix, gmm_pval
#' @examples
#' set.seed(1)
#' n <- 200; J <- 5; K <- 2
#' Z <- matrix(rbinom(n * J, 2, 0.3), n, J)
#' X <- Z[, 1:K] + matrix(rnorm(n * K, sd = 0.5), n, K)
#' lin_pred <- X %*% c(0.8, -0.4)
#' Y <- rbinom(n, 1, plogis(lin_pred))
#' fit <- cf_logit(X, Y, Z)
#' fit$gmm_est
#' @export
cf_logit <- function(X, Y, Z, alpha = 1, nfolds = 10, standardize = TRUE, use_lasso = FALSE) {
  n <- nrow(Z)
  J <- ncol(Z)
  K <- ncol(X)
  
  # Center variables
  X_centered <- scale(X, center = TRUE, scale = FALSE)
  Z_centered <- scale(Z, center = TRUE, scale = FALSE)
  
  # Initialize matrix for first-stage residuals
  X_hat_res <- matrix(NA, nrow = n, ncol = K)
  models <- list()
  
  # First stage: Predict each X using regularized regression or linear regression on Z
  for (i in 1:K) {
    if (use_lasso) {
      cv_model <- suppressWarnings(
        suppressMessages(
          cv.glmnet(
            Z_centered, X_centered[, i],
            alpha = alpha,
            nfolds = nfolds,
            standardize = standardize,
            intercept = TRUE
          )
        )
      )
      
      
      models[[i]] <- cv_model
      X_hat <- predict(cv_model, newx = Z_centered, s = "lambda.min")
      X_hat_res[, i] <- X_centered[, i] - X_hat
    } else {
      lm_model <- stats::lm(X_centered[, i] ~ Z_centered)
      models[[i]] <- lm_model
      X_hat_res[, i] <- stats::residuals(lm_model)
    }
  }
  
  # Second stage: Include original X and residuals
  cf_data <- data.frame(Y = Y, X_centered, X_hat_res)
  colnames(cf_data) <- c("Y", paste0("PC", 1:K), paste0("Resid", 1:K))
  
  # Fit logistic regression
  cf_model <- stats::glm(Y ~ ., data = cf_data, family = stats::binomial(link = "logit"))
  cf_summary <- summary(cf_model)
  
  res.list <- list(
    "gmm_est" = cf_summary$coefficients[2:(K + 1), 1],
    "gmm_se" = cf_summary$coefficients[2:(K + 1), 2],
    "variance_matrix" = stats::vcov(cf_model)[2:(K + 1), 2:(K + 1)],
    "gmm_pval" = cf_summary$coefficients[2:(K + 1), 4]
  )
  
  return(res.list)
}


# ============= GMM ESTIMATION - TWO SAMPLE # =============
#' Two-sample GMM
#' 
#' @param bx Matrix J x K of first-stage coefficients (SNP -> PC associations)
#' @param by Vector length J of outcome GWAS betas
#' @param sy Vector length J of outcome GWAS standard errors  
#' @param ny Outcome GWAS sample size
#' @return  List with gmm_est, gmm_se, variance_matrix, gmm_pval, Q_stat, Q_df, Q_pval
#' @examples
#' set.seed(1)
#' J <- 10; K <- 2
#' bx <- matrix(rnorm(J * K, sd = 0.3), J, K)
#' by <- bx %*% c(0.5, -0.2) + rnorm(J, sd = 0.05)
#' sy <- runif(J, 0.02, 0.05)
#' fit <- gmm_twosample_simple(bx, by, sy, ny = 50000)
#' fit$gmm_est
#' @export
gmm_twosample_simple <- function(bx, by, sy, ny) {
  
  J <- nrow(bx)  # Number of instruments
  K <- ncol(bx)  # Number of PCs/exposures
  
  # Inverse variance weighting
  W <- diag(1 / sy^2)
  
  # Two-sample GMM: beta = (bx' W bx)^{-1} bx' W by
  beta_hat <- solve(t(bx) %*% W %*% bx) %*% t(bx) %*% W %*% by
  
  # Variance
  Sigma_y <- diag(sy^2)
  var_beta <- solve(t(bx) %*% W %*% bx) %*%  t(bx) %*% W %*% Sigma_y %*% W %*% bx %*% solve(t(bx) %*% W %*% bx)
  
  # Standard errors and p-values
  se_beta <- sqrt(diag(var_beta))
  p_values <- 2 * stats::pnorm(-abs(beta_hat / se_beta))
  
  # Q-statistic
  residuals <- by - bx %*% beta_hat
  Q_stat <- as.numeric(t(residuals) %*% W %*% residuals)
  Q_df <- J - K
  Q_pval <- stats::pchisq(Q_stat, df = Q_df, lower.tail = FALSE)
  
  list(
    gmm_est = as.vector(beta_hat),
    gmm_se = se_beta,
    variance_matrix = var_beta,
    gmm_pval = p_values,
    Q_stat = Q_stat,
    Q_df = Q_df,
    Q_pval = Q_pval
  )
}
