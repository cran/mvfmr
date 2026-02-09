# ============= AUTOMATIC MULTIVARIABLE MVFMR - JOINT ESTIMATION # =============
# Joint estimation of multiple time-varying exposures with automatic 
# component selection via cross-validation

#' Automatic Multivariable Functional MR with joint estimation (internal)
#'
#' Core function that performs joint estimation of time-varying causal effects
#' from multiple correlated exposures using automatic component selection.
#'
#' @param Gmatrix Genetic instrument matrix (N × J)
#' @param res1 FPCA result for exposure 1 
#' @param res2 FPCA result for exposure 2 
#' @param Yvector Outcome vector
#' @param IDmatch Optional index vector to match rows of Gmatrix and Yvector (default: 1:N)
#' @param nPC1_selected Fixed number of principal components to retain for exposure 1 (NA = select automatically)
#' @param max_nPC1 Maximum number of principal components to consider for exposure 1 during selection
#' @param nPC2_selected Fixed number of principal components to retain for exposure 2 (NA = select automatically)
#' @param max_nPC2 Maximum number of principal components to consider for exposure 2 during selection
#' @param X1_true Optional: true X1 curves (simulation only)
#' @param X2_true Optional: true X2 curves (simulation only)
#' @param method Estimation method: "gmm" (Generalized Method of Moments), "cf" (control function), or "cf-lasso" (control function with Lasso)
#' @param basis Basis type for functional representation: "eigenfunction" or "polynomial"
#' @param outcome Outcome type: "continuous" for numeric or "binary" for 0/1 outcomes
#' @param bootstrap Logical; whether to perform bootstrap inference for confidence intervals
#' @param n_B Number of bootstrap iterations (used only if bootstrap = TRUE)
#' @param improvement_threshold Minimum cross-validation improvement required to add an additional principal component
#' @param X1Ymodel Optional: true effect model for X1 on Y (simulation only)
#' @param X2Ymodel Optional: true effect model for X2 on Y (simulation only)
#' @param num_cores_set Number of CPU cores to use for parallel processing
#' @param verbose Print progress messages and diagnostics during computation
#'
#' @return List with estimation results, selected components, performance metrics
#' @keywords internal
AUTOMATIC_Multi_MVFMR <- function(Gmatrix, 
                                  res1 = NA, 
                                  res2 = NA, 
                                  Yvector, 
                                  IDmatch = NA,
                                  nPC1_selected = NA, 
                                  max_nPC1 = NA,
                                  nPC2_selected = NA,
                                  max_nPC2 = NA,
                                  X1_true = matrix(),
                                  X2_true = matrix(),
                                  method = "gmm",
                                  basis = "eigenfunction",
                                  outcome = "continuous",
                                  bootstrap = FALSE,
                                  n_B = 10,
                                  improvement_threshold = 0.01,
                                  X1Ymodel = NA, 
                                  X2Ymodel = NA,
                                  num_cores_set = NA,
                                  verbose = FALSE) {
  
  
  # 1. VALIDATION
  
  
  if (!basis %in% c("eigenfunction", "polynomial")) {
    stop("The 'basis' parameter must be either 'eigenfunction' or 'polynomial'")
  }
  
  if (!method %in% c("gmm", "cf", "cf-lasso")) {
    stop("The 'method' parameter must be one of 'gmm', 'cf', or 'cf-lasso'")
  }
  
  if (!outcome %in% c("continuous", "binary")) {
    stop("The 'outcome' parameter must be either 'continuous' or 'binary'")
  }
  
  
  # 2. INITIALIZATION
  
  
  fitRES <- list()
  Gmatrix <- as.data.frame(Gmatrix)
  
  max_npc1 <- ifelse(is.na(max_nPC1) | max_nPC1 > ncol(res1$xiEst), 
                     ncol(res1$xiEst), max_nPC1)
  max_npc2 <- ifelse(is.na(max_nPC2) | max_nPC2 > ncol(res2$xiEst), 
                     ncol(res2$xiEst), max_nPC2)
  
  J <- ncol(Gmatrix)
  
  workGrid_diff1 <- res1$workGrid[3] - res1$workGrid[2]
  workGrid_diff2 <- res2$workGrid[3] - res2$workGrid[2]
  
  # Prepare fitted curves
  if (any(is.na(X1_true))) {
    X1_curves <- fitted(res1, ciOptns = list(alpha = 0.05, kernelType = 'gauss'))$fitted
  } else {
    X1_curves <- t(apply(X1_true, 1, function(x) {
      spline(x = seq(1, 50, length.out = ncol(X1_true)), 
             y = x, xout = res1$workGrid)$y
    }))
  }
  
  if (any(is.na(X2_true))) {
    X2_curves <- fitted(res2, ciOptns = list(alpha = 0.05, kernelType = 'gauss'))$fitted
  } else {
    X2_curves <- t(apply(X2_true, 1, function(x) {
      spline(x = seq(1, 50, length.out = ncol(X2_true)), 
             y = x, xout = res2$workGrid)$y
    }))
  }
  
  
  # 3. HELPER FUNCTIONS
  
  
  get_polynomial_basis <- function(grid, n_components) {
    outer(grid, 0:(n_components - 1), `^`)
  }
  
  run_model <- function(X, Y, Z, method_type, outcome_type) {
    if (method_type == "gmm" && outcome_type == "continuous") {
      return(gmm_lm_onesample(X = X, Y = Y, Z = Z))
    } else if (method_type == "cf" && outcome_type == "binary") {
      return(cf_logit(X = X, Y = Y, Z = Z, use_lasso = FALSE))
    } else if (method_type == "cf-lasso" && outcome_type == "binary") {
      return(cf_logit(X = X, Y = Y, Z = Z, use_lasso = TRUE))
    }
  }
  
  
  # 4. PARALLEL SETUP
  
  
  if (is.na(num_cores_set)) {
    num_cores <- parallel::detectCores() - 1 
  } else {
    num_cores <- num_cores_set
  } 
  
  cl <- parallel::makeCluster(num_cores)
  doParallel::registerDoParallel(cl)
  
  if(verbose == TRUE) print("***Feature selection start***")
  
  
  # 5. CROSS-VALIDATION COMPONENT SELECTION
  
  
  if (basis == "polynomial") {
    bbb1 <- get_polynomial_basis(res1$workGrid, max_npc1)
    bbb2 <- get_polynomial_basis(res2$workGrid, max_npc2)
    phi_transposed1 <- t(res1$phi[-1, 1:max_npc1])
    phi_transposed2 <- t(res2$phi[-1, 1:max_npc2])
  }
  
  n_folds <- 3
  if(verbose == TRUE) print("***Feature selection starting***")
  
  valid_indices <- which(!is.na(IDmatch))
  n_valid <- length(valid_indices)
  
  set.seed(123)
  fold_indices <- sample(rep(1:n_folds, length.out = n_valid))
  
  # Helper: evaluate component combination
  evaluate_components <- function(nPC1, nPC2) {
    fold_metrics <- numeric(n_folds)
    
    for (fold in 1:n_folds) {
      val_idx <- valid_indices[fold_indices == fold]
      train_idx <- valid_indices[fold_indices != fold]
      
      PC_1 <- na.omit(as.matrix(res1$xiEst[, 1:nPC1]))  
      PC_2 <- na.omit(as.matrix(res2$xiEst[, 1:nPC2]))  
      PC_ <- cbind(PC_1, PC_2)
      
      Z_train <- Gmatrix[train_idx, , drop = FALSE]
      X_train <- PC_[train_idx, , drop = FALSE]
      Y_train <- Yvector[IDmatch[train_idx]]
      
      Z_val <- Gmatrix[val_idx, , drop = FALSE]
      X_val <- PC_[val_idx, , drop = FALSE]
      Y_val <- Yvector[IDmatch[val_idx]]
      
      if (basis == "polynomial") {
        B1 <- phi_transposed1[1:nPC1, ] %*% bbb1[-1, 1:nPC1] * workGrid_diff1  
        B2 <- phi_transposed2[1:nPC2, ] %*% bbb2[-1, 1:nPC2] * workGrid_diff2
        
        X1_train_model <- PC_1[train_idx, ] %*% B1
        X2_train_model <- PC_2[train_idx, ] %*% B2
        X_train_model <- cbind(X1_train_model, X2_train_model)
        
        X1_val_model <- PC_1[val_idx, ] %*% B1
        X2_val_model <- PC_2[val_idx, ] %*% B2
        X_val_model <- cbind(X1_val_model, X2_val_model)
      } else {
        X_train_model <- X_train
        X_val_model <- X_val
      }
      
      tryCatch({
        suppressWarnings({
          gmm_res <- run_model(X = X_train_model, Y = Y_train, Z = Z_train, 
                               method_type = method, outcome_type = outcome)
        })
        
        if (basis == "polynomial") {
          X1_curves_val <- X1_curves[val_idx, ]
          X2_curves_val <- X2_curves[val_idx, ]
          
          val_predictions <- X1_curves_val %*% bbb1[, 1:nPC1] %*% gmm_res$gmm_est[1:nPC1] +
            X2_curves_val %*% bbb2[, 1:nPC2] %*% gmm_res$gmm_est[(nPC1 + 1):(nPC1 + nPC2)]
        } else {
          X1_curves_val <- X1_curves[val_idx, ]
          X2_curves_val <- X2_curves[val_idx, ]
          
          val_predictions <- X1_curves_val %*% (res1$phi)[, 1:nPC1] %*% gmm_res$gmm_est[1:nPC1] +
            X2_curves_val %*% (res2$phi)[, 1:nPC2] %*% gmm_res$gmm_est[(nPC1 + 1):(nPC1 + nPC2)]
        }
        
        if (outcome == "continuous") {
          fold_metrics[fold] <- mean((Y_val - val_predictions)^2)
        } else if (outcome == "binary") {
          val_predictions_prob <- 1 / (1 + exp(-val_predictions))
          fold_metrics[fold] <- pROC::auc(pROC::roc(response  = Y_val, predictor = as.numeric(val_predictions_prob), levels    = c(0, 1), direction = "<", quiet     = TRUE))
          
        }
      }, error = function(e) {
        if (outcome == "continuous") {
          fold_metrics[fold] <- Inf
        } else {
          fold_metrics[fold] <- 0
        }
        warning(paste("Model fitting failed for nPC1 =", nPC1, "and nPC2 =", nPC2, "in fold", fold))
      })
    }
    
    mean(fold_metrics)
  }
  
  is_better <- function(new_metric, current_best) {
    if (outcome == "continuous") {
      return(new_metric < current_best)
    } else {
      return(new_metric > current_best)
    }
  }
  
  # Sequential component selection
  component_results <- data.frame(
    nPC1 = integer(),
    nPC2 = integer(),
    metric = numeric(),
    step = integer()
  )
  
  best_nPC1 <- 2
  best_nPC2 <- 2
  
  best_metric <- suppressWarnings({evaluate_components(best_nPC1, best_nPC2)})
  
  component_results <- rbind(component_results, data.frame(
    nPC1 = best_nPC1, 
    nPC2 = best_nPC2,
    metric = best_metric,
    step = 0
  ))
  
  max_components <- min(max_npc1 + max_npc2, 12)
  
  improved <- TRUE
  step <- 1
  
  while (improved && (best_nPC1 < max_npc1 || best_nPC2 < max_npc2) && 
         (best_nPC1 + best_nPC2) < max_components) {
    improved <- FALSE
    
    if (best_nPC1 < max_npc1) {
      if(verbose == TRUE) print(paste("Trying nPC1 =", best_nPC1 + 1, "and nPC2 =", best_nPC2))
      metric1 <- suppressWarnings({evaluate_components(best_nPC1 + 1, best_nPC2)})
      
      component_results <- rbind(component_results, data.frame(
        nPC1 = best_nPC1 + 1, 
        nPC2 = best_nPC2,
        metric = metric1,
        step = step
      ))
    } else {
      metric1 <- if (outcome == "continuous") Inf else 0
    }
    
    if (best_nPC2 < max_npc2) {
      if(verbose == TRUE) print(paste("Trying nPC1 =", best_nPC1, "and nPC2 =", best_nPC2 + 1))
      metric2 <- suppressWarnings({evaluate_components(best_nPC1, best_nPC2 + 1)})
      
      component_results <- rbind(component_results, data.frame(
        nPC1 = best_nPC1, 
        nPC2 = best_nPC2 + 1,
        metric = metric2,
        step = step
      ))
    } else {
      metric2 <- if (outcome == "continuous") Inf else 0
    }
    
    if (outcome == "continuous") {
      improvement1 <- (best_metric - metric1) / best_metric
      improvement2 <- (best_metric - metric2) / best_metric
    } else {
      improvement1 <- (metric1 - best_metric) / best_metric
      improvement2 <- (metric2 - best_metric) / best_metric
    }
    
    if (improvement1 > improvement_threshold || improvement2 > improvement_threshold) {
      improved <- TRUE
      
      if (is_better(metric1, metric2)) {
        best_nPC1 <- best_nPC1 + 1
        best_metric <- metric1
        if(verbose == TRUE) print(paste("Improvement found: increasing nPC1 to", best_nPC1, "with metric", best_metric))
      } else {
        best_nPC2 <- best_nPC2 + 1
        best_metric <- metric2
        if(verbose == TRUE) print(paste("Improvement found: increasing nPC2 to", best_nPC2, "with metric", best_metric))
      }
    } else {
      if(verbose == TRUE) print("No significant improvement found. Stopping component selection.")
    }
    
    step <- step + 1
  }
  
  if(verbose == TRUE) print(paste("Final selected components: nPC1 =", best_nPC1, "and nPC2 =", best_nPC2))
  if(verbose == TRUE) print(paste("Final metric:", best_metric))
  
  fitRES$component_selection_results <- component_results
  
  
  # 6. FINAL MODEL WITH SELECTED COMPONENTS
  
  if(!is.na(nPC1_selected) & !is.na(nPC2_selected))
  {
    best_nPC1 = nPC1_selected
    best_nPC2 = nPC2_selected
  }
  PC_1 <- na.omit(as.matrix(res1$xiEst[, 1:best_nPC1]))  
  PC_2 <- na.omit(as.matrix(res2$xiEst[, 1:best_nPC2]))  
  PC_ <- cbind(PC_1, PC_2)
  
  Z_GMMused <- Gmatrix[!is.na(IDmatch), , drop = FALSE]
  X_GMMused <- PC_[!is.na(IDmatch), , drop = FALSE]
  Y_GMMused <- Yvector[IDmatch][!is.na(IDmatch)]
  
  if (basis == "polynomial") {
    B1 <- phi_transposed1[1:best_nPC1, ] %*% bbb1[-1, 1:best_nPC1] * workGrid_diff1
    B2 <- phi_transposed2[1:best_nPC2, ] %*% bbb2[-1, 1:best_nPC2] * workGrid_diff2
    X1_model <- PC_1 %*% B1
    X2_model <- PC_2 %*% B2
    X_model <- cbind(X1_model, X2_model)
  } else {
    X_model <- X_GMMused
  }
  
  suppressWarnings({
    gmm_res <- run_model(X = X_model, Y = Y_GMMused, Z = Z_GMMused,
                         method_type = method, outcome_type = outcome)
  })
  
  
  if (basis == "polynomial") {
    y_estimated <- t(X1_curves %*% bbb1[, 1:best_nPC1] %*% gmm_res$gmm_est[1:best_nPC1] +
                       X2_curves %*% bbb2[, 1:best_nPC2] %*% gmm_res$gmm_est[(best_nPC1 + 1):(best_nPC1 + best_nPC2)])
  } else {
    y_estimated <- t(X1_curves %*% (res1$phi)[, 1:best_nPC1] %*% gmm_res$gmm_est[1:best_nPC1] +
                       X2_curves %*% (res2$phi)[, 1:best_nPC2] %*% gmm_res$gmm_est[(best_nPC1 + 1):(best_nPC1 + best_nPC2)])
  }
  
  if (outcome == "continuous") {
    gmm_res$mse <- mean((Yvector - y_estimated)^2)
    fitRES$final_mse <- gmm_res$mse
  } else if (outcome == "binary") {
    y_estimated_prob <- 1 / (1 + exp(-t(y_estimated)))
    y_binary_pred <- ifelse(y_estimated_prob > 0.5, 1, 0)
    gmm_res$auc <- pROC::auc(pROC::roc(response  = Yvector, predictor = as.numeric(y_estimated_prob), levels    = c(0, 1), direction = "<", quiet     = TRUE))
    
    fitRES$final_auc <- gmm_res$auc
  }
  
  parallel::stopCluster(cl)
  
  
  # 7. STORE RESULTS
  
  
  nPC1 <- best_nPC1
  nPC2 <- best_nPC2
  
  fitRES$nPC1_used <- nPC1
  fitRES$nPC2_used <- nPC2
  
  # Add true shapes if available
  if (!is.na(X1Ymodel) | !is.na(X2Ymodel)) {
    true_shape1 <- get_true_shape_values(res1$workGrid, as.character(X1Ymodel))
    true_shape2 <- get_true_shape_values(res2$workGrid, as.character(X2Ymodel))
  } else {
    true_shape1 <- NULL
    true_shape2 <- NULL
  }
  
  # Prepare data for plots
  data <- data.frame(
    Beta = gmm_res$gmm_est,
    CI_Lower_Trad = gmm_res$gmm_est - 1.96 * gmm_res$gmm_se,
    CI_Upper_Trad = gmm_res$gmm_est + 1.96 * gmm_res$gmm_se
  )
  data$X_axis <- 1:nrow(data)
  
  final_phi <- cbind((res1$phi)[, 1:nPC1], (res2$phi)[, 1:nPC2])
  
  if (!is.null(true_shape1) & !is.null(true_shape2)) {
    if (basis == "polynomial") {
      true_beta_k <- c(
        t((res1$phi)[, 1:nPC1]) %*% bbb1[, 1:nPC1] %*% 
          solve(t(bbb1[, 1:nPC1]) %*% bbb1[, 1:nPC1]) %*% t(bbb1[, 1:nPC1]) %*% true_shape1,
        t((res2$phi)[, 1:nPC2]) %*% bbb2[, 1:nPC2] %*% 
          solve(t(bbb2[, 1:nPC2]) %*% bbb1[, 1:nPC2]) %*% t(bbb1[, 1:nPC2]) %*% true_shape2
      )
    } else {
      true_beta_k <- c(
        t((res1$phi)[, 1:nPC1]) %*% (c(true_shape1) * workGrid_diff1),
        t((res2$phi)[, 1:nPC2]) %*% (c(true_shape2) * workGrid_diff2)
      )
    }
    
    plot <- ggplot2::ggplot(data, ggplot2::aes(x = factor(X_axis))) +
      ggplot2::geom_point(ggplot2::aes(y = Beta, color = "Estimated beta_k"), size = 2) +
      ggplot2::geom_point(ggplot2::aes(y = true_beta_k, color = "True beta_k"), size = 2) +
      ggplot2::geom_segment(ggplot2::aes(x = X_axis, xend = X_axis, 
                                         y = CI_Lower_Trad, yend = CI_Upper_Trad, 
                                         color = "Traditional CI"), size = 1) +
      ggplot2::geom_text(ggplot2::aes(y = Beta, label = round(Beta, 4)), 
                         color = "darkgreen", hjust = -0.4, vjust = -0.1, size = 3) +  
      ggplot2::theme_minimal() +
      ggplot2::labs(title = "", y = "Estimate (95% CI)") +
      ggplot2::scale_color_manual(values = c("Estimated beta_k" = "darkgreen", 
                                              "True beta_k" = "blue")) +
      ggplot2::scale_x_discrete(labels = c(paste0("Beta1_", 1:nPC1), paste0("Beta2_", 1:nPC2))) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(size = 12), 
                     axis.title.x = ggplot2::element_blank(), 
                     legend.title = ggplot2::element_blank())
    
  } else {
    plot <- ggplot2::ggplot(data, ggplot2::aes(x = factor(X_axis))) +
      ggplot2::geom_point(ggplot2::aes(y = Beta, color = "Estimated beta_k"), size = 2) +
      ggplot2::geom_segment(ggplot2::aes(x = X_axis, xend = X_axis, 
                                         y = CI_Lower_Trad, yend = CI_Upper_Trad, 
                                         color = "Traditional CI"), size = 1) +
      ggplot2::geom_text(ggplot2::aes(y = Beta, label = round(Beta, 4)), 
                         color = "darkgreen", hjust = -0.4, vjust = -0.1, size = 3) +  
      ggplot2::theme_minimal() +
      ggplot2::labs(title = "", y = "Estimate (95% CI)") +
      ggplot2::scale_color_manual(values = c("Estimated beta_k" = "darkgreen", 
                                              "True beta_k" = "blue")) +
      ggplot2::scale_x_discrete(labels = c(paste0("Beta1_", 1:nPC1), paste0("Beta2_", 1:nPC2))) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(size = 12), 
                     axis.title.x = ggplot2::element_blank(), 
                     legend.title = ggplot2::element_blank())
  }
  
  fitRES$plot_beta <- plot
  fitRES$MPCMRest <- gmm_res$gmm_est
  fitRES$MPCMRvar <- gmm_res$variance_matrix
  
  # Define projection matrices
  if (basis == "polynomial") {
    phi1 <- bbb1[, 1:nPC1]
    phi2 <- bbb2[, 1:nPC2]
    phi <- cbind(phi1, phi2)
  } else {
    phi <- cbind((res1$phi)[, 1:nPC1], (res2$phi)[, 1:nPC2])
  }
  
  pointwise_shape_var <- diag(phi %*% gmm_res$variance_matrix %*% t(phi))
  fitRES$pointwise_shape_var <- pointwise_shape_var
  fitRES$pointwise_estimates <- phi %*% gmm_res$gmm_est
  
  # Prepare effect curves
  if (basis == "polynomial") {
    effect1 <- bbb1[, 1:nPC1] %*% gmm_res$gmm_est[1:nPC1]
    effect2 <- bbb2[, 1:nPC2] %*% gmm_res$gmm_est[(nPC1 + 1):(nPC1 + nPC2)]
  } else {
    effect1 <- (res1$phi)[, 1:nPC1] %*% gmm_res$gmm_est[1:nPC1]
    effect2 <- (res2$phi)[, 1:nPC2] %*% gmm_res$gmm_est[(nPC1 + 1):(nPC1 + nPC2)]
  }
  
  ggdata <- data.frame(
    time = res1$workGrid,
    effect1 = effect1,
    effect1_low = effect1 - 1.96 * sqrt(pointwise_shape_var),
    effect1_up = effect1 + 1.96 * sqrt(pointwise_shape_var),
    effect2 = effect2,
    effect2_low = effect2 - 1.96 * sqrt(pointwise_shape_var),
    effect2_up = effect2 + 1.96 * sqrt(pointwise_shape_var),
    true_shape = NaN
  )
  
  if (!is.na(X1Ymodel) | !is.na(X2Ymodel)) {
    ggdata$true_shape1 <- true_shape1
    ggdata$true_shape2 <- true_shape2
  }
  
  fitRES$ggdata <- ggdata
  plotdif <- max(c(ggdata$effect1_up, ggdata$effect2_up)) - 
    min(c(ggdata$effect1_low, ggdata$effect2_low))
  
  # Bootstrap procedure
  if(bootstrap == TRUE) {
    if(verbose == TRUE) print("Starting bootstrap")
    
    # Setup storage for bootstrap results
    causal_effect1_boot <- matrix(NA, nrow = n_B, ncol = nrow(ggdata))
    causal_effect2_boot <- matrix(NA, nrow = n_B, ncol = nrow(ggdata))
    causal_beta1_boot <- matrix(NA, nrow = n_B, ncol = nPC1)
    causal_beta2_boot <- matrix(NA, nrow = n_B, ncol = nPC2)
    
    # Prepare data for bootstrap
    PC_1 <- na.omit(as.matrix(res1$xiEst[, 1:nPC1]))
    PC_2 <- na.omit(as.matrix(res2$xiEst[, 1:nPC2]))
    PC_ <- cbind(PC_1, PC_2)
    
    Z_GMMused <- Gmatrix[!is.na(IDmatch), , drop = FALSE]
    X_GMMused <- PC_[!is.na(IDmatch), , drop = FALSE]
    Y_GMMused <- Yvector[IDmatch][!is.na(IDmatch)]
    
    # Apply basis transformation if needed
    if(basis == "polynomial") {
      # Create transformation matrices
      B1 <- phi_transposed1[1:nPC1, ] %*% bbb1[-1, 1:nPC1] * workGrid_diff1
      B2 <- phi_transposed2[1:nPC2, ] %*% bbb2[-1, 1:nPC2] * workGrid_diff2
      
      # Apply transformations
      X1_model <- PC_1 %*% B1
      X2_model <- PC_2 %*% B2
      X_model <- cbind(X1_model, X2_model)
    } else {
      # Use eigenfunction basis
      X_model <- X_GMMused
    }
    
    # Setup progress bar
    pb <- progress::progress_bar$new(
      format = "  processing [:bar] :percent eta: :eta",
      total = n_B,
      clear = FALSE)
    
    # Bootstrap loop
    for (b in 1:n_B) {
      # Sample with replacement
      bootstrap_ind <- sample(1:nrow(X_model), replace = TRUE)
      
      X.b <- X_model[bootstrap_ind, ]
      Y.b <- Y_GMMused[bootstrap_ind]
      Z.b <- Z_GMMused[bootstrap_ind, ]
      
      # Fit model
      suppressWarnings({
      gmm_res_boot <- run_model(X = X.b, Y = Y.b, Z = Z.b, 
                                method_type = method, outcome_type = outcome)
      })
      
      # Store bootstrap causal effect surfaces
      causal_beta1_boot[b, ] <- gmm_res_boot$gmm_est[1:nPC1]
      causal_beta2_boot[b, ] <- gmm_res_boot$gmm_est[(nPC1+1):(nPC1+nPC2)]
      
      # Transform to get effects based on basis type
      if(basis == "polynomial") {
        causal_effect1_boot[b, ] <- bbb1[, 1:nPC1] %*% gmm_res_boot$gmm_est[1:nPC1]
        causal_effect2_boot[b, ] <- bbb2[, 1:nPC2] %*% gmm_res_boot$gmm_est[(nPC1+1):(nPC1+nPC2)]
      } else {
        causal_effect1_boot[b, ] <- res1$phi[, 1:nPC1] %*% gmm_res_boot$gmm_est[1:nPC1]
        causal_effect2_boot[b, ] <- res2$phi[, 1:nPC2] %*% gmm_res_boot$gmm_est[(nPC1+1):(nPC1+nPC2)]
      }
      
      pb$tick()
    }
    
    alpha = 0.05
    
    # Calculate confidence intervals for beta coefficients
    CI_beta1_k <- data.frame(
      lwr = apply(causal_beta1_boot, 2, function(x) quantile(x, alpha / 2)),
      obs = gmm_res$gmm_est[1:nPC1], 
      upr = apply(causal_beta1_boot, 2, function(x) quantile(x, 1 - alpha / 2))
    )
    
    CI_beta2_k <- data.frame(
      lwr = apply(causal_beta2_boot, 2, function(x) quantile(x, alpha / 2)),
      obs = gmm_res$gmm_est[(nPC1+1):(nPC1+nPC2)], 
      upr = apply(causal_beta2_boot, 2, function(x) quantile(x, 1 - alpha / 2))
    )
    
    # Apply exponentiation for binary outcomes
    if(outcome == "binary") {
      CI_beta1_k <- exp(CI_beta1_k)
      CI_beta2_k <- exp(CI_beta2_k)
    }
    
    # Calculate confidence intervals for effects
    CI_beta1_t <- data.frame(
      lwr = apply(causal_effect1_boot, 2, function(x) quantile(x, alpha / 2)),
      obs = effect1, 
      upr = apply(causal_effect1_boot, 2, function(x) quantile(x, 1 - alpha / 2))
    )
    
    CI_beta2_t <- data.frame(
      lwr = apply(causal_effect2_boot, 2, function(x) quantile(x, alpha / 2)),
      obs = effect2, 
      upr = apply(causal_effect2_boot, 2, function(x) quantile(x, 1 - alpha / 2))
    )
    
    # Update ggdata with bootstrap confidence intervals
    ggdata$effect1_low <- CI_beta1_t$lwr
    ggdata$effect1_up <- CI_beta1_t$upr
    ggdata$effect2_low <- CI_beta2_t$lwr
    ggdata$effect2_up <- CI_beta2_t$upr
    
    # Save bootstrap results
    fitRES$CI_beta1_k <- CI_beta1_k
    fitRES$CI_beta2_k <- CI_beta2_k
    fitRES$CI_beta1_t <- CI_beta1_t
    fitRES$CI_beta2_t <- CI_beta2_t
    
  } else {
    if(verbose == TRUE) print("No bootstrap")
    # Standard asymptotic CIs are already calculated above
  }
  
  
  # 8. CREATE PLOTS
  p1 <- ggplot2::ggplot(ggdata,
                        ggplot2::aes(x = time, y = effect1)) +
    
    # Reference line
    ggplot2::geom_hline(
      yintercept = 0,
      linewidth = 0.4,
      linetype = "dashed",
      colour = "grey50"
    ) +
    
    # Confidence band
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = effect1_low, ymax = effect1_up),
      fill = "blue4",
      alpha = 0.2
    ) +
    
    # Estimated effect
    ggplot2::geom_line(
      linewidth = 1.2,
      colour = "blue4"
    ) +
    
    ggplot2::labs(
      x = "Age",
      y = "Time-varying effect",
      title = expression(beta[1](t))
    ) +
    
    ggplot2::coord_cartesian(
      ylim = c(
        min(ggdata$effect1_low) - 0.5 * plotdif,
        max(ggdata$effect1_up)  + 0.5 * plotdif
      )
    ) +
    
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold")
    )
  if (!is.null(true_shape1)) {
    ggdata$true_shape1 <- true_shape1
    
    p1 <- p1 +
      ggplot2::geom_line(
        ggplot2::aes(y = true_shape1),
        linewidth = 1,
        linetype = "longdash",
        colour = "#E15759"
      )
  }
  p2 <- ggplot2::ggplot(ggdata,
                        ggplot2::aes(x = time, y = effect2)) +
    
    # Reference line
    ggplot2::geom_hline(
      yintercept = 0,
      linewidth = 0.4,
      linetype = "dashed",
      colour = "grey50"
    ) +
    
    # Confidence band
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = effect2_low, ymax = effect2_up),
      fill = "blue4",
      alpha = 0.2
    ) +
    
    # Estimated effect
    ggplot2::geom_line(
      linewidth = 1.2,
      colour = "blue4"
    ) +
    
    ggplot2::labs(
      x = "Age",
      y = "Time-varying effect",
      title = expression(beta[2](t))
    ) +
    
    ggplot2::coord_cartesian(
      ylim = c(
        min(ggdata$effect2_low) - 0.5 * plotdif,
        max(ggdata$effect2_up)  + 0.5 * plotdif
      )
    ) +
    
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold")
    )
  if (!is.null(true_shape2)) {
    ggdata$true_shape2 <- true_shape2
    
    p2 <- p2 +
      ggplot2::geom_line(
        ggplot2::aes(y = true_shape2),
        linewidth = 1,
        linetype = "longdash",
        colour = "#E15759"
      )
  }
  
  fitRES$p1 <- p1
  fitRES$p2 <- p2
  
  
  # 9. PERFORMANCE METRICS
  
  
  if (!is.na(X1Ymodel) | !is.na(X2Ymodel)) {
    if (!is.na(X1Ymodel)) {
      SE1 <- (effect1 - ggdata$true_shape1)^2
      fitRES$SE1 <- SE1
      fitRES$MISE1 <- mean(SE1)
      Co1 <- (ggdata$true_shape1 >= ggdata$effect1_low) & 
        (ggdata$true_shape1 <= ggdata$effect1_up)
      fitRES$Coverage_rate1 <- mean(Co1)
    }
    
    if (!is.na(X2Ymodel)) {
      SE2 <- (effect2 - ggdata$true_shape2)^2
      fitRES$SE2 <- SE2
      fitRES$MISE2 <- mean(SE2)
      Co2 <- (ggdata$true_shape2 >= ggdata$effect2_low) & 
        (ggdata$true_shape2 <= ggdata$effect2_up)
      fitRES$Coverage_rate2 <- mean(Co2)
    }
  }
  
  return(fitRES)
}


# ============= SEPARATE UNIVARIABLE FMVMR ESTIMATION # =============
# Separate estimation of time-varying effects for each exposure independently

#' Separate univariable functional MR estimation (internal)
#'
#' Performs separate estimation of time-varying causal effects for each
#' exposure independently with automatic component selection.
#'
#' @param Gmatrix1 Genetic instrument matrix for exposure 1 (N × J1)
#' @param Gmatrix2 Genetic instrument matrix from exposure 2 (N × J2) or NULL for single exposure
#' @param res1 FPCA result for exposure 1 (from fdapace)
#' @param res2 FPCA result for exposure 2 (from fdapace)
#' @param Yvector Outcome vector (length N)
#' @param IDmatch Optional index vector to match rows of Gmatrix, Gmatrix2 (if used), and Yvector (default: 1:N)
#' @param nPC1_selected Fixed number of principal components to retain for exposure 1 (NA = select automatically)
#' @param max_nPC1 Maximum number of principal components to consider for exposure 1 during selection
#' @param nPC2_selected Fixed number of principal components to retain for exposure 2 (NA = select automatically)
#' @param max_nPC2 Maximum number of principal components to consider for exposure 2 during selection
#' @param method Estimation method: "gmm" (Generalized Method of Moments), "cf" (control function), or "cf-lasso" (control function with Lasso)
#' @param basis Basis type for functional representation: "eigenfunction" or "polynomial"
#' @param outcome Outcome type: "continuous" for numeric or "binary" for 0/1 outcomes
#' @param bootstrap Logical; whether to perform bootstrap inference for confidence intervals
#' @param n_B Number of bootstrap iterations (used only if bootstrap = TRUE)
#' @param improvement_threshold Minimum cross-validation improvement required to add an additional principal component
#' @param X1Ymodel Optional: true effect model for X1 on Y (simulation only)
#' @param X2Ymodel Optional: true effect model for X2 on Y (simulation only)
#' @param Gmatrix2 Genetic instrument matrix for exposure 2 (required if separate_G = TRUE)
#' @param num_cores_set Number of CPU cores to use for parallel processing
#' @param verbose Print progress messages and diagnostics during computation
#'
#' @return List with separate estimation results for both exposures
#' @keywords internal
Separate_Multi_MVFMR <- function(Gmatrix1, 
                                 Gmatrix2 = NULL, 
                                 res1 = NA, 
                                 res2 = NA, 
                                 Yvector, 
                                 IDmatch = NA,
                                 nPC1_selected = NA, 
                                 max_nPC1 = NA,
                                 nPC2_selected = NA, 
                                 max_nPC2 = NA,
                                 method = "gmm",
                                 basis = "eigenfunction",
                                 outcome = "continuous",
                                 bootstrap = FALSE,
                                 n_B = 10,
                                 improvement_threshold = 0.01,
                                 X1Ymodel = NA, 
                                 X2Ymodel = NA,
                                 num_cores_set = NA,
                                 verbose = FALSE) {
  
  
  # 1. VALIDATION
  
  
  if (!basis %in% c("eigenfunction", "polynomial")) {
    stop("The 'basis' parameter must be either 'eigenfunction' or 'polynomial'")
  }
  
  if (!method %in% c("gmm", "cf", "cf-lasso")) {
    stop("The 'method' parameter must be one of 'gmm', 'cf', or 'cf-lasso'")
  }
  
  if (!outcome %in% c("continuous", "binary")) {
    stop("The 'outcome' parameter must be either 'continuous' or 'binary'")
  }
  

  G1 <- as.data.frame(Gmatrix1)
  G2 <- as.data.frame(Gmatrix2)
  
  
  # 2. CORE FUNCTION: PROCESS SINGLE EXPOSURE
  
  
  process_single_exposure <- function(Gmatrix, 
                                      res, 
                                      Yvector,
                                      IDmatch,
                                      selected_nPC = NA,
                                      max_nPC = NA,
                                      method = "gmm",
                                      basis = "eigenfunction",
                                      outcome = "continuous",
                                      bootstrap = FALSE,
                                      n_B = 10,
                                      improvement_threshold = 0.01,
                                      true_shape_model = NA,
                                      exposure_name = "X",
                                      X_curves = NULL) {
    
    
    result <- list()
    
    max_npc <- ifelse(is.na(max_nPC) | max_nPC > ncol(res$xiEst), 
                      ncol(res$xiEst), max_nPC)
    
    workGrid_diff <- res$workGrid[3] - res$workGrid[2]
    
    if (is.null(X_curves)) {
      X_curves <- fitted(res, ciOptns = list(alpha = 0.05, kernelType = 'gauss'))$fitted
    }
    
    get_polynomial_basis <- function(grid, n_components) {
      outer(grid, 0:(n_components - 1), `^`)
    }
    
    if (basis == "polynomial") {
      bbb <- get_polynomial_basis(res$workGrid, max_npc)
      phi_transposed <- t(res$phi[-1, 1:max_npc])
    }
    
    # Cross-validation component selection
    n_folds <- 3
    valid_indices <- which(!is.na(IDmatch))
    n_valid <- length(valid_indices)
    
    set.seed(123)
    fold_indices <- sample(rep(1:n_folds, length.out = n_valid))
    
    evaluate_components <- function(nPC) {
      fold_metrics <- numeric(n_folds)
      
      for (fold in 1:n_folds) {
        val_idx <- valid_indices[fold_indices == fold]
        train_idx <- valid_indices[fold_indices != fold]
        
        PC_ <- na.omit(as.matrix(res$xiEst[, 1:nPC]))
        
        Z_train <- Gmatrix[train_idx, , drop = FALSE]
        X_train <- PC_[train_idx, , drop = FALSE]
        Y_train <- Yvector[IDmatch[train_idx]]
        
        Z_val <- Gmatrix[val_idx, , drop = FALSE]
        X_val <- PC_[val_idx, , drop = FALSE]
        Y_val <- Yvector[IDmatch[val_idx]]
        
        if (basis == "polynomial") {
          B <- phi_transposed[1:nPC, ] %*% bbb[-1, 1:nPC] * workGrid_diff
          X_train_model <- X_train %*% B
          X_val_model <- X_val %*% B
        } else {
          X_train_model <- X_train
          X_val_model <- X_val
        }
        
        tryCatch({
          suppressWarnings({
          model_res <- run_model(X = X_train_model, Y = Y_train, Z = Z_train, 
                                 method_type = method, outcome_type = outcome)
          })
          
          if (basis == "polynomial") {
            X_curves_val <- X_curves[val_idx, ]
            val_predictions <- X_curves_val %*% bbb[, 1:nPC] %*% model_res$gmm_est
          } else {
            X_curves_val <- X_curves[val_idx, ]
            val_predictions <- X_curves_val %*% (res$phi)[, 1:nPC] %*% model_res$gmm_est
          }
          
          if (outcome == "continuous") {
            fold_metrics[fold] <- mean((Y_val - val_predictions)^2)
          } else if (outcome == "binary") {
            val_predictions_prob <- 1 / (1 + exp(-val_predictions))
            fold_metrics[fold] <- pROC::auc(pROC::roc(response  = Y_val, predictor = as.numeric(val_predictions_prob), levels    = c(0, 1), direction = "<", quiet     = TRUE))
          }
        }, error = function(e) {
          if (outcome == "continuous") {
            fold_metrics[fold] <- Inf
          } else {
            fold_metrics[fold] <- 0
          }
          warning(paste("Model fitting failed for nPC =", nPC, "in fold", fold))
          warning(paste("Model fitting failed for nPC =", nPC, "in fold", fold, 
                        "| Error:", e$message))
          print(e)
        })
      }
      
      mean(fold_metrics)
    }
    
    is_better <- function(new_metric, current_best) {
      if (outcome == "continuous") {
        return(new_metric < current_best)
      } else {
        return(new_metric > current_best)
      }
    }
    
    # Sequential component selection
    component_results <- data.frame(
      nPC = integer(),
      metric = numeric(),
      step = integer()
    )
    
    best_nPC <- 2
    
    if(verbose == TRUE) print(paste("Evaluating initial model with nPC =", best_nPC, "for", exposure_name))
    best_metric <- evaluate_components(best_nPC)
    
    component_results <- rbind(component_results, data.frame(
      nPC = best_nPC, 
      metric = best_metric,
      step = 0
    ))
    
    max_components <- min(max_npc, 10)
    
    improved <- TRUE
    step <- 1
    
    while (improved && best_nPC < max_components) {
      improved <- FALSE
      
      if (best_nPC < max_npc) {
        if(verbose == TRUE) print(paste("Trying nPC =", best_nPC + 1, "for", exposure_name))
        new_metric <- evaluate_components(best_nPC + 1)
        
        component_results <- rbind(component_results, data.frame(
          nPC = best_nPC + 1, 
          metric = new_metric,
          step = step
        ))
        
        if (outcome == "continuous") {
          improvement <- (best_metric - new_metric) / best_metric
        } else {
          improvement <- (new_metric - best_metric) / best_metric
        }
        
        if (improvement > improvement_threshold) {
          improved <- TRUE
          best_nPC <- best_nPC + 1
          best_metric <- new_metric
          if(verbose == TRUE) print(paste("Improvement found: increasing nPC to", best_nPC, "with metric", best_metric, "for", exposure_name))
        } else {
          if(verbose == TRUE) print(paste("No significant improvement found. Stopping component selection for", exposure_name))
        }
      } else {
        if(verbose == TRUE) print(paste("Reached maximum number of components for", exposure_name))
        improved <- FALSE
      }
      
      step <- step + 1
    }
    
    result$component_selection_results <- component_results
    
    if(!is.na(selected_nPC)) {best_nPC = selected_nPC}
    result$nPC_used <- best_nPC
    
    # Train final model
    PC_ <- na.omit(as.matrix(res$xiEst[, 1:best_nPC]))  
    
    Z_GMMused <- Gmatrix[!is.na(IDmatch), , drop = FALSE]
    X_GMMused <- PC_[!is.na(IDmatch), , drop = FALSE]
    Y_GMMused <- Yvector[IDmatch][!is.na(IDmatch)]
    
    if (basis == "polynomial") {
      B <- phi_transposed[1:best_nPC, ] %*% bbb[-1, 1:best_nPC] * workGrid_diff
      X_model <- X_GMMused %*% B
    } else {
      X_model <- X_GMMused
    }
    
    suppressWarnings({
      gmm_res <- run_model(X = X_model, Y = Y_GMMused, Z = Z_GMMused,
                           method_type = method, outcome_type = outcome)
    })
    
    if (basis == "polynomial") {
      y_estimated <- t(X_curves %*% bbb[, 1:best_nPC] %*% gmm_res$gmm_est)
    } else {
      y_estimated <- t(X_curves %*% (res$phi)[, 1:best_nPC] %*% gmm_res$gmm_est)
    }
    
    if (outcome == "continuous") {
      gmm_res$mse <- mean((Yvector - y_estimated)^2)
      result$final_mse <- gmm_res$mse
    } else if (outcome == "binary") {
      y_estimated_prob <- 1 / (1 + exp(-t(y_estimated)))
      gmm_res$auc <- pROC::auc(pROC::roc(Yvector, as.numeric(y_estimated_prob)))
      result$final_auc <- gmm_res$auc
    }
    
    result$MPCMRest <- gmm_res$gmm_est
    result$MPCMRvar <- gmm_res$variance_matrix
    
    if (!is.na(true_shape_model)) {
      true_shape <- get_true_shape_values(res$workGrid, as.character(true_shape_model))
      
      if (basis == "polynomial") {
        true_beta_k <- c(t((res$phi)[, 1:best_nPC]) %*% bbb[, 1:best_nPC] %*% 
                           solve(t(bbb[, 1:best_nPC]) %*% bbb[, 1:best_nPC]) %*% 
                           t(bbb[, 1:best_nPC]) %*% true_shape)
      } else {
        true_beta_k <- c(t((res$phi)[, 1:best_nPC]) %*% (c(true_shape) * workGrid_diff))
      }
      
      data <- data.frame(
        Beta = gmm_res$gmm_est,
        CI_Lower_Trad = gmm_res$gmm_est - 1.96 * gmm_res$gmm_se,
        CI_Upper_Trad = gmm_res$gmm_est + 1.96 * gmm_res$gmm_se
      )
      data$X_axis <- 1:nrow(data)
      
      plot_beta <- ggplot2::ggplot(data, ggplot2::aes(x = factor(X_axis))) +
        ggplot2::geom_point(ggplot2::aes(y = Beta), size = 2, color = "darkgreen") +
        ggplot2::geom_point(ggplot2::aes(y = true_beta_k, color = "True beta_k"), size = 2) +
        ggplot2::geom_segment(ggplot2::aes(x = X_axis, xend = X_axis, 
                                           y = CI_Lower_Trad, yend = CI_Upper_Trad), size = 1) +
        ggplot2::geom_text(ggplot2::aes(y = Beta, label = round(Beta, 4)), 
                           color = "darkgreen", hjust = -0.4, vjust = -0.1, size = 3) +
        ggplot2::theme_minimal() +
        ggplot2::labs(title = paste(exposure_name, "Beta Coefficients"), y = "Estimate (95% CI)") +
        ggplot2::scale_x_discrete(labels = paste0("Beta_", 1:best_nPC)) +
        ggplot2::scale_color_manual(values = c("Estimated beta_k" = "darkgreen", 
                                               "True beta_k" = "blue")) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(size = 12), 
                       axis.title.x = ggplot2::element_blank(), 
                       legend.title = ggplot2::element_blank())
      
      result$plot_beta <- plot_beta
    }
    
    # Calculate effects
    if (basis == "polynomial") {
      phi <- bbb[, 1:best_nPC]
    } else {
      phi <- (res$phi)[, 1:best_nPC]
    }
    
    pointwise_shape_var <- diag(phi %*% gmm_res$variance_matrix %*% t(phi))
    
    if (basis == "polynomial") {
      effect <- bbb[, 1:best_nPC] %*% gmm_res$gmm_est
    } else {
      effect <- (res$phi)[, 1:best_nPC] %*% gmm_res$gmm_est
    }
    
    ggdata <- data.frame(
      time = res$workGrid,
      effect = effect,
      effect_low = effect - 1.96 * sqrt(pointwise_shape_var),
      effect_up = effect + 1.96 * sqrt(pointwise_shape_var),
      true_shape = NaN
    )
    
    if (!is.na(true_shape_model)) {
      ggdata$true_shape <- true_shape
    }
    
    result$ggdata <- ggdata
    
    # Bootstrap procedure
    if (bootstrap == TRUE) {
      if(verbose == TRUE) print(paste("Starting bootstrap for", exposure_name))
      
      causal_effect_boot <- matrix(NA, nrow = n_B, ncol = nrow(ggdata))
      causal_beta_boot <- matrix(NA, nrow = n_B, ncol = best_nPC)
      
      PC_ <- na.omit(as.matrix(res$xiEst[, 1:best_nPC]))
      
      Z_GMMused <- Gmatrix[!is.na(IDmatch), , drop = FALSE]
      X_GMMused <- PC_[!is.na(IDmatch), , drop = FALSE]
      Y_GMMused <- Yvector[IDmatch][!is.na(IDmatch)]
      
      if (basis == "polynomial") {
        B <- phi_transposed[1:best_nPC, ] %*% bbb[-1, 1:best_nPC] * workGrid_diff
        X_model <- X_GMMused %*% B
      } else {
        X_model <- X_GMMused
      }
      
      pb <- progress::progress_bar$new(
        format = paste("  processing", exposure_name, "bootstrap [:bar] :percent eta: :eta"),
        total = n_B,
        clear = FALSE
      )
      
      for (b in 1:n_B) {
        bootstrap_ind <- sample(1:nrow(X_model), replace = TRUE)
        
        X.b <- X_model[bootstrap_ind, ]
        Y.b <- Y_GMMused[bootstrap_ind]
        Z.b <- Z_GMMused[bootstrap_ind, ]
        
        suppressWarnings({
          gmm_res_boot <- run_model(X = X.b, Y = Y.b, Z = Z.b, 
                                    method_type = method, outcome_type = outcome)
        })
        
        
        causal_beta_boot[b, ] <- gmm_res_boot$gmm_est
        
        if (basis == "polynomial") {
          causal_effect_boot[b, ] <- bbb[, 1:best_nPC] %*% gmm_res_boot$gmm_est
        } else {
          causal_effect_boot[b, ] <- res$phi[, 1:best_nPC] %*% gmm_res_boot$gmm_est
        }
        
        pb$tick()
      }
      
      alpha <- 0.05
      
      CI_beta_k <- data.frame(
        lwr = apply(causal_beta_boot, 2, function(x) quantile(x, alpha / 2)),
        obs = gmm_res$gmm_est, 
        upr = apply(causal_beta_boot, 2, function(x) quantile(x, 1 - alpha / 2))
      )
      
      if (outcome == "binary") {
        CI_beta_k <- exp(CI_beta_k)
      }
      
      CI_beta_t <- data.frame(
        lwr = apply(causal_effect_boot, 2, function(x) quantile(x, alpha / 2)),
        obs = effect, 
        upr = apply(causal_effect_boot, 2, function(x) quantile(x, 1 - alpha / 2))
      )
      
      ggdata$effect_low <- CI_beta_t$lwr
      ggdata$effect_up <- CI_beta_t$upr
      
      result$CI_beta_k <- CI_beta_k
      result$CI_beta_t <- CI_beta_t
    }
    
    # Performance metrics
    if (!is.na(true_shape_model)) {
      SE <- (effect - ggdata$true_shape)^2
      result$SE <- SE
      result$MISE <- mean(SE)
      Co <- (ggdata$true_shape >= ggdata$effect_low) & 
        (ggdata$true_shape <= ggdata$effect_up)
      result$Coverage_rate <- mean(Co)
    }
    
    return(list(result = result, best_nPC = best_nPC, ggdata = ggdata))
  }
  
  
  # 3. HELPER FUNCTION: RUN MODEL
  
  
  run_model <- function(X, Y, Z, method_type, outcome_type) {
    if (method_type == "gmm" && outcome_type == "continuous") {
      return(gmm_lm_onesample(X = X, Y = Y, Z = Z))
    } else if (method_type == "cf" && outcome_type == "binary") {
      return(cf_logit(X = X, Y = Y, Z = Z, use_lasso = FALSE))
    } else if (method_type == "cf-lasso" && outcome_type == "binary") {
      return(cf_logit(X = X, Y = Y, Z = Z, use_lasso = TRUE))
    }
  }
  
  
  # 4. INITIALIZATION
  
  
  fitRES <- list()
  
  X1_curves <- fitted(res1, ciOptns = list(alpha = 0.05, kernelType = 'gauss'))$fitted
  
  if(nrow(G2)>0){
    X2_curves <- fitted(res2, ciOptns = list(alpha = 0.05, kernelType = 'gauss'))$fitted
  }
  
  if (is.na(num_cores_set)) {
    num_cores <- parallel::detectCores() - 1 
  } else {
    num_cores <- num_cores_set
  } 
  
  cl <- parallel::makeCluster(num_cores)
  doParallel::registerDoParallel(cl) 
  
  
  # 5. PROCESS X1
  
  
  print("Processing X1")
  
  X1_results <- process_single_exposure(
    Gmatrix = G1,
    res = res1,
    Yvector = Yvector,
    IDmatch = IDmatch,
    selected_nPC = nPC1_selected,
    max_nPC = max_nPC1,
    method = method,
    basis = basis,
    outcome = outcome,
    bootstrap = bootstrap,
    n_B = n_B,
    improvement_threshold = improvement_threshold,
    true_shape_model = X1Ymodel,
    exposure_name = "X1",
    X_curves = X1_curves
  )
  
  fitRES$nPC_used1 <- X1_results$best_nPC
  fitRES$component_selection_results1 <- X1_results$result$component_selection_results
  fitRES$plot_beta1 <- X1_results$result$plot_beta
  fitRES$MPCMRest1 <- X1_results$result$MPCMRest
  fitRES$MPCMRvar1 <- X1_results$result$MPCMRvar
  fitRES$ggdata1 <- X1_results$ggdata
  
  if (!is.na(X1Ymodel)) {
    fitRES$SE1 <- X1_results$result$SE
    fitRES$MISE1 <- X1_results$result$MISE
    fitRES$Coverage_rate1 <- X1_results$result$Coverage_rate
  }
  
  if (bootstrap) {
    fitRES$CI_beta1_k <- X1_results$result$CI_beta_k
    fitRES$CI_beta1_t <- X1_results$result$CI_beta_t
  }
  
  if (outcome == "continuous") {
    fitRES$final_mse1 <- X1_results$result$final_mse
  } else {
    fitRES$final_auc1 <- X1_results$result$final_auc
  }
  
  
  # 6. PROCESS X2
  
  
  if(nrow(G2)>0){
    print("Processing X2")
    
    X2_results <- process_single_exposure(
      Gmatrix = G2,
      res = res2,
      Yvector = Yvector,
      IDmatch = IDmatch,
      selected_nPC = nPC2_selected,
      max_nPC = max_nPC2,
      method = method,
      basis = basis,
      outcome = outcome,
      bootstrap = bootstrap,
      n_B = n_B,
      improvement_threshold = improvement_threshold,
      true_shape_model = X2Ymodel,
      exposure_name = "X2",
      X_curves = X2_curves
    )
    
    fitRES$nPC_used2 <- X2_results$best_nPC
    fitRES$component_selection_results2 <- X2_results$result$component_selection_results
    fitRES$plot_beta2 <- X2_results$result$plot_beta
    fitRES$MPCMRest2 <- X2_results$result$MPCMRest
    fitRES$MPCMRvar2 <- X2_results$result$MPCMRvar
    fitRES$ggdata2 <- X2_results$ggdata
    
    if (!is.na(X2Ymodel)) {
      fitRES$SE2 <- X2_results$result$SE
      fitRES$MISE2 <- X2_results$result$MISE
      fitRES$Coverage_rate2 <- X2_results$result$Coverage_rate
    }
    
    if (bootstrap) {
      fitRES$CI_beta2_k <- X2_results$result$CI_beta_k
      fitRES$CI_beta2_t <- X2_results$result$CI_beta_t
    }
    
    if (outcome == "continuous") {
      fitRES$final_mse2 <- X2_results$result$final_mse
    } else {
      fitRES$final_auc2 <- X2_results$result$final_auc
    }
  }
  
  
  
  # 7. CREATE PLOTS
  
  
  if(nrow(G2)>0){
    plotdif <- max(c(X1_results$ggdata$effect_up, X2_results$ggdata$effect_up)) - 
      min(c(X1_results$ggdata$effect_low, X2_results$ggdata$effect_low))
  }else{
    plotdif <- max(X1_results$ggdata$effect_up) - min(X1_results$ggdata$effect_low)
  }
  
  
  p1 <- ggplot2::ggplot(X1_results$ggdata,
                        ggplot2::aes(x = time, y = effect)) +
    
    # Reference line
    ggplot2::geom_hline(
      yintercept = 0,
      linewidth = 0.4,
      linetype = "dashed",
      colour = "grey50"
    ) +
    
    # Confidence band
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = effect_low, ymax = effect_up),
      fill = "blue4",
      alpha = 0.2
    ) +
    
    # Estimated effect
    ggplot2::geom_line(
      linewidth = 1.2,
      colour = "blue4"
    ) +
    
    ggplot2::labs(
      x = "Age",
      y = "Time-varying effect",
      title = expression(beta[1](t))
    ) +
    
    ggplot2::coord_cartesian(
      ylim = c(
        min(X1_results$ggdata$effect_low) - 0.5 * plotdif,
        max(X1_results$ggdata$effect_up)  + 0.5 * plotdif
      )
    ) +
    
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold")
    )
  if (!is.na(X1Ymodel)) {
    true_shape1 <- get_true_shape_values(res1$workGrid, as.character(X1Ymodel))
    X1_results$ggdata$true_shape <- true_shape1
    
    p1 <- p1 +
      ggplot2::geom_line(
        ggplot2::aes(y = true_shape),
        linewidth = 1,
        linetype = "longdash",
        colour = "#E15759"
      )
  }
  
  fitRES$p1 <- p1
  if(nrow(G2)>0){
    p2 <- ggplot2::ggplot(X2_results$ggdata,
                          ggplot2::aes(x = time, y = effect)) +
      
      # Reference line at zero
      ggplot2::geom_hline(
        yintercept = 0,
        linewidth = 0.4,
        linetype = "dashed",
        colour = "grey50"
      ) +
      
      # Confidence band
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = effect_low, ymax = effect_up),
        fill = "blue4",
        alpha = 0.2
      ) +
      
      # Estimated effect
      ggplot2::geom_line(
        linewidth = 1.2,
        colour = "blue4"
      ) +
      
      ggplot2::labs(
        x = "Age",
        y = "Time-varying effect",
        title = expression(beta[2](t))
      ) +
      
      ggplot2::coord_cartesian(
        ylim = c(
          min(X2_results$ggdata$effect_low) - 0.5 * plotdif,
          max(X2_results$ggdata$effect_up)  + 0.5 * plotdif
        )
      ) +
      
      ggplot2::theme_minimal(base_size = 13) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(face = "bold"),
        axis.title = ggplot2::element_text(face = "bold")
      )
    if (!is.na(X2Ymodel)) {
      true_shape2 <- get_true_shape_values(res2$workGrid, as.character(X2Ymodel))
      X2_results$ggdata$true_shape <- true_shape2
      
      p2 <- p2 +
        ggplot2::geom_line(
          ggplot2::aes(y = true_shape),
          linewidth = 1,
          linetype = "dash",
          colour = "#E15759"
        )
    }
    
    
    fitRES$p2 <- p2
  }
  
  
  parallel::stopCluster(cl)
  
  return(fitRES)
}

# ============= TWO-SAMPLE ESTIMATION FUNCTIONS # =============
# Functions for two-sample MR using summary statistics
# Exposure data (G, FPCA) is used to get first-stage, then combined with outcome summary stats

#' Two-sample joint multivariable FMR (internal)
#' 
#' @param Gmatrix Genetic instrument matrix from the exposure sample (N × J)
#' @param res1 FPCA result for exposure 1 
#' @param res2 FPCA result for exposure 2 
#' @param by_used Vector of SNP-outcome effect estimates (betas) from the outcome GWAS, length J
#' @param sy_used Vector of standard errors for SNP-outcome effects, length J
#' @param ny_used Sample size of the outcome GWAS
#' @param max_nPC1 Maximum number of principal components to retain for exposure 1 (NA = select automatically)
#' @param max_nPC2 Maximum number of principal components to retain for exposure 2 (NA = select automatically)
#' @param X1Ymodel True effect model for X1 on Y (for simulation only)
#' @param X2Ymodel True effect model for X2 on Y (for simulation only)
#' @param basis Basis type for functional representation: "eigenfunction" or "polynomial"
#' 
#' @return List with separate estimation results for both exposures
AUTOMATIC_Multi_FMVMR_twosample_simple <- function(Gmatrix,
                                                   res1,
                                                   res2,
                                                   by_used,
                                                   sy_used,
                                                   ny_used,
                                                   max_nPC1 = NA,
                                                   max_nPC2 = NA,
                                                   X1Ymodel = NA,
                                                   X2Ymodel = NA,
                                                   basis = "eigenfunction") {
  
  fitRES <- list()
  
  # Component selection (95% variance)
  if (is.na(max_nPC1)) {
    cumvar1 <- cumsum(res1$lambda) / sum(res1$lambda)
    nPC1 <- which(cumvar1 >= 0.95)[1]
  } else {
    nPC1 <- max_nPC1
  }
  
  if (is.na(max_nPC2)) {
    cumvar2 <- cumsum(res2$lambda) / sum(res2$lambda)
    nPC2 <- which(cumvar2 >= 0.95)[1]
  } else {
    nPC2 <- max_nPC2
  }
  
  print(paste("Selected: nPC1 =", nPC1, ", nPC2 =", nPC2))
  
  # Get PCs
  PC1 <- na.omit(as.matrix(res1$xiEst[, 1:nPC1]))
  PC2 <- na.omit(as.matrix(res2$xiEst[, 1:nPC2]))
  
  # Calculate first-stage bx (SNP -> PC associations)
  J <- ncol(Gmatrix)
  K <- nPC1 + nPC2
  bx <- matrix(NA, J, K)
  
  for (j in 1:J) {
    for (k in 1:nPC1) {
      fit <- stats::lm(PC1[, k] ~ Gmatrix[, j])
      bx[j, k] <- coef(fit)[2]
    }
    for (k in 1:nPC2) {
      fit <- stats::lm(PC2[, k] ~ Gmatrix[, j])
      bx[j, nPC1 + k] <- coef(fit)[2]
    }
  }
  
  # Two-sample GMM
  gmm_res <- gmm_twosample_simple(
    bx = bx,
    by = by_used,
    sy = sy_used,
    ny = ny_used
  )
  
  fitRES$nPC1_used <- nPC1
  fitRES$nPC2_used <- nPC2
  fitRES$MPCMRest <- gmm_res$gmm_est
  fitRES$MPCMRvar <- gmm_res$variance_matrix
  fitRES$Q_stat <- gmm_res$Q_stat
  fitRES$Q_pval <- gmm_res$Q_pval
  
  # Calculate time-varying effects
  effect1 <- (res1$phi)[, 1:nPC1] %*% gmm_res$gmm_est[1:nPC1]
  effect2 <- (res2$phi)[, 1:nPC2] %*% gmm_res$gmm_est[(nPC1 + 1):K]
  
  phi <- cbind((res1$phi)[, 1:nPC1], (res2$phi)[, 1:nPC2])
  pointwise_var <- diag(phi %*% gmm_res$variance_matrix %*% t(phi))
  
  ggdata <- data.frame(
    time = res1$workGrid,
    effect1 = effect1,
    effect1_low = effect1 - 1.96 * sqrt(pointwise_var),
    effect1_up = effect1 + 1.96 * sqrt(pointwise_var),
    effect2 = effect2,
    effect2_low = effect2 - 1.96 * sqrt(pointwise_var),
    effect2_up = effect2 + 1.96 * sqrt(pointwise_var)
  )
  
  # Add true shapes if provided
  if (!is.na(X1Ymodel)) {
    ggdata$true_shape1 <- get_true_shape_values(res1$workGrid, as.character(X1Ymodel))
  }
  if (!is.na(X2Ymodel)) {
    ggdata$true_shape2 <- get_true_shape_values(res2$workGrid, as.character(X2Ymodel))
  }
  
  fitRES$ggdata <- ggdata
  
  # Plots
  p1 <- ggplot2::ggplot(ggdata, ggplot2::aes(time, effect1)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.5, linetype = 2, col = 'grey') +
    ggplot2::geom_line(ggplot2::aes(time, effect1), linewidth = 1) +
    ggplot2::geom_line(ggplot2::aes(time, effect1_low), linewidth = 1, linetype = 2) +
    ggplot2::geom_line(ggplot2::aes(time, effect1_up), linewidth = 1, linetype = 2) +
    ggplot2::labs(x = 'Age', y = 'Effect') +
    ggplot2::ggtitle("Beta1(t) - Two-Sample") +
    ggplot2::theme_bw()
  
  p2 <- ggplot2::ggplot(ggdata, ggplot2::aes(time, effect2)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.5, linetype = 2, col = 'grey') +
    ggplot2::geom_line(ggplot2::aes(time, effect2), linewidth = 1) +
    ggplot2::geom_line(ggplot2::aes(time, effect2_low), linewidth = 1, linetype = 2) +
    ggplot2::geom_line(ggplot2::aes(time, effect2_up), linewidth = 1, linetype = 2) +
    ggplot2::labs(x = 'Age', y = 'Effect') +
    ggplot2::ggtitle("Beta2(t) - Two-Sample") +
    ggplot2::theme_bw()
  
  if (!is.na(X1Ymodel)) {
    p1 <- p1 + ggplot2::geom_line(ggplot2::aes(time, true_shape1), linewidth = 1, col = 'blue')
  }
  if (!is.na(X2Ymodel)) {
    p2 <- p2 + ggplot2::geom_line(ggplot2::aes(time, true_shape2), linewidth = 1, col = 'blue')
  }
  
  fitRES$p1 <- p1
  fitRES$p2 <- p2
  
  # Performance if available
  if (!is.na(X1Ymodel)) {
    SE1 <- (effect1 - ggdata$true_shape1)^2
    fitRES$MISE1 <- mean(SE1)
    fitRES$Coverage_rate1 <- mean((ggdata$true_shape1 >= ggdata$effect1_low) & 
                                    (ggdata$true_shape1 <= ggdata$effect1_up))
  }
  if (!is.na(X2Ymodel)) {
    SE2 <- (effect2 - ggdata$true_shape2)^2
    fitRES$MISE2 <- mean(SE2)
    fitRES$Coverage_rate2 <- mean((ggdata$true_shape2 >= ggdata$effect2_low) & 
                                    (ggdata$true_shape2 <= ggdata$effect2_up))
  }
  
  return(fitRES)
}


#' Separate univariable two-sample FMR (internal)
#' @param Gmatrix1 Genetic instrument matrix from exposure 1 (N × J1)
#' @param Gmatrix2 Genetic instrument matrix from exposure 2 (N × J2) or NULL
#' @param res1 FPCA result for exposure 1 
#' @param res2 FPCA result for exposure 2 
#' @param by_used1 Vector of SNP-outcome effect estimates (betas for X1) from the outcome GWAS, length J
#' @param by_used2 Vector of SNP-outcome effect estimates (betas for X2) from the outcome GWAS, length J or NULL
#' @param sy_used1 Vector of standard errors for SNP-outcome effects for X1, length J
#' @param sy_used2 Vector of standard errors for SNP-outcome effects for X2, length J or NULL
#' @param ny_used Sample size of the outcome GWAS
#' @param max_nPC1 Maximum number of principal components to retain for exposure 1 (NA = select automatically)
#' @param max_nPC2 Maximum number of principal components to retain for exposure 2 (NA = select automatically)
#' @param X1Ymodel True effect model for X1 on Y (for simulation only)
#' @param X2Ymodel True effect model for X2 on Y (for simulation only)
#' @param basis Basis type for functional representation: "eigenfunction" or "polynomial"
#' 
#' @return List with separate estimation results for both exposures
Separate_Multi_FMVMR_twosample_simple <- function(Gmatrix1,
                                                  Gmatrix2 = NULL,
                                                  res1,
                                                  res2,
                                                  by_used1,
                                                  by_used2 = NULL,
                                                  sy_used1,
                                                  sy_used2 = NULL,
                                                  ny_used,
                                                  max_nPC1 = NA,
                                                  max_nPC2 = NA,
                                                  X1Ymodel = NA,
                                                  X2Ymodel = NA, 
                                                  basis = "eigenfunction") {
  
  fitRES <- list()
  
  # === EXPOSURE 1 ===
  print("*** Processing Exposure 1 ***")
  
  if (is.na(max_nPC1)) {
    cumvar1 <- cumsum(res1$lambda) / sum(res1$lambda)
    nPC1 <- which(cumvar1 >= 0.95)[1]
  } else {
    nPC1 <- max_nPC1
  }
  
  PC1 <- na.omit(as.matrix(res1$xiEst[, 1:nPC1]))
  J1 <- ncol(Gmatrix1)
  bx1 <- matrix(NA, J1, nPC1)
  
  for (j in 1:J1) {
    for (k in 1:nPC1) {
      fit <- stats::lm(PC1[, k] ~ Gmatrix1[, j])
      bx1[j, k] <- coef(fit)[2]
    }
  }
  
  gmm_res1 <- gmm_twosample_simple(
    bx = bx1,
    by = by_used1,
    sy = sy_used1,
    ny = ny_used
  )
  
  fitRES$nPC_used1 <- nPC1
  fitRES$MPCMRest1 <- gmm_res1$gmm_est
  fitRES$MPCMRvar1 <- gmm_res1$variance_matrix
  
  effect1 <- (res1$phi)[, 1:nPC1] %*% gmm_res1$gmm_est
  pointwise_var1 <- diag((res1$phi)[, 1:nPC1] %*% gmm_res1$variance_matrix %*% 
                           t((res1$phi)[, 1:nPC1]))
  
  ggdata1 <- data.frame(
    time = res1$workGrid,
    effect = effect1,
    effect_low = effect1 - 1.96 * sqrt(pointwise_var1),
    effect_up = effect1 + 1.96 * sqrt(pointwise_var1)
  )
  
  if (!is.na(X1Ymodel)) {
    ggdata1$true_shape <- get_true_shape_values(res1$workGrid, as.character(X1Ymodel))
    SE1 <- (effect1 - ggdata1$true_shape)^2
    fitRES$MISE1 <- mean(SE1)
    fitRES$Coverage_rate1 <- mean((ggdata1$true_shape >= ggdata1$effect_low) & 
                                    (ggdata1$true_shape <= ggdata1$effect_up))
  }
  
  fitRES$ggdata1 <- ggdata1
  
  p1 <- ggplot2::ggplot(ggdata1,
                        ggplot2::aes(x = time, y = effect)) +
    
    # Reference line at zero
    ggplot2::geom_hline(
      yintercept = 0,
      linewidth = 0.4,
      linetype = "dashed",
      colour = "grey50"
    ) +
    
    # Confidence band
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = effect_low, ymax = effect_up),
      fill = "blue4",
      alpha = 0.2
    ) +
    
    # Estimated effect
    ggplot2::geom_line(
      linewidth = 1.2,
      colour = "blue4"
    ) +
    
    ggplot2::labs(
      x = "Age",
      y = "Effect",
      title = expression(beta[1](t))
    ) +
    
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold")
    )
  if (!is.na(X1Ymodel)) {
    p1 <- p1 +
      ggplot2::geom_line(
        ggplot2::aes(y = true_shape),
        linewidth = 1,
        linetype = "longdash",
        colour = "#E15759"
      )
  }
  
  
  fitRES$p1 <- p1
  
  # === EXPOSURE 2 (if provided) ===
  if (!is.null(Gmatrix2) && !is.null(by_used2)) {
    print("*** Processing Exposure 2 ***")
    
    if (is.na(max_nPC2)) {
      cumvar2 <- cumsum(res2$lambda) / sum(res2$lambda)
      nPC2 <- which(cumvar2 >= 0.95)[1]
    } else {
      nPC2 <- max_nPC2
    }
    
    PC2 <- na.omit(as.matrix(res2$xiEst[, 1:nPC2]))
    J2 <- ncol(Gmatrix2)
    bx2 <- matrix(NA, J2, nPC2)
    
    for (j in 1:J2) {
      for (k in 1:nPC2) {
        fit <- stats::lm(PC2[, k] ~ Gmatrix2[, j])
        bx2[j, k] <- coef(fit)[2]
      }
    }
    
    gmm_res2 <- gmm_twosample_simple(
      bx = bx2,
      by = by_used2,
      sy = sy_used2,
      ny = ny_used
    )
    
    fitRES$nPC_used2 <- nPC2
    fitRES$MPCMRest2 <- gmm_res2$gmm_est
    fitRES$MPCMRvar2 <- gmm_res2$variance_matrix
    
    effect2 <- (res2$phi)[, 1:nPC2] %*% gmm_res2$gmm_est
    pointwise_var2 <- diag((res2$phi)[, 1:nPC2] %*% gmm_res2$variance_matrix %*% 
                             t((res2$phi)[, 1:nPC2]))
    
    ggdata2 <- data.frame(
      time = res2$workGrid,
      effect = effect2,
      effect_low = effect2 - 1.96 * sqrt(pointwise_var2),
      effect_up = effect2 + 1.96 * sqrt(pointwise_var2)
    )
    
    if (!is.na(X2Ymodel)) {
      ggdata2$true_shape <- get_true_shape_values(res2$workGrid, as.character(X2Ymodel))
      SE2 <- (effect2 - ggdata2$true_shape)^2
      fitRES$MISE2 <- mean(SE2)
      fitRES$Coverage_rate2 <- mean((ggdata2$true_shape >= ggdata2$effect_low) & 
                                      (ggdata2$true_shape <= ggdata2$effect_up))
    }
    
    fitRES$ggdata2 <- ggdata2
    
    p2 <- ggplot2::ggplot(ggdata2,
                          ggplot2::aes(x = time, y = effect)) +
      
      # Reference line at zero
      ggplot2::geom_hline(
        yintercept = 0,
        linewidth = 0.4,
        linetype = "dashed",
        colour = "grey50"
      ) +
      
      # Confidence band
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = effect_low, ymax = effect_up),
        fill = "blue4",
        alpha = 0.2
      ) +
      
      # Estimated effect
      ggplot2::geom_line(
        linewidth = 1.2,
        colour = "blue4"
      ) +
      
      ggplot2::labs(
        x = "Age",
        y = "Effect",
        title = expression(beta[2](t))
      ) +
      
      ggplot2::theme_minimal(base_size = 13) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(face = "bold"),
        axis.title = ggplot2::element_text(face = "bold")
      )
    if (!is.na(X2Ymodel)) {
      p2 <- p2 +
        ggplot2::geom_line(
          ggplot2::aes(y = true_shape),
          linewidth = 1,
          linetype = "longdash",
          colour = "#E15759"
        )
    }
    
    
    fitRES$p2 <- p2
  }
  
  return(fitRES)
}