# ============= AUTOMATIC MULTIVARIABLE MVFMR - JOINT ESTIMATION # =============
# Joint estimation of multiple time-varying exposures with automatic 
# component selection via cross-validation

#' Automatic Multivariable Functional MR with joint estimation (internal)
#'
#' Core function that performs joint estimation of time-varying causal effects
#' from m correlated exposures using automatic component selection.
#'
#' @param Gmatrix Genetic instrument matrix (N × J)
#' @param res_list List of length m of FPCA results, one per exposure
#' @param Yvector Outcome vector
#' @param IDmatch Optional index vector to match rows of Gmatrix and Yvector (default: 1:N)
#' @param nPC_selected Length-m vector: fixed number of principal components to retain per exposure (NA = select automatically)
#' @param max_nPC Length-m vector: maximum number of principal components to consider per exposure during selection
#' @param X_true Length-m list: optional true X curves per exposure (simulation only), NULL entries allowed
#' @param method Estimation method: "gmm" (Generalized Method of Moments), "cf" (control function), or "cf-lasso" (control function with Lasso)
#' @param basis Basis type for functional representation: "eigenfunction" or "polynomial"
#' @param outcome Outcome type: "continuous" for numeric or "binary" for 0/1 outcomes
#' @param bootstrap Logical; whether to perform bootstrap inference for confidence intervals
#' @param n_B Number of bootstrap iterations (used only if bootstrap = TRUE)
#' @param improvement_threshold Minimum cross-validation improvement required to add an additional principal component
#' @param XYmodels Length-m vector: optional true effect model for each exposure on Y (simulation only)
#' @param num_cores_set Number of CPU cores to use for parallel processing
#' @param verbose Print progress messages and diagnostics during computation
#'
#' @return List with estimation results, selected components, performance metrics
#' @keywords internal
AUTOMATIC_Multi_MVFMR <- function(Gmatrix,
                                  res_list,
                                  Yvector,
                                  IDmatch = NA,
                                  nPC_selected = NA,
                                  max_nPC = NA,
                                  X_true = NULL,
                                  method = "gmm",
                                  basis = "eigenfunction",
                                  outcome = "continuous",
                                  bootstrap = FALSE,
                                  n_B = 10,
                                  improvement_threshold = 0.01,
                                  XYmodels = NA,
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

  m <- length(res_list)
  nPC_selected <- recycle_arg(nPC_selected, m, default = NA)
  max_nPC <- recycle_arg(max_nPC, m, default = NA)
  XYmodels <- recycle_arg(XYmodels, m, default = NA)
  if (is.null(X_true)) X_true <- vector("list", m)

  max_npc_vec <- sapply(seq_len(m), function(k) {
    ifelse(is.na(max_nPC[k]) | max_nPC[k] > ncol(res_list[[k]]$xiEst),
           ncol(res_list[[k]]$xiEst), max_nPC[k])
  })

  J <- ncol(Gmatrix)

  workGrid_diff_vec <- sapply(seq_len(m), function(k) res_list[[k]]$workGrid[3] - res_list[[k]]$workGrid[2])

  # Prepare fitted curves
  X_curves_list <- vector("list", m)
  for (k in seq_len(m)) {
    if (is.null(X_true[[k]]) || any(is.na(X_true[[k]]))) {
      X_curves_list[[k]] <- fitted(res_list[[k]], ciOptns = list(alpha = 0.05, kernelType = 'gauss'))$fitted
    } else {
      X_curves_list[[k]] <- t(apply(X_true[[k]], 1, function(x) {
        spline(x = seq(1, 50, length.out = ncol(X_true[[k]])),
               y = x, xout = res_list[[k]]$workGrid)$y
      }))
    }
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
    bbb_list <- lapply(seq_len(m), function(k) get_polynomial_basis(res_list[[k]]$workGrid, max_npc_vec[k]))
    phi_transposed_list <- lapply(seq_len(m), function(k) t(res_list[[k]]$phi[-1, 1:max_npc_vec[k]]))
  }

  n_folds <- 3
  if(verbose == TRUE) print("***Feature selection starting***")

  valid_indices <- which(!is.na(IDmatch))
  n_valid <- length(valid_indices)

  set.seed(123)
  fold_indices <- sample(rep(1:n_folds, length.out = n_valid))

  # Helper: evaluate component combination
  evaluate_components <- function(nPC_vec) {
    fold_metrics <- numeric(n_folds)
    offsets_local <- compute_offsets(nPC_vec)

    for (fold in 1:n_folds) {
      val_idx <- valid_indices[fold_indices == fold]
      train_idx <- valid_indices[fold_indices != fold]

      PC_list <- lapply(seq_len(m), function(k) na.omit(as.matrix(res_list[[k]]$xiEst[, 1:nPC_vec[k]])))
      PC_ <- do.call(cbind, PC_list)

      Z_train <- Gmatrix[train_idx, , drop = FALSE]
      X_train <- PC_[train_idx, , drop = FALSE]
      Y_train <- Yvector[IDmatch[train_idx]]

      Z_val <- Gmatrix[val_idx, , drop = FALSE]
      X_val <- PC_[val_idx, , drop = FALSE]
      Y_val <- Yvector[IDmatch[val_idx]]

      if (basis == "polynomial") {
        B_list <- lapply(seq_len(m), function(k) {
          phi_transposed_list[[k]][1:nPC_vec[k], ] %*% bbb_list[[k]][-1, 1:nPC_vec[k]] * workGrid_diff_vec[k]
        })

        X_train_model <- do.call(cbind, lapply(seq_len(m), function(k) PC_list[[k]][train_idx, ] %*% B_list[[k]]))
        X_val_model <- do.call(cbind, lapply(seq_len(m), function(k) PC_list[[k]][val_idx, ] %*% B_list[[k]]))
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
          val_predictions <- Reduce(`+`, lapply(seq_len(m), function(k) {
            X_curves_list[[k]][val_idx, ] %*% bbb_list[[k]][, 1:nPC_vec[k]] %*% gmm_res$gmm_est[block_idx(offsets_local, k)]
          }))
        } else {
          val_predictions <- Reduce(`+`, lapply(seq_len(m), function(k) {
            X_curves_list[[k]][val_idx, ] %*% (res_list[[k]]$phi)[, 1:nPC_vec[k]] %*% gmm_res$gmm_est[block_idx(offsets_local, k)]
          }))
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
        warning(paste("Model fitting failed for nPC =", paste(nPC_vec, collapse = ","), "in fold", fold))
      })
    }

    mean(fold_metrics)
  }

  # Helper: record a component-selection trial
  record_step <- function(component_results, nPC_vec, metric, step) {
    row <- as.data.frame(as.list(nPC_vec))
    names(row) <- paste0("nPC", seq_len(m))
    row$metric <- metric
    row$step <- step
    rbind(component_results, row)
  }

  # Sequential component selection
  component_results <- as.data.frame(matrix(nrow = 0, ncol = m + 2))
  names(component_results) <- c(paste0("nPC", seq_len(m)), "metric", "step")

  best_nPC_vec <- rep(2, m)

  best_metric <- suppressWarnings({evaluate_components(best_nPC_vec)})

  component_results <- record_step(component_results, best_nPC_vec, best_metric, 0)

  max_components <- min(sum(max_npc_vec), 6 * m)

  improved <- TRUE
  step <- 1

  while (improved && any(best_nPC_vec < max_npc_vec) && sum(best_nPC_vec) < max_components) {
    improved <- FALSE

    candidate_metric <- rep(if (outcome == "continuous") Inf else 0, m)

    for (k in seq_len(m)) {
      if (best_nPC_vec[k] < max_npc_vec[k]) {
        trial_nPC_vec <- best_nPC_vec
        trial_nPC_vec[k] <- trial_nPC_vec[k] + 1
        if(verbose == TRUE) print(paste("Trying nPC =", paste(trial_nPC_vec, collapse = ",")))
        candidate_metric[k] <- suppressWarnings({evaluate_components(trial_nPC_vec)})

        component_results <- record_step(component_results, trial_nPC_vec, candidate_metric[k], step)
      }
    }

    if (outcome == "continuous") {
      improvement <- (best_metric - candidate_metric) / best_metric
    } else {
      improvement <- (candidate_metric - best_metric) / best_metric
    }

    if (max(improvement) > improvement_threshold) {
      improved <- TRUE

      k_star <- which.max(improvement)
      best_nPC_vec[k_star] <- best_nPC_vec[k_star] + 1
      best_metric <- candidate_metric[k_star]
      if(verbose == TRUE) print(paste("Improvement found: increasing nPC", k_star, "to", best_nPC_vec[k_star], "with metric", best_metric))
    } else {
      if(verbose == TRUE) print("No significant improvement found. Stopping component selection.")
    }

    step <- step + 1
  }

  if(verbose == TRUE) print(paste("Final selected components:", paste(paste0("nPC", seq_len(m), " = ", best_nPC_vec), collapse = ", ")))
  if(verbose == TRUE) print(paste("Final metric:", best_metric))

  fitRES$component_selection_results <- component_results


  # 6. FINAL MODEL WITH SELECTED COMPONENTS

  if (all(!is.na(nPC_selected))) {
    best_nPC_vec <- nPC_selected
  }

  offsets <- compute_offsets(best_nPC_vec)

  PC_list <- lapply(seq_len(m), function(k) na.omit(as.matrix(res_list[[k]]$xiEst[, 1:best_nPC_vec[k]])))
  PC_ <- do.call(cbind, PC_list)

  Z_GMMused <- Gmatrix[!is.na(IDmatch), , drop = FALSE]
  X_GMMused <- PC_[!is.na(IDmatch), , drop = FALSE]
  Y_GMMused <- Yvector[IDmatch][!is.na(IDmatch)]

  if (basis == "polynomial") {
    B_list <- lapply(seq_len(m), function(k) {
      phi_transposed_list[[k]][1:best_nPC_vec[k], ] %*% bbb_list[[k]][-1, 1:best_nPC_vec[k]] * workGrid_diff_vec[k]
    })
    X_model <- do.call(cbind, lapply(seq_len(m), function(k) PC_list[[k]] %*% B_list[[k]]))
  } else {
    X_model <- X_GMMused
  }

  suppressWarnings({
    gmm_res <- run_model(X = X_model, Y = Y_GMMused, Z = Z_GMMused,
                         method_type = method, outcome_type = outcome)
  })


  if (basis == "polynomial") {
    y_estimated <- t(Reduce(`+`, lapply(seq_len(m), function(k) {
      X_curves_list[[k]] %*% bbb_list[[k]][, 1:best_nPC_vec[k]] %*% gmm_res$gmm_est[block_idx(offsets, k)]
    })))
  } else {
    y_estimated <- t(Reduce(`+`, lapply(seq_len(m), function(k) {
      X_curves_list[[k]] %*% (res_list[[k]]$phi)[, 1:best_nPC_vec[k]] %*% gmm_res$gmm_est[block_idx(offsets, k)]
    })))
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


  fitRES$nPC_used <- best_nPC_vec
  fitRES$offsets <- offsets

  # Add true shapes if available
  if (any(!is.na(XYmodels))) {
    true_shape_list <- lapply(seq_len(m), function(k) get_true_shape_values(res_list[[k]]$workGrid, as.character(XYmodels[k])))
  } else {
    true_shape_list <- vector("list", m)
  }
  all_true_shapes_known <- all(!sapply(true_shape_list, is.null))

  # Prepare data for plots
  data <- data.frame(
    Beta = gmm_res$gmm_est,
    CI_Lower_Trad = gmm_res$gmm_est - 1.96 * gmm_res$gmm_se,
    CI_Upper_Trad = gmm_res$gmm_est + 1.96 * gmm_res$gmm_se
  )
  data$X_axis <- 1:nrow(data)

  final_phi <- do.call(cbind, lapply(seq_len(m), function(k) (res_list[[k]]$phi)[, 1:best_nPC_vec[k]]))

  beta_labels <- unlist(lapply(seq_len(m), function(k) paste0("Beta", k, "_", 1:best_nPC_vec[k])))

  if (all_true_shapes_known) {
    if (basis == "polynomial") {
      true_beta_k <- unlist(lapply(seq_len(m), function(k) {
        t((res_list[[k]]$phi)[, 1:best_nPC_vec[k]]) %*% bbb_list[[k]][, 1:best_nPC_vec[k]] %*%
          solve(t(bbb_list[[k]][, 1:best_nPC_vec[k]]) %*% bbb_list[[k]][, 1:best_nPC_vec[k]]) %*%
          t(bbb_list[[k]][, 1:best_nPC_vec[k]]) %*% true_shape_list[[k]]
      }))
    } else {
      true_beta_k <- unlist(lapply(seq_len(m), function(k) {
        t((res_list[[k]]$phi)[, 1:best_nPC_vec[k]]) %*% (c(true_shape_list[[k]]) * workGrid_diff_vec[k])
      }))
    }

    plot <- ggplot2::ggplot(data, ggplot2::aes(x = factor(X_axis))) +
      ggplot2::geom_point(ggplot2::aes(y = Beta, color = "Estimated beta_k"), size = 2) +
      ggplot2::geom_point(ggplot2::aes(y = true_beta_k, color = "True beta_k"), size = 2) +
      ggplot2::geom_segment(ggplot2::aes(x = X_axis, xend = X_axis,
                                         y = CI_Lower_Trad, yend = CI_Upper_Trad,
                                         color = "Traditional CI"), linewidth = 1) +
      ggplot2::geom_text(ggplot2::aes(y = Beta, label = round(Beta, 4)),
                         color = "darkgreen", hjust = -0.4, vjust = -0.1, size = 3) +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = "", y = "Estimate (95% CI)") +
      ggplot2::scale_color_manual(values = c("Estimated beta_k" = "darkgreen",
                                              "True beta_k" = "blue")) +
      ggplot2::scale_x_discrete(labels = beta_labels) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(size = 12),
                     axis.title.x = ggplot2::element_blank(),
                     legend.title = ggplot2::element_blank())

  } else {
    plot <- ggplot2::ggplot(data, ggplot2::aes(x = factor(X_axis))) +
      ggplot2::geom_point(ggplot2::aes(y = Beta, color = "Estimated beta_k"), size = 2) +
      ggplot2::geom_segment(ggplot2::aes(x = X_axis, xend = X_axis,
                                         y = CI_Lower_Trad, yend = CI_Upper_Trad,
                                         color = "Traditional CI"), linewidth = 1) +
      ggplot2::geom_text(ggplot2::aes(y = Beta, label = round(Beta, 4)),
                         color = "darkgreen", hjust = -0.4, vjust = -0.1, size = 3) +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = "", y = "Estimate (95% CI)") +
      ggplot2::scale_color_manual(values = c("Estimated beta_k" = "darkgreen",
                                              "True beta_k" = "blue")) +
      ggplot2::scale_x_discrete(labels = beta_labels) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(size = 12),
                     axis.title.x = ggplot2::element_blank(),
                     legend.title = ggplot2::element_blank())
  }

  fitRES$plot_beta <- plot
  fitRES$MPCMRest <- gmm_res$gmm_est
  fitRES$MPCMRvar <- gmm_res$variance_matrix

  # Define projection matrices
  if (basis == "polynomial") {
    phi_list <- lapply(seq_len(m), function(k) bbb_list[[k]][, 1:best_nPC_vec[k]])
  } else {
    phi_list <- lapply(seq_len(m), function(k) (res_list[[k]]$phi)[, 1:best_nPC_vec[k]])
  }
  phi <- do.call(cbind, phi_list)

  pointwise_shape_var <- diag(phi %*% gmm_res$variance_matrix %*% t(phi))
  fitRES$pointwise_shape_var <- pointwise_shape_var
  fitRES$pointwise_estimates <- phi %*% gmm_res$gmm_est

  # Prepare effect curves
  effect_list <- lapply(seq_len(m), function(k) phi_list[[k]] %*% gmm_res$gmm_est[block_idx(offsets, k)])

  ggdata_list <- vector("list", m)
  for (k in seq_len(m)) {
    ggdata_list[[k]] <- data.frame(
      time = res_list[[k]]$workGrid,
      effect = effect_list[[k]],
      effect_low = effect_list[[k]] - 1.96 * sqrt(pointwise_shape_var),
      effect_up = effect_list[[k]] + 1.96 * sqrt(pointwise_shape_var),
      true_shape = NaN
    )
  }

  if (all_true_shapes_known) {
    for (k in seq_len(m)) {
      ggdata_list[[k]]$true_shape <- true_shape_list[[k]]
    }
  }

  fitRES$ggdata <- ggdata_list
  plotdif <- max(sapply(ggdata_list, function(g) max(g$effect_up))) -
    min(sapply(ggdata_list, function(g) min(g$effect_low)))

  # Bootstrap procedure
  if(bootstrap == TRUE) {
    if(verbose == TRUE) print("Starting bootstrap")

    # Setup storage for bootstrap results
    causal_effect_boot <- lapply(seq_len(m), function(k) matrix(NA, nrow = n_B, ncol = nrow(ggdata_list[[k]])))
    causal_beta_boot <- lapply(seq_len(m), function(k) matrix(NA, nrow = n_B, ncol = best_nPC_vec[k]))

    # Prepare data for bootstrap
    PC_list <- lapply(seq_len(m), function(k) na.omit(as.matrix(res_list[[k]]$xiEst[, 1:best_nPC_vec[k]])))
    PC_ <- do.call(cbind, PC_list)

    Z_GMMused <- Gmatrix[!is.na(IDmatch), , drop = FALSE]
    X_GMMused <- PC_[!is.na(IDmatch), , drop = FALSE]
    Y_GMMused <- Yvector[IDmatch][!is.na(IDmatch)]

    # Apply basis transformation if needed
    if(basis == "polynomial") {
      # Create transformation matrices
      B_list <- lapply(seq_len(m), function(k) {
        phi_transposed_list[[k]][1:best_nPC_vec[k], ] %*% bbb_list[[k]][-1, 1:best_nPC_vec[k]] * workGrid_diff_vec[k]
      })

      # Apply transformations
      X_model <- do.call(cbind, lapply(seq_len(m), function(k) PC_list[[k]] %*% B_list[[k]]))
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
      for (k in seq_len(m)) {
        idx <- block_idx(offsets, k)
        causal_beta_boot[[k]][b, ] <- gmm_res_boot$gmm_est[idx]

        # Transform to get effects based on basis type
        if(basis == "polynomial") {
          causal_effect_boot[[k]][b, ] <- bbb_list[[k]][, 1:best_nPC_vec[k]] %*% gmm_res_boot$gmm_est[idx]
        } else {
          causal_effect_boot[[k]][b, ] <- res_list[[k]]$phi[, 1:best_nPC_vec[k]] %*% gmm_res_boot$gmm_est[idx]
        }
      }

      pb$tick()
    }

    alpha = 0.05

    CI_beta_k_list <- vector("list", m)
    CI_beta_t_list <- vector("list", m)

    for (k in seq_len(m)) {
      idx <- block_idx(offsets, k)

      # Calculate confidence intervals for beta coefficients
      CI_beta_k_list[[k]] <- data.frame(
        lwr = apply(causal_beta_boot[[k]], 2, function(x) quantile(x, alpha / 2)),
        obs = gmm_res$gmm_est[idx],
        upr = apply(causal_beta_boot[[k]], 2, function(x) quantile(x, 1 - alpha / 2))
      )

      # Apply exponentiation for binary outcomes
      if(outcome == "binary") {
        CI_beta_k_list[[k]] <- exp(CI_beta_k_list[[k]])
      }

      # Calculate confidence intervals for effects
      CI_beta_t_list[[k]] <- data.frame(
        lwr = apply(causal_effect_boot[[k]], 2, function(x) quantile(x, alpha / 2)),
        obs = effect_list[[k]],
        upr = apply(causal_effect_boot[[k]], 2, function(x) quantile(x, 1 - alpha / 2))
      )

      # Update ggdata with bootstrap confidence intervals
      ggdata_list[[k]]$effect_low <- CI_beta_t_list[[k]]$lwr
      ggdata_list[[k]]$effect_up <- CI_beta_t_list[[k]]$upr
    }

    fitRES$ggdata <- ggdata_list

    # Save bootstrap results
    fitRES$CI_beta_k <- CI_beta_k_list
    fitRES$CI_beta_t <- CI_beta_t_list

  } else {
    if(verbose == TRUE) print("No bootstrap")
    # Standard asymptotic CIs are already calculated above
  }


  # 8. CREATE PLOTS
  fitRES$p <- vector("list", m)

  for (k in seq_len(m)) {
    ggdata_k <- fitRES$ggdata[[k]]

    p_k <- ggplot2::ggplot(ggdata_k,
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
        ggplot2::aes(ymin = effect_low, ymax = effect_up, fill = "95% CI"),
        alpha = 0.2
      ) +

      # Estimated effect
      ggplot2::geom_line(
        ggplot2::aes(colour = "Estimated effect"),
        linewidth = 1.2
      ) +

      ggplot2::labs(
        x = "Age",
        y = "Time-varying effect",
        title = bquote(beta[.(k)](t))
      ) +

      ggplot2::coord_cartesian(
        ylim = c(
          min(ggdata_k$effect_low) - 0.5 * plotdif,
          max(ggdata_k$effect_up)  + 0.5 * plotdif
        )
      ) +

      ggplot2::scale_colour_manual(name = NULL, values = c("Estimated effect" = "blue4", "True effect" = "#E15759")) +
      ggplot2::scale_fill_manual(name = NULL, values = c("95% CI" = "blue4")) +

      ggplot2::theme_minimal(base_size = 13) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(face = "bold"),
        axis.title = ggplot2::element_text(face = "bold"),
        legend.position = "bottom",
        legend.title = ggplot2::element_blank()
      )
    if (!is.null(true_shape_list[[k]])) {
      ggdata_k$true_shape <- true_shape_list[[k]]
      fitRES$ggdata[[k]] <- ggdata_k

      p_k <- p_k +
        ggplot2::geom_line(
          data = ggdata_k,
          ggplot2::aes(y = true_shape, colour = "True effect"),
          linewidth = 1,
          linetype = "longdash"
        )
    }

    fitRES$p[[k]] <- p_k
  }


  # 9. PERFORMANCE METRICS


  if (any(!is.na(XYmodels))) {
    fitRES$SE <- vector("list", m)
    fitRES$MISE <- vector("list", m)
    fitRES$Coverage_rate <- vector("list", m)

    for (k in seq_len(m)) {
      if (!is.na(XYmodels[k])) {
        ggdata_k <- fitRES$ggdata[[k]]
        SE_k <- (effect_list[[k]] - ggdata_k$true_shape)^2
        fitRES$SE[[k]] <- SE_k
        fitRES$MISE[[k]] <- mean(SE_k)
        Co_k <- (ggdata_k$true_shape >= ggdata_k$effect_low) &
          (ggdata_k$true_shape <= ggdata_k$effect_up)
        fitRES$Coverage_rate[[k]] <- mean(Co_k)
      }
    }
  }

  return(fitRES)
}


# ============= SEPARATE UNIVARIABLE FMVMR ESTIMATION # =============
# Separate estimation of time-varying effects for each exposure independently

#' Separate univariable functional MR estimation (internal)
#'
#' Performs separate estimation of time-varying causal effects for each
#' of m exposures independently with automatic component selection.
#'
#' @param Gmatrix_list List of length m of genetic instrument matrices, one per exposure (N x J_k)
#' @param res_list List of length m of FPCA results (from fdapace), one per exposure
#' @param Yvector Outcome vector (length N)
#' @param IDmatch Optional index vector to match rows of the Gmatrix_list entries and Yvector (default: 1:N)
#' @param nPC_selected Length-m vector: fixed number of principal components to retain per exposure (NA = select automatically)
#' @param max_nPC Length-m vector: maximum number of principal components to consider per exposure during selection
#' @param X_true Length-m list: optional true X curves per exposure (simulation only), NULL entries allowed
#' @param method Estimation method: "gmm" (Generalized Method of Moments), "cf" (control function), or "cf-lasso" (control function with Lasso)
#' @param basis Basis type for functional representation: "eigenfunction" or "polynomial"
#' @param outcome Outcome type: "continuous" for numeric or "binary" for 0/1 outcomes
#' @param bootstrap Logical; whether to perform bootstrap inference for confidence intervals
#' @param n_B Number of bootstrap iterations (used only if bootstrap = TRUE)
#' @param improvement_threshold Minimum cross-validation improvement required to add an additional principal component
#' @param XYmodels Length-m vector: optional true effect model for each exposure on Y (simulation only)
#' @param num_cores_set Number of CPU cores to use for parallel processing
#' @param verbose Print progress messages and diagnostics during computation
#'
#' @return List with separate estimation results for each of the m exposures
#' @keywords internal
Separate_Multi_MVFMR <- function(Gmatrix_list,
                                 res_list,
                                 Yvector,
                                 IDmatch = NA,
                                 nPC_selected = NA,
                                 max_nPC = NA,
                                 X_true = NULL,
                                 method = "gmm",
                                 basis = "eigenfunction",
                                 outcome = "continuous",
                                 bootstrap = FALSE,
                                 n_B = 10,
                                 improvement_threshold = 0.01,
                                 XYmodels = NA,
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


  m <- length(Gmatrix_list)
  G_list <- lapply(Gmatrix_list, as.data.frame)
  nPC_selected <- recycle_arg(nPC_selected, m, default = NA)
  max_nPC <- recycle_arg(max_nPC, m, default = NA)
  XYmodels <- recycle_arg(XYmodels, m, default = NA)
  if (is.null(X_true)) X_true <- vector("list", m)


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
                                           y = CI_Lower_Trad, yend = CI_Upper_Trad), linewidth = 1) +
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

  # Prepare fitted curves
  X_curves_list <- vector("list", m)
  for (k in seq_len(m)) {
    if (is.null(X_true[[k]]) || any(is.na(X_true[[k]]))) {
      X_curves_list[[k]] <- fitted(res_list[[k]], ciOptns = list(alpha = 0.05, kernelType = 'gauss'))$fitted
    } else {
      X_curves_list[[k]] <- t(apply(X_true[[k]], 1, function(x) {
        spline(x = seq(1, 50, length.out = ncol(X_true[[k]])),
               y = x, xout = res_list[[k]]$workGrid)$y
      }))
    }
  }

  if (is.na(num_cores_set)) {
    num_cores <- parallel::detectCores() - 1
  } else {
    num_cores <- num_cores_set
  }

  cl <- parallel::makeCluster(num_cores)
  doParallel::registerDoParallel(cl)


  # 5. PROCESS EACH EXPOSURE


  exposure_results <- vector("list", m)

  for (k in seq_len(m)) {
    print(paste0("Processing X", k))

    exposure_results[[k]] <- process_single_exposure(
      Gmatrix = G_list[[k]],
      res = res_list[[k]],
      Yvector = Yvector,
      IDmatch = IDmatch,
      selected_nPC = nPC_selected[k],
      max_nPC = max_nPC[k],
      method = method,
      basis = basis,
      outcome = outcome,
      bootstrap = bootstrap,
      n_B = n_B,
      improvement_threshold = improvement_threshold,
      true_shape_model = XYmodels[k],
      exposure_name = paste0("X", k),
      X_curves = X_curves_list[[k]]
    )

    fitRES$nPC_used[[k]] <- exposure_results[[k]]$best_nPC
    fitRES$component_selection_results[[k]] <- exposure_results[[k]]$result$component_selection_results
    fitRES$plot_beta[[k]] <- exposure_results[[k]]$result$plot_beta
    fitRES$MPCMRest[[k]] <- exposure_results[[k]]$result$MPCMRest
    fitRES$MPCMRvar[[k]] <- exposure_results[[k]]$result$MPCMRvar
    fitRES$ggdata[[k]] <- exposure_results[[k]]$ggdata

    if (!is.na(XYmodels[k])) {
      fitRES$SE[[k]] <- exposure_results[[k]]$result$SE
      fitRES$MISE[[k]] <- exposure_results[[k]]$result$MISE
      fitRES$Coverage_rate[[k]] <- exposure_results[[k]]$result$Coverage_rate
    }

    if (bootstrap) {
      fitRES$CI_beta_k[[k]] <- exposure_results[[k]]$result$CI_beta_k
      fitRES$CI_beta_t[[k]] <- exposure_results[[k]]$result$CI_beta_t
    }

    if (outcome == "continuous") {
      fitRES$final_mse[[k]] <- exposure_results[[k]]$result$final_mse
    } else {
      fitRES$final_auc[[k]] <- exposure_results[[k]]$result$final_auc
    }
  }



  # 6. CREATE PLOTS


  effect_up_all <- sapply(exposure_results, function(r) max(r$ggdata$effect_up))
  effect_low_all <- sapply(exposure_results, function(r) min(r$ggdata$effect_low))
  plotdif <- max(effect_up_all) - min(effect_low_all)

  fitRES$p <- vector("list", m)

  for (k in seq_len(m)) {
    p_k <- ggplot2::ggplot(exposure_results[[k]]$ggdata,
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
        ggplot2::aes(ymin = effect_low, ymax = effect_up, fill = "95% CI"),
        alpha = 0.2
      ) +

      # Estimated effect
      ggplot2::geom_line(
        ggplot2::aes(colour = "Estimated effect"),
        linewidth = 1.2
      ) +

      ggplot2::labs(
        x = "Age",
        y = "Time-varying effect",
        title = bquote(beta[.(k)](t))
      ) +

      ggplot2::coord_cartesian(
        ylim = c(
          min(exposure_results[[k]]$ggdata$effect_low) - 0.5 * plotdif,
          max(exposure_results[[k]]$ggdata$effect_up)  + 0.5 * plotdif
        )
      ) +

      ggplot2::scale_colour_manual(name = NULL, values = c("Estimated effect" = "blue4", "True effect" = "#E15759")) +
      ggplot2::scale_fill_manual(name = NULL, values = c("95% CI" = "blue4")) +

      ggplot2::theme_minimal(base_size = 13) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(face = "bold"),
        axis.title = ggplot2::element_text(face = "bold"),
        legend.position = "bottom",
        legend.title = ggplot2::element_blank()
      )
    if (!is.na(XYmodels[k])) {
      true_shape_k <- get_true_shape_values(res_list[[k]]$workGrid, as.character(XYmodels[k]))
      exposure_results[[k]]$ggdata$true_shape <- true_shape_k

      p_k <- p_k +
        ggplot2::geom_line(
          data = exposure_results[[k]]$ggdata,
          ggplot2::aes(y = true_shape, colour = "True effect"),
          linewidth = 1,
          linetype = "longdash"
        )
    }

    fitRES$p[[k]] <- p_k
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
#' @param res_list List of length m of FPCA results, one per exposure
#' @param by_used Vector of SNP-outcome effect estimates (betas) from the outcome GWAS, length J
#' @param sy_used Vector of standard errors for SNP-outcome effects, length J
#' @param ny_used Sample size of the outcome GWAS
#' @param max_nPC Length-m vector: maximum number of principal components to retain per exposure (NA = select automatically)
#' @param XYmodels Length-m vector: true effect model for each exposure on Y (for simulation only)
#' @param basis Basis type for functional representation: "eigenfunction" or "polynomial"
#'
#' @return List with separate estimation results for each of the m exposures
AUTOMATIC_Multi_FMVMR_twosample_simple <- function(Gmatrix,
                                                   res_list,
                                                   by_used,
                                                   sy_used,
                                                   ny_used,
                                                   max_nPC = NA,
                                                   XYmodels = NA,
                                                   basis = "eigenfunction") {

  m <- length(res_list)
  max_nPC <- recycle_arg(max_nPC, m, default = NA)
  XYmodels <- recycle_arg(XYmodels, m, default = NA)

  fitRES <- list()

  # Component selection (95% variance)
  nPC_vec <- integer(m)
  for (k in seq_len(m)) {
    if (is.na(max_nPC[k])) {
      cumvar <- cumsum(res_list[[k]]$lambda) / sum(res_list[[k]]$lambda)
      nPC_vec[k] <- which(cumvar >= 0.95)[1]
    } else {
      nPC_vec[k] <- max_nPC[k]
    }
  }

  print(paste("Selected:", paste(paste0("nPC", seq_len(m), " = ", nPC_vec), collapse = ", ")))

  # Get PCs
  PC_list <- lapply(seq_len(m), function(k) na.omit(as.matrix(res_list[[k]]$xiEst[, 1:nPC_vec[k]])))

  # Calculate first-stage bx (SNP -> PC associations)
  J <- ncol(Gmatrix)
  offsets <- compute_offsets(nPC_vec)
  K <- sum(nPC_vec)
  bx <- matrix(NA, J, K)

  for (j in 1:J) {
    for (kexp in seq_len(m)) {
      for (k in 1:nPC_vec[kexp]) {
        fit <- stats::lm(PC_list[[kexp]][, k] ~ Gmatrix[, j])
        bx[j, offsets[kexp] + k] <- coef(fit)[2]
      }
    }
  }

  # Two-sample GMM
  gmm_res <- gmm_twosample_simple(
    bx = bx,
    by = by_used,
    sy = sy_used,
    ny = ny_used
  )

  fitRES$nPC_used <- nPC_vec
  fitRES$offsets <- offsets
  fitRES$MPCMRest <- gmm_res$gmm_est
  fitRES$MPCMRvar <- gmm_res$variance_matrix
  fitRES$Q_stat <- gmm_res$Q_stat
  fitRES$Q_pval <- gmm_res$Q_pval

  # Calculate time-varying effects
  phi_list <- lapply(seq_len(m), function(k) (res_list[[k]]$phi)[, 1:nPC_vec[k]])
  effect_list <- lapply(seq_len(m), function(k) phi_list[[k]] %*% gmm_res$gmm_est[block_idx(offsets, k)])

  phi <- do.call(cbind, phi_list)
  pointwise_var <- diag(phi %*% gmm_res$variance_matrix %*% t(phi))

  ggdata_list <- vector("list", m)
  for (k in seq_len(m)) {
    ggdata_k <- data.frame(
      time = res_list[[k]]$workGrid,
      effect = effect_list[[k]],
      effect_low = effect_list[[k]] - 1.96 * sqrt(pointwise_var),
      effect_up = effect_list[[k]] + 1.96 * sqrt(pointwise_var)
    )

    # Add true shape if provided
    if (!is.na(XYmodels[k])) {
      ggdata_k$true_shape <- get_true_shape_values(res_list[[k]]$workGrid, as.character(XYmodels[k]))
    }

    ggdata_list[[k]] <- ggdata_k
  }

  fitRES$ggdata <- ggdata_list

  # Plots
  fitRES$p <- vector("list", m)
  for (k in seq_len(m)) {
    ggdata_k <- ggdata_list[[k]]

    p_k <- ggplot2::ggplot(ggdata_k, ggplot2::aes(time, effect)) +
      ggplot2::geom_hline(yintercept = 0, linewidth = 0.5, linetype = 2, col = 'grey') +
      ggplot2::geom_line(ggplot2::aes(time, effect, colour = "Estimated effect", linetype = "Estimated effect"), linewidth = 1) +
      ggplot2::geom_line(ggplot2::aes(time, effect_low, colour = "95% CI", linetype = "95% CI"), linewidth = 1) +
      ggplot2::geom_line(ggplot2::aes(time, effect_up, colour = "95% CI", linetype = "95% CI"), linewidth = 1) +
      ggplot2::labs(x = 'Age', y = 'Effect') +
      ggplot2::ggtitle(paste0("Beta", k, "(t) - Two-Sample")) +
      ggplot2::scale_colour_manual(name = NULL, values = c("Estimated effect" = "black", "95% CI" = "black", "True effect" = "blue")) +
      ggplot2::scale_linetype_manual(name = NULL, values = c("Estimated effect" = "solid", "95% CI" = "dashed", "True effect" = "solid")) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "bottom", legend.title = ggplot2::element_blank())

    if (!is.na(XYmodels[k])) {
      p_k <- p_k + ggplot2::geom_line(ggplot2::aes(time, true_shape, colour = "True effect", linetype = "True effect"), linewidth = 1)
    }

    fitRES$p[[k]] <- p_k
  }

  # Performance if available
  fitRES$MISE <- vector("list", m)
  fitRES$Coverage_rate <- vector("list", m)
  for (k in seq_len(m)) {
    if (!is.na(XYmodels[k])) {
      ggdata_k <- ggdata_list[[k]]
      SE_k <- (effect_list[[k]] - ggdata_k$true_shape)^2
      fitRES$MISE[[k]] <- mean(SE_k)
      fitRES$Coverage_rate[[k]] <- mean((ggdata_k$true_shape >= ggdata_k$effect_low) &
                                          (ggdata_k$true_shape <= ggdata_k$effect_up))
    }
  }

  return(fitRES)
}


#' Separate univariable two-sample FMR (internal)
#' @param Gmatrix_list List of length m of genetic instrument matrices, one per exposure (N x J_k)
#' @param res_list List of length m of FPCA results, one per exposure
#' @param by_used_list List of length m of SNP-outcome effect estimate vectors, one per exposure
#' @param sy_used_list List of length m of SNP-outcome standard error vectors, one per exposure
#' @param ny_used Sample size of the outcome GWAS
#' @param max_nPC Length-m vector: maximum number of principal components to retain per exposure (NA = select automatically)
#' @param XYmodels Length-m vector: true effect model for each exposure on Y (for simulation only)
#' @param basis Basis type for functional representation: "eigenfunction" or "polynomial"
#'
#' @return List with separate estimation results for each of the m exposures
Separate_Multi_FMVMR_twosample_simple <- function(Gmatrix_list,
                                                  res_list,
                                                  by_used_list,
                                                  sy_used_list,
                                                  ny_used,
                                                  max_nPC = NA,
                                                  XYmodels = NA,
                                                  basis = "eigenfunction") {

  m <- length(Gmatrix_list)
  max_nPC <- recycle_arg(max_nPC, m, default = NA)
  XYmodels <- recycle_arg(XYmodels, m, default = NA)

  fitRES <- list()
  fitRES$nPC_used <- vector("list", m)
  fitRES$MPCMRest <- vector("list", m)
  fitRES$MPCMRvar <- vector("list", m)
  fitRES$ggdata <- vector("list", m)
  fitRES$p <- vector("list", m)
  fitRES$MISE <- vector("list", m)
  fitRES$Coverage_rate <- vector("list", m)

  # === PROCESS EACH EXPOSURE ===
  for (k in seq_len(m)) {
    print(paste0("*** Processing Exposure ", k, " ***"))

    if (is.na(max_nPC[k])) {
      cumvar_k <- cumsum(res_list[[k]]$lambda) / sum(res_list[[k]]$lambda)
      nPC_k <- which(cumvar_k >= 0.95)[1]
    } else {
      nPC_k <- max_nPC[k]
    }

    PC_k <- na.omit(as.matrix(res_list[[k]]$xiEst[, 1:nPC_k]))
    J_k <- ncol(Gmatrix_list[[k]])
    bx_k <- matrix(NA, J_k, nPC_k)

    for (j in 1:J_k) {
      for (kk in 1:nPC_k) {
        fit <- stats::lm(PC_k[, kk] ~ Gmatrix_list[[k]][, j])
        bx_k[j, kk] <- coef(fit)[2]
      }
    }

    gmm_res_k <- gmm_twosample_simple(
      bx = bx_k,
      by = by_used_list[[k]],
      sy = sy_used_list[[k]],
      ny = ny_used
    )

    fitRES$nPC_used[[k]] <- nPC_k
    fitRES$MPCMRest[[k]] <- gmm_res_k$gmm_est
    fitRES$MPCMRvar[[k]] <- gmm_res_k$variance_matrix

    effect_k <- (res_list[[k]]$phi)[, 1:nPC_k] %*% gmm_res_k$gmm_est
    pointwise_var_k <- diag((res_list[[k]]$phi)[, 1:nPC_k] %*% gmm_res_k$variance_matrix %*%
                             t((res_list[[k]]$phi)[, 1:nPC_k]))

    ggdata_k <- data.frame(
      time = res_list[[k]]$workGrid,
      effect = effect_k,
      effect_low = effect_k - 1.96 * sqrt(pointwise_var_k),
      effect_up = effect_k + 1.96 * sqrt(pointwise_var_k)
    )

    if (!is.na(XYmodels[k])) {
      ggdata_k$true_shape <- get_true_shape_values(res_list[[k]]$workGrid, as.character(XYmodels[k]))
      SE_k <- (effect_k - ggdata_k$true_shape)^2
      fitRES$MISE[[k]] <- mean(SE_k)
      fitRES$Coverage_rate[[k]] <- mean((ggdata_k$true_shape >= ggdata_k$effect_low) &
                                      (ggdata_k$true_shape <= ggdata_k$effect_up))
    }

    fitRES$ggdata[[k]] <- ggdata_k

    p_k <- ggplot2::ggplot(ggdata_k,
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
        ggplot2::aes(ymin = effect_low, ymax = effect_up, fill = "95% CI"),
        alpha = 0.2
      ) +

      # Estimated effect
      ggplot2::geom_line(
        ggplot2::aes(colour = "Estimated effect"),
        linewidth = 1.2
      ) +

      ggplot2::labs(
        x = "Age",
        y = "Effect",
        title = bquote(beta[.(k)](t))
      ) +

      ggplot2::scale_colour_manual(name = NULL, values = c("Estimated effect" = "blue4", "True effect" = "#E15759")) +
      ggplot2::scale_fill_manual(name = NULL, values = c("95% CI" = "blue4")) +

      ggplot2::theme_minimal(base_size = 13) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(face = "bold"),
        axis.title = ggplot2::element_text(face = "bold"),
        legend.position = "bottom",
        legend.title = ggplot2::element_blank()
      )
    if (!is.na(XYmodels[k])) {
      p_k <- p_k +
        ggplot2::geom_line(
          ggplot2::aes(y = true_shape, colour = "True effect"),
          linewidth = 1,
          linetype = "longdash"
        )
    }

    fitRES$p[[k]] <- p_k
  }

  return(fitRES)
}