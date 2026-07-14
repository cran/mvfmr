# ============= MVFMR - MAIN FUNCTIONS # =============

#' Joint Multivariable Functional Mendelian Randomization
#'
#' @param G Genetic instrument matrix (N x J)
#' @param fpca_results List of length m of FPCA objects from fdapace, one per exposure
#' @param Y Outcome vector
#' @param outcome_type Type of outcome: "continuous" for numeric outcomes, "binary" for 0/1 outcomes
#' @param method Estimation method: "gmm" (Generalized Method of Moments), "cf" (control function), or "cf-lasso" (control function with Lasso)
#' @param nPC Fixed number of principal components to retain per exposure (length 1 or m; NA = select automatically)
#' @param max_nPC Maximum number of principal components to retain per exposure (length 1 or m; NA = automatically determined)
#' @param improvement_threshold Minimum cross-validation improvement required to add an additional principal component
#' @param bootstrap Whether to compute confidence intervals using bootstrap resampling
#' @param n_bootstrap Number of bootstrap replicates (only used if bootstrap = TRUE)
#' @param n_cores Number of CPU cores to use for parallel computations
#' @param true_effects Length-m vector of true effect model codes, one per exposure (simulation only)
#' @param X_true Length-m list of true X curves, one per exposure (simulation only)
#' @param verbose Print progress and diagnostic messages during computation
#'
#' @return mvfmr object with:
#' \itemize{
#'   \item coefficients - Estimated beta coefficients
#'   \item vcov - Variance-covariance matrix
#'   \item effects - List of length m with the estimated time-varying effect curves
#'   \item nPC_used - Components selected per exposure
#'   \item diagnostics - F-statistics, instrument diagnostics
#'   \item performance - MISE, coverage (if true effects provided)
#' }
#' @examples
#' set.seed(1)
#' sim_data <- getX_multi_exposure(N = 60, J = 8, nSparse = 5, n_exposures = 2)
#' outcome_data <- getY_multi_exposure(sim_data, XYmodels = c("2", "8"))
#' fpca_results <- lapply(sim_data$exposures, function(exp_k) {
#'   fdapace::FPCA(exp_k$Ly_sim, exp_k$Lt_sim,
#'                 list(dataType = "Sparse", error = TRUE, verbose = FALSE))
#' })
#' result <- mvfmr(
#'   G = sim_data$details$G,
#'   fpca_results = fpca_results,
#'   Y = outcome_data$Y,
#'   max_nPC = c(2, 2),
#'   n_cores = 1,
#'   verbose = FALSE
#' )
#' coef(result)
#' @export
mvfmr <- function(G,
                  fpca_results,
                  Y,
                  outcome_type = c("continuous", "binary"),
                  method = c("gmm", "cf", "cf-lasso"),
                  nPC = NA,
                  max_nPC = NA,
                  improvement_threshold = 0.001,
                  bootstrap = FALSE,
                  n_bootstrap = 100,
                  n_cores = parallel::detectCores() - 1,
                  true_effects = NULL,
                  X_true = NULL,
                  verbose = FALSE) {

  outcome_type <- match.arg(outcome_type)
  method <- match.arg(method)

  # ---- Check outcome if binary ----
  if (outcome_type == "binary") {
    if (!is.numeric(Y)) {
      stop("For binary outcome, Y must be numeric and coded as 0/1.")
    }

    Y_unique <- unique(na.omit(Y))

    if (!all(Y_unique %in% c(0, 1))) {
      stop(
        "For binary outcome, Y must contain only values 0 and 1. "
      )
    }
  }

  # Extract FPCA objects
  if (!is.list(fpca_results) || length(fpca_results) < 1) {
    stop("fpca_results must be a list of m FPCA objects")
  }

  m <- length(fpca_results)

  # Extract true effects if provided
  XYmodels <- recycle_arg(true_effects, m, default = NA)
  X_true_list <- if (!is.null(X_true)) X_true else vector("list", m)

  N <- nrow(G)
  IDmatch <- 1:N

  if (verbose) {
    cat("\n=== Functional Multivariable MR ===\n")
    cat("Sample size:", N, "\n")
    cat("Exposures:", m, "\n")
    cat("Outcome:", outcome_type, "\n")
    cat("Method:", method, "\n\n")
  }

  # Call original AUTOMATIC_Multi_MVFMR function
  result <- AUTOMATIC_Multi_MVFMR(
    Gmatrix = G,
    res_list = fpca_results,
    Yvector = Y,
    IDmatch = IDmatch,
    nPC_selected = NA,
    max_nPC = max_nPC,
    X_true = X_true_list,
    method = method,
    basis = "eigenfunction",
    outcome = outcome_type,
    bootstrap = bootstrap,
    n_B = n_bootstrap,
    improvement_threshold = improvement_threshold,
    XYmodels = XYmodels,
    num_cores_set = n_cores,
    verbose = verbose
  )

  # Reformat output to package standard
  out <- structure(
    list(
      coefficients = result$MPCMRest,
      vcov = result$MPCMRvar,
      effects = lapply(result$ggdata, function(g) g$effect),
      time_grid = result$ggdata[[1]]$time,
      confidence_intervals = list(
        lower = lapply(result$ggdata, function(g) g$effect_low),
        upper = lapply(result$ggdata, function(g) g$effect_up)
      ),
      nPC_used = result$nPC_used,
      offsets = result$offsets,
      performance = if (!is.null(result$MISE)) {
        list(MISE = result$MISE, Coverage = result$Coverage_rate)
      } else NULL,
      plots = list(effects = result$p, plot_beta = result$plot_beta),
      raw_result = result,
      n_exposures = m,
      outcome_type = outcome_type,
      method = method,
      n_observations = N
    ),
    class = "mvfmr"
  )

  if (verbose) cat("=== Estimation complete ===\n\n")

  return(out)
}

#' Separate Univariable Functional Mendelian Randomization
#'
#' @param G_list List of length m of genetic instrument matrices, one per exposure. Use a
#'   list of length 1 to analyze a single exposure.
#' @param fpca_results List of length m of FPCA objects, same length as G_list
#' @param Y Outcome vector
#' @param outcome_type Type of outcome: "continuous" for numeric outcomes, "binary" for 0/1 outcomes
#' @param method Estimation method: "gmm" (Generalized Method of Moments), "cf" (control function), or "cf-lasso" (control function with Lasso)
#' @param nPC Fixed number of principal components to retain per exposure (length 1 or m; NA = select automatically)
#' @param max_nPC Maximum number of principal components to retain per exposure (length 1 or m; NA = automatically determined)
#' @param improvement_threshold Minimum cross-validation improvement required to add an additional principal component
#' @param bootstrap Whether to compute confidence intervals using bootstrap resampling
#' @param n_bootstrap Number of bootstrap replicates (only used if bootstrap = TRUE)
#' @param n_cores Number of CPU cores to use for parallel computations
#' @param true_effects Length-m vector of true effect model codes, one per exposure (simulation only)
#' @param X_true Length-m list of true X curves, one per exposure (simulation only)
#' @param verbose Print progress and diagnostic messages during computation
#'
#' @return mvfmr_separate object
#' @examples
#' set.seed(1)
#' sim_data <- getX_multi_exposure(N = 60, J = 8, nSparse = 5, n_exposures = 2)
#' outcome_data <- getY_multi_exposure(sim_data, XYmodels = c("2", "8"))
#' fpca_results <- lapply(sim_data$exposures, function(exp_k) {
#'   fdapace::FPCA(exp_k$Ly_sim, exp_k$Lt_sim,
#'                 list(dataType = "Sparse", error = TRUE, verbose = FALSE))
#' })
#' result <- mvfmr_separate(
#'   G_list = list(sim_data$details$G, sim_data$details$G),
#'   fpca_results = fpca_results,
#'   Y = outcome_data$Y,
#'   max_nPC = c(2, 2),
#'   n_cores = 1,
#'   verbose = FALSE
#' )
#' coef(result, exposure = 1)
#' @export
mvfmr_separate <- function(G_list,
                           fpca_results,
                           Y,
                           outcome_type = c("continuous", "binary"),
                           method = c("gmm", "cf", "cf-lasso"),
                           nPC = NA,
                           max_nPC = NA,
                           improvement_threshold = 0.001,
                           bootstrap = FALSE,
                           n_bootstrap = 100,
                           n_cores = parallel::detectCores() - 1,
                           true_effects = NULL,
                           X_true = NULL,
                           verbose = FALSE) {

  outcome_type <- match.arg(outcome_type)
  method <- match.arg(method)

  # ---- Check outcome if binary ----
  if (outcome_type == "binary") {
    if (!is.numeric(Y)) {
      stop("For binary outcome, Y must be numeric and coded as 0/1.")
    }

    Y_unique <- unique(na.omit(Y))

    if (!all(Y_unique %in% c(0, 1))) {
      stop(
        "For binary outcome, Y must contain only values 0 and 1. "
      )
    }
  }

  if (!is.list(G_list) || length(G_list) < 1) {
    stop("G_list must be a list of m genetic instrument matrices")
  }

  m <- length(G_list)

  if (!is.list(fpca_results) || length(fpca_results) != m) {
    stop("fpca_results must be a list with the same length as G_list (", m, ")")
  }

  # Extract true effects if provided
  XYmodels <- recycle_arg(true_effects, m, default = NA)
  X_true_list <- if (!is.null(X_true)) X_true else vector("list", m)

  N <- nrow(G_list[[1]])
  IDmatch <- 1:N

  if (verbose) {
    cat("\n=== Separate Univariable MR ===\n")
    cat("Exposures:", m, "\n")
    cat("Sample size:", N, "\n\n")
  }

  # Call original Separate_Multi_MVFMR function
  result <- Separate_Multi_MVFMR(
    Gmatrix_list = G_list,
    res_list = fpca_results,
    nPC_selected = nPC,
    max_nPC = max_nPC,
    X_true = X_true_list,
    Yvector = Y,
    IDmatch = IDmatch,
    method = method,
    basis = "eigenfunction",
    outcome = outcome_type,
    bootstrap = bootstrap,
    n_B = n_bootstrap,
    improvement_threshold = improvement_threshold,
    XYmodels = XYmodels,
    num_cores_set = n_cores,
    verbose = verbose
  )

  # Reformat output
  exposures <- lapply(seq_len(m), function(k) {
    list(
      coefficients = result$MPCMRest[[k]],
      vcov = result$MPCMRvar[[k]],
      effect = result$ggdata[[k]]$effect,
      nPC_used = result$nPC_used[[k]],
      performance = if (!is.null(result$MISE[[k]])) {
        list(MISE = result$MISE[[k]], Coverage = result$Coverage_rate[[k]])
      } else NULL
    )
  })

  out <- structure(
    list(
      exposures = exposures,
      plots = list(effects = result$p),
      raw_result = result,
      n_exposures = m,
      separate_instruments = TRUE,
      outcome_type = outcome_type,
      method = method
    ),
    class = "mvfmr_separate"
  )

  if (verbose) cat("=== Univariable estimation complete ===\n\n")

  return(out)
}



# ============= MVFMR - TWO SAMPLE # =============
# Functions for two-sample MR design using summary statistics
#' Two-Sample Joint Multivariable Functional MR
#'
#' Joint estimation using outcome GWAS summary statistics.
#' Simplified approach: only needs by, sy, ny (not individual outcome data).
#'
#' @param G_exposure Genetic instrument matrix from the exposure sample (N × J)
#' @param fpca_results List of length m of FPCA objects, one per exposure
#' @param by_outcome Vector of SNP-outcome effect estimates (betas) from the outcome GWAS, length J
#' @param sy_outcome VVector of standard errors for SNP-outcome effects, length J
#' @param ny_outcome Sample size of the outcome GWAS
#' @param max_nPC Maximum number of principal components to retain per exposure (length 1 or m; NA = automatically determined)
#' @param true_effects Length-m vector of true effect model codes, one per exposure (simulation only)
#' @param verbose Print progress messages and diagnostics during computation
#'
#' @return fmvmr_twosample object
#' @examples
#' set.seed(1)
#' sim_data <- getX_multi_exposure(N = 60, J = 8, nSparse = 5, n_exposures = 2)
#' outcome_data <- getY_multi_exposure(sim_data, XYmodels = c("2", "8"))
#' fpca_results <- lapply(sim_data$exposures, function(exp_k) {
#'   fdapace::FPCA(exp_k$Ly_sim, exp_k$Lt_sim,
#'                 list(dataType = "Sparse", error = TRUE, verbose = FALSE))
#' })
#' # Simulate outcome GWAS summary statistics for the two-sample design
#' by_outcome <- sapply(1:8, function(j) {
#'   coef(lm(outcome_data$Y ~ sim_data$details$G[, j]))[2]
#' })
#' sy_outcome <- sapply(1:8, function(j) {
#'   summary(lm(outcome_data$Y ~ sim_data$details$G[, j]))$coefficients[2, 2]
#' })
#' result <- fmvmr_twosample(
#'   G_exposure = sim_data$details$G,
#'   fpca_results = fpca_results,
#'   by_outcome = by_outcome,
#'   sy_outcome = sy_outcome,
#'   ny_outcome = 60,
#'   max_nPC = c(2, 2),
#'   verbose = FALSE
#' )
#' coef(result)
#' @export
fmvmr_twosample <- function(G_exposure,
                            fpca_results,
                            by_outcome,
                            sy_outcome,
                            ny_outcome,
                            max_nPC = NA,
                            true_effects = NULL,
                            verbose = TRUE) {

  if (!is.list(fpca_results) || length(fpca_results) < 1) {
    stop("fpca_results must be a list of m FPCA objects")
  }

  m <- length(fpca_results)

  if (length(by_outcome) != ncol(G_exposure)) {
    stop("by_outcome length must equal number of instruments")
  }

  if (length(sy_outcome) != length(by_outcome)) {
    stop("sy_outcome must have same length as by_outcome")
  }

  XYmodels <- recycle_arg(true_effects, m, default = NA)

  if (verbose) {
    cat("\n=== Two-Sample MV-FMR ===\n")
    cat("Exposures:", m, "\n")
    cat("Exposure N:", nrow(G_exposure), "\n")
    cat("Outcome N:", ny_outcome, "\n\n")
  }

  result <- AUTOMATIC_Multi_FMVMR_twosample_simple(
    Gmatrix = G_exposure,
    res_list = fpca_results,
    by_used = by_outcome,
    sy_used = sy_outcome,
    ny_used = ny_outcome,
    max_nPC = max_nPC,
    XYmodels = XYmodels,
    basis = "eigenfunction"
  )

  out <- structure(
    list(
      coefficients = result$MPCMRest,
      vcov = result$MPCMRvar,
      effects = lapply(result$ggdata, function(g) g$effect),
      time_grid = result$ggdata[[1]]$time,
      confidence_intervals = list(
        lower = lapply(result$ggdata, function(g) g$effect_low),
        upper = lapply(result$ggdata, function(g) g$effect_up)
      ),
      nPC_used = result$nPC_used,
      offsets = result$offsets,
      Q_stat = result$Q_stat,
      Q_pval = result$Q_pval,
      performance = if (!is.null(result$MISE) && any(!sapply(result$MISE, is.null))) {
        list(MISE = result$MISE, Coverage = result$Coverage_rate)
      } else NULL,
      plots = list(effects = result$p),
      raw_result = result,
      design = "twosample",
      n_exposure = nrow(G_exposure),
      n_outcome = ny_outcome,
      n_exposures = m
    ),
    class = c("fmvmr_twosample", "fmvmr")
  )

  if (verbose) cat("=== Complete ===\n\n")

  return(out)
}


#' Two-Sample Separate Univariable Functional MR
#'
#' Separate estimation for each exposure using outcome GWAS summary statistics.
#'
#' @param G_list List of length m of genetic instrument matrices, one per exposure (N x J_k)
#' @param fpca_results List of length m of FPCA objects, same length as G_list
#' @param by_outcome_list List of length m of SNP-outcome beta vectors, one per exposure
#' @param sy_outcome_list List of length m of SNP-outcome standard error vectors, one per exposure
#' @param ny_outcome Outcome GWAS sample size
#' @param max_nPC Maximum number of principal components to retain per exposure (length 1 or m; NA = automatically determined)
#' @param true_effects Length-m vector of true effect model codes, one per exposure (simulation only)
#' @param verbose Print progress messages and diagnostics during computation
#'
#' @return fmvmr_separate_twosample object
#' @examples
#' set.seed(1)
#' sim_data <- getX_multi_exposure(N = 60, J = 8, nSparse = 5, n_exposures = 2)
#' outcome_data <- getY_multi_exposure(sim_data, XYmodels = c("2", "8"))
#' fpca_results <- lapply(sim_data$exposures, function(exp_k) {
#'   fdapace::FPCA(exp_k$Ly_sim, exp_k$Lt_sim,
#'                 list(dataType = "Sparse", error = TRUE, verbose = FALSE))
#' })
#' # Simulate outcome GWAS summary statistics for the two-sample design
#' by_outcome <- sapply(1:8, function(j) {
#'   coef(lm(outcome_data$Y ~ sim_data$details$G[, j]))[2]
#' })
#' sy_outcome <- sapply(1:8, function(j) {
#'   summary(lm(outcome_data$Y ~ sim_data$details$G[, j]))$coefficients[2, 2]
#' })
#' result <- fmvmr_separate_twosample(
#'   G_list = list(sim_data$details$G, sim_data$details$G),
#'   fpca_results = fpca_results,
#'   by_outcome_list = list(by_outcome, by_outcome),
#'   sy_outcome_list = list(sy_outcome, sy_outcome),
#'   ny_outcome = 60,
#'   max_nPC = c(2, 2),
#'   verbose = FALSE
#' )
#' result$exposures[[1]]$coefficients
#' @export
fmvmr_separate_twosample <- function(G_list,
                                     fpca_results,
                                     by_outcome_list,
                                     sy_outcome_list,
                                     ny_outcome,
                                     max_nPC = NA,
                                     true_effects = NULL,
                                     verbose = TRUE) {

  if (!is.list(G_list) || length(G_list) < 1) {
    stop("G_list must be a list of m genetic instrument matrices")
  }

  m <- length(G_list)

  if (!is.list(fpca_results) || length(fpca_results) != m) {
    stop("fpca_results must be a list with the same length as G_list (", m, ")")
  }

  if (!is.list(by_outcome_list) || length(by_outcome_list) != m) {
    stop("by_outcome_list must be a list with the same length as G_list (", m, ")")
  }

  if (!is.list(sy_outcome_list) || length(sy_outcome_list) != m) {
    stop("sy_outcome_list must be a list with the same length as G_list (", m, ")")
  }

  if (length(by_outcome_list[[1]]) != ncol(G_list[[1]])) {
    stop("by_outcome_list[[1]] length must equal number of instruments in G_list[[1]]")
  }

  if (length(sy_outcome_list[[1]]) != length(by_outcome_list[[1]])) {
    stop("sy_outcome_list[[1]] must have same length as by_outcome_list[[1]]")
  }

  XYmodels <- recycle_arg(true_effects, m, default = NA)

  if (verbose) {
    cat("\n=== Two-Sample Separate U-FMR ===\n")
    cat("Exposures:", m, "\n")
    cat("Exposure N:", nrow(G_list[[1]]), "\n")
    cat("Outcome N:", ny_outcome, "\n\n")
  }

  result <- Separate_Multi_FMVMR_twosample_simple(
    Gmatrix_list = G_list,
    res_list = fpca_results,
    by_used_list = by_outcome_list,
    sy_used_list = sy_outcome_list,
    ny_used = ny_outcome,
    max_nPC = max_nPC,
    XYmodels = XYmodels,
    basis = "eigenfunction"
  )

  exposures <- lapply(seq_len(m), function(k) {
    list(
      coefficients = result$MPCMRest[[k]],
      vcov = result$MPCMRvar[[k]],
      effect = result$ggdata[[k]]$effect,
      nPC_used = result$nPC_used[[k]],
      performance = if (!is.null(result$MISE[[k]])) {
        list(MSE = result$MISE[[k]], Coverage = result$Coverage_rate[[k]])
      } else NULL
    )
  })

  out <- structure(
    list(
      exposures = exposures,
      plots = list(effects = result$p),
      raw_result = result,
      design = "twosample",
      n_exposure = nrow(G_list[[1]]),
      n_outcome = ny_outcome,
      n_exposures = m
    ),
    class = c("fmvmr_separate_twosample", "fmvmr_separate")
  )

  if (verbose) cat("=== Complete ===\n\n")

  return(out)
}
