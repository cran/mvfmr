# ============= MVFMR - MAIN FUNCTIONS # =============

#' Joint Multivariable Functional Mendelian Randomization
#'
#' @param G Genetic instrument matrix (N x J)
#' @param fpca_results List with two FPCA objects from fdapace (res1 and res2)
#' @param Y Outcome vector
#' @param outcome_type Type of outcome: "continuous" for numeric outcomes, "binary" for 0/1 outcomes
#' @param method Estimation method: "gmm" (Generalized Method of Moments), "cf" (control function), or "cf-lasso" (control function with Lasso)
#' @param nPC1 Fixed number of principal components to retain for exposure 1 (NA = select automatically)
#' @param max_nPC1 Maximum number of principal components to retain for exposure 1 (NA = automatically determined)
#' @param nPC2 Fixed number of principal components to retain for exposure 2 (NA = select automatically)
#' @param max_nPC2 Maximum number of principal components to retain for exposure 2 (NA = automatically determined)
#' @param improvement_threshold Minimum cross-validation improvement required to add an additional principal component
#' @param bootstrap Whether to compute confidence intervals using bootstrap resampling
#' @param n_bootstrap Number of bootstrap replicates (only used if bootstrap = TRUE)
#' @param n_cores Number of CPU cores to use for parallel computations
#' @param true_effects List with true_effect1 and true_effect2 (simulation only)
#' @param X_true List with X1_true and X2_true curves (simulation only)
#' @param verbose Print progress and diagnostic messages during computation
#'
#' @return mvfmr object with:
#' \itemize{
#'   \item coefficients - Estimated beta coefficients
#'   \item vcov - Variance-covariance matrix
#'   \item effects - List with effect1 and effect2 curves
#'   \item nPC_used - Components selected (nPC1, nPC2)
#'   \item diagnostics - F-statistics, instrument diagnostics
#'   \item performance - MISE, coverage (if true effects provided)
#' }
#' @export 
mvfmr <- function(G,
                  fpca_results,
                  Y,
                  outcome_type = c("continuous", "binary"),
                  method = c("gmm", "cf", "cf-lasso"),
                  nPC1 = NA,
                  max_nPC1 = NA,
                  nPC2 = NA,
                  max_nPC2 = NA,
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
  if (!is.list(fpca_results) || length(fpca_results) != 2) {
    stop("fpca_results must be a list with 2 FPCA objects (res1, res2)")
  }
  
  res1 <- fpca_results[[1]]
  res2 <- fpca_results[[2]]
  
  # Extract true effects if provided
  X1Ymodel <- NA
  X2Ymodel <- NA
  X1_true_mat <- matrix()
  X2_true_mat <- matrix()
  
  if (!is.null(true_effects)) {
    X1Ymodel <- true_effects$model1
    X2Ymodel <- true_effects$model2
  }
  
  if (!is.null(X_true)) {
    X1_true_mat <- X_true$X1_true
    X2_true_mat <- X_true$X2_true
  }
  
  N <- nrow(G)
  IDmatch <- 1:N
  
  if (verbose) {
    cat("\n=== Functional Multivariable MR ===\n")
    cat("Sample size:", N, "\n")
    cat("Outcome:", outcome_type, "\n")
    cat("Method:", method, "\n\n")
  }
  
  # Call original AUTOMATIC_Multi_MVFMR function
  result <- AUTOMATIC_Multi_MVFMR(
    Gmatrix = G,
    res1 = res1,
    res2 = res2,
    Yvector = Y,
    IDmatch = IDmatch,
    nPC1_selected = NA,
    max_nPC1 = max_nPC1,
    nPC2_selected = NA,
    max_nPC2 = max_nPC2,
    X1_true = X1_true_mat,
    X2_true = X2_true_mat,
    method = method,
    basis = "eigenfunction",
    outcome = outcome_type,
    bootstrap = bootstrap,
    n_B = n_bootstrap,
    improvement_threshold = improvement_threshold,
    X1Ymodel = X1Ymodel,
    X2Ymodel = X2Ymodel,
    num_cores_set = n_cores,
    verbose = verbose
  )
  
  # Reformat output to package standard
  out <- structure(
    list(
      coefficients = result$MPCMRest,
      vcov = result$MPCMRvar,
      effects = list(
        effect1 = if (!is.null(result$ggdata)) result$ggdata$effect1 else NULL,
        effect2 = if (!is.null(result$ggdata)) result$ggdata$effect2 else NULL,
        time_grid = if (!is.null(result$ggdata)) result$ggdata$time else NULL
      ),
      confidence_intervals = list(
        effect1_lower = if (!is.null(result$ggdata)) result$ggdata$effect1_low else NULL,
        effect1_upper = if (!is.null(result$ggdata)) result$ggdata$effect1_up else NULL,
        effect2_lower = if (!is.null(result$ggdata)) result$ggdata$effect2_low else NULL,
        effect2_upper = if (!is.null(result$ggdata)) result$ggdata$effect2_up else NULL
      ),
      nPC_used = list(nPC1 = result$nPC1_used, nPC2 = result$nPC2_used),
      performance = if (!is.null(X1Ymodel)) {
        list(
          MISE1 = result$MISE1,
          MISE2 = result$MISE2,
          Coverage1 = result$Coverage_rate1,
          Coverage2 = result$Coverage_rate2
        )
      } else NULL,
      plots = list(p1 = result$p1, p2 = result$p2, plot_beta = result$plot_beta),
      raw_result = result,
      n_exposures = 2,
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
#' @param G1  Genetic instrument matrix for exposure 1
#' @param G2  Genetic instrument matrix for exposure 2, or NULL if only a single exposure is analyzed
#' @param fpca_results List of FPCA objects
#' @param Y Outcome vector
#' @param outcome_type Type of outcome: "continuous" for numeric outcomes, "binary" for 0/1 outcomes
#' @param method Estimation method: "gmm" (Generalized Method of Moments), "cf" (control function), or "cf-lasso" (control function with Lasso)
#' @param nPC1 Fixed number of principal components to retain for exposure 1 (NA = select automatically)
#' @param max_nPC1 Maximum number of principal components to retain for exposure 1 (NA = automatically determined)
#' @param nPC2 Fixed number of principal components to retain for exposure 2 (NA = select automatically)
#' @param max_nPC2 Maximum number of principal components to retain for exposure 2 (NA = automatically determined)
#' @param improvement_threshold Minimum cross-validation improvement required to add an additional principal component
#' @param bootstrap Whether to compute confidence intervals using bootstrap resampling
#' @param n_bootstrap Number of bootstrap replicates (only used if bootstrap = TRUE)
#' @param n_cores Number of CPU cores to use for parallel computations
#' @param true_effects List with true_effect1 and true_effect2 (simulation only)
#' @param X_true List with X1_true and X2_true curves (simulation only)
#' @param verbose Print progress and diagnostic messages during computation
#'
#' @return fmvmr_separate object
#' @export 
mvfmr_separate <- function(G1,
                           G2,
                           fpca_results,
                           Y,
                           outcome_type = c("continuous", "binary"),
                           method = c("gmm", "cf", "cf-lasso"),
                           nPC1 = NA,
                           max_nPC1 = NA,
                           nPC2 = NA,
                           max_nPC2 = NA,
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
  if (!is.list(fpca_results) || (!is.null(G2) & length(fpca_results) != 2)) {
    stop("fpca_results must be a list with 2 FPCA objects")
  }
  # Extract FPCA objects
  if (!is.list(fpca_results) || (is.null(G2) & length(fpca_results) != 1)) {
    stop("fpca_results must be a list with 1 FPCA objects")
  }
  
  
  # Handle G: single matrix or list
  if (!is.null(G1)) {
    separate_G <- TRUE
    Gmatrix1 <- G1
    res1 <- fpca_results[[1]]
  } else{
    Gmatrix1 <- NULL
    res1 <- NULL
  }
  if (!is.null(G2)) {
    separate_G <- TRUE
    Gmatrix2 <- G2
    res2 <- fpca_results[[2]]
  }else{
    Gmatrix2 <- NULL
    res2 <- NULL
  }
  
  # Extract true effects if provided
  X1Ymodel <- NA
  X2Ymodel <- NA
  X1_true_mat <- matrix()
  X2_true_mat <- matrix()
  
  if (!is.null(true_effects)) {
    X1Ymodel <- true_effects$model1
    X2Ymodel <- true_effects$model2
  }
  
  N <- nrow(Gmatrix1)
  IDmatch <- 1:N
  
  if (verbose) {
    cat("\n=== Separate Univariable MR ===\n")
    cat("Separate instruments:", separate_G, "\n")
    cat("Sample size:", N, "\n\n")
  }
  
  # Call original Separate_Multi_MVFMR function
  result <- Separate_Multi_MVFMR(
    Gmatrix1 = Gmatrix1,
    Gmatrix2 = Gmatrix2,
    res1 = res1,
    res2 = res2,
    nPC1_selected = nPC1,
    max_nPC1 = max_nPC1,
    nPC2_selected = nPC2,
    max_nPC2 = max_nPC2,
    Yvector = Y,
    IDmatch = IDmatch,
    method = method,
    basis = "eigenfunction",
    outcome = outcome_type,
    bootstrap = bootstrap,
    n_B = n_bootstrap,
    improvement_threshold = improvement_threshold,
    X1Ymodel = X1Ymodel,
    X2Ymodel = X2Ymodel,
    num_cores_set = n_cores,
    verbose = verbose
  )
  
  # Reformat output
  out <- structure(
    list(
      exposure1 = list(
        coefficients = result$MPCMRest1,
        vcov = result$MPCMRvar1,
        effect = if (!is.null(result$ggdata1)) result$ggdata1$effect else NULL,
        nPC_used = result$nPC_used1,
        performance = if (!is.null(X1Ymodel)) {
          list(MISE = result$MISE1, Coverage = result$Coverage_rate1)
        } else NULL
      ),
      exposure2 = list(
        coefficients = result$MPCMRest2,
        vcov = result$MPCMRvar2,
        effect = if (!is.null(result$ggdata2)) result$ggdata2$effect else NULL,
        nPC_used = result$nPC_used2,
        performance = if (!is.null(X2Ymodel)) {
          list(MISE = result$MISE2, Coverage = result$Coverage_rate2)
        } else NULL
      ),
      plots = list(p1 = result$p1, p2 = result$p2),
      raw_result = result,
      n_exposures = sum(!is.null(G1), !is.null(G2)),
      separate_instruments = separate_G,
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
#' @param fpca_results List of 2 FPCA objects 
#' @param by_outcome Vector of SNP-outcome effect estimates (betas) from the outcome GWAS, length J
#' @param sy_outcome VVector of standard errors for SNP-outcome effects, length J
#' @param ny_outcome Sample size of the outcome GWAS
#' @param max_nPC1 Maximum number of principal components to retain for exposure 1 (NA = automatically determined)
#' @param max_nPC2 Maximum number of principal components to retain for exposure 2 (NA = automatically determined)
#' @param true_effects List containing true effects for exposure 1 and exposure 2 (simulation only)
#' @param verbose Print progress messages and diagnostics during computation
#'
#' @return fmvmr_twosample object
#' @export 
fmvmr_twosample <- function(G_exposure,
                            fpca_results,
                            by_outcome,
                            sy_outcome,
                            ny_outcome,
                            max_nPC1 = NA,
                            max_nPC2 = NA,
                            true_effects = NULL,
                            verbose = TRUE) {
  
  if (!is.list(fpca_results) || length(fpca_results) != 2) {
    stop("fpca_results must be a list with 2 FPCA objects")
  }
  
  if (length(by_outcome) != ncol(G_exposure)) {
    stop("by_outcome length must equal number of instruments")
  }
  
  if (length(sy_outcome) != length(by_outcome)) {
    stop("sy_outcome must have same length as by_outcome")
  }
  
  res1 <- fpca_results[[1]]
  res2 <- fpca_results[[2]]
  
  X1Ymodel <- NA
  X2Ymodel <- NA
  if (!is.null(true_effects)) {
    X1Ymodel <- true_effects$model1
    X2Ymodel <- true_effects$model2
  }
  
  if (verbose) {
    cat("\n=== Two-Sample MV-FMR ===\n")
    cat("Exposure N:", nrow(G_exposure), "\n")
    cat("Outcome N:", ny_outcome, "\n\n")
  }
  
  result <- AUTOMATIC_Multi_FMVMR_twosample_simple(
    Gmatrix = G_exposure,
    res1 = res1,
    res2 = res2,
    by_used = by_outcome,
    sy_used = sy_outcome,
    ny_used = ny_outcome,
    max_nPC1 = max_nPC1,
    max_nPC2 = max_nPC2,
    X1Ymodel = X1Ymodel,
    X2Ymodel = X2Ymodel,
    basis = "eigenfunction"
  )
  
  out <- structure(
    list(
      coefficients = result$MPCMRest,
      vcov = result$MPCMRvar,
      effects = list(
        effect1 = if (!is.null(result$ggdata)) result$ggdata$effect1 else NULL,
        effect2 = if (!is.null(result$ggdata)) result$ggdata$effect2 else NULL,
        time_grid = if (!is.null(result$ggdata)) result$ggdata$time else NULL
      ),
      confidence_intervals = list(
        effect1_lower = if (!is.null(result$ggdata)) result$ggdata$effect1_low else NULL,
        effect1_upper = if (!is.null(result$ggdata)) result$ggdata$effect1_up else NULL,
        effect2_lower = if (!is.null(result$ggdata)) result$ggdata$effect2_low else NULL,
        effect2_upper = if (!is.null(result$ggdata)) result$ggdata$effect2_up else NULL
      ),
      nPC_used = list(nPC1 = result$nPC1_used, nPC2 = result$nPC2_used),
      Q_stat = result$Q_stat,
      Q_pval = result$Q_pval,
      performance = if (!is.null(X1Ymodel) && !is.na(X1Ymodel)) {
        list(
          MISE1 = result$MISE1,
          MISE2 = result$MISE2,
          Coverage1 = result$Coverage_rate1,
          Coverage2 = result$Coverage_rate2
        )
      } else NULL,
      plots = list(p1 = result$p1, p2 = result$p2),
      raw_result = result,
      design = "twosample",
      n_exposure = nrow(G_exposure),
      n_outcome = ny_outcome
    ),
    class = c("fmvmr_twosample", "fmvmr")
  )
  
  if (verbose) cat("=== Complete ===\n\n")
  
  return(out)
}


#' Two-Sample Separate Univariable Functional MR
#'
#' Separate estimation for each exposure using outcome GWAS summary statistics.
#' For single exposure: set G2 = NULL, by2 = NULL, sy2 = NULL.
#'
#' @param G1_exposure Genetic instrument matrix from exposure 1 (N × J1)
#' @param G2_exposure Genetic instrument matrix from exposure 2 (N × J2) or NULL for single exposure
#' @param fpca_results List of 2 FPCA objects
#' @param by_outcome1 SNP-outcome betas for exposure 1 instruments
#' @param by_outcome2 SNP-outcome betas for exposure 2 instruments or NULL
#' @param sy_outcome1 Standard errors for exposure 1
#' @param sy_outcome2 Standard errors for exposure 2 or NULL
#' @param ny_outcome Outcome GWAS sample size
#' @param max_nPC1 Maximum number of principal components to retain for exposure 1 (NA = automatically determined)
#' @param max_nPC2 Maximum number of principal components to retain for exposure 2 (NA = automatically determined)
#' @param true_effects List containing true effects for exposure 1 and exposure 2 (simulation only)
#' @param verbose Print progress messages and diagnostics during computation

#'
#' @return fmvmr_separate_twosample object
#' @export 
fmvmr_separate_twosample <- function(G1_exposure,
                                     G2_exposure = NULL,
                                     fpca_results,
                                     by_outcome1,
                                     by_outcome2 = NULL,
                                     sy_outcome1,
                                     sy_outcome2 = NULL,
                                     ny_outcome,
                                     max_nPC1 = NA,
                                     max_nPC2 = NA,
                                     true_effects = NULL,
                                     verbose = TRUE) {
  
  if (!is.list(fpca_results) || length(fpca_results) != 2) {
    stop("fpca_results must be a list with 2 FPCA objects")
  }
  
  if (length(by_outcome1) != ncol(G1_exposure)) {
    stop("by_outcome1 length must equal number of instruments in G1")
  }
  
  if (length(sy_outcome1) != length(by_outcome1)) {
    stop("sy_outcome1 must have same length as by_outcome1")
  }
  
  res1 <- fpca_results[[1]]
  res2 <- fpca_results[[2]]
  
  X1Ymodel <- NA
  X2Ymodel <- NA
  if (!is.null(true_effects)) {
    X1Ymodel <- true_effects$model1
    X2Ymodel <- true_effects$model2
  }
  
  if (verbose) {
    cat("\n=== Two-Sample Separate U-FMR ===\n")
    cat("Exposure N:", nrow(G1_exposure), "\n")
    cat("Outcome N:", ny_outcome, "\n\n")
  }
  
  result <- Separate_Multi_FMVMR_twosample_simple(
    Gmatrix1 = G1_exposure,
    Gmatrix2 = G2_exposure,
    res1 = res1,
    res2 = res2,
    by_used1 = by_outcome1,
    by_used2 = by_outcome2,
    sy_used1 = sy_outcome1,
    sy_used2 = sy_outcome2,
    ny_used = ny_outcome,
    max_nPC1 = max_nPC1,
    max_nPC2 = max_nPC2,
    X1Ymodel = X1Ymodel,
    X2Ymodel = X2Ymodel
  )
  
  out <- structure(
    list(
      exposure1 = list(
        coefficients = result$MPCMRest1,
        vcov = result$MPCMRvar1,
        effect = if (!is.null(result$ggdata1)) result$ggdata1$effect else NULL,
        nPC_used = result$nPC_used1,
        performance = if (!is.null(X1Ymodel) && !is.na(X1Ymodel)) {
          list(MSE = result$MISE1, Coverage = result$Coverage_rate1)
        } else NULL
      ),
      exposure2 = if (!is.null(G2_exposure)) {
        list(
          coefficients = result$MPCMRest2,
          vcov = result$MPCMRvar2,
          effect = if (!is.null(result$ggdata2)) result$ggdata2$effect else NULL,
          nPC_used = result$nPC_used2,
          performance = if (!is.null(X2Ymodel) && !is.na(X2Ymodel)) {
            list(MSE = result$MISE2, Coverage = result$Coverage_rate2)
          } else NULL
        )
      } else NULL,
      plots = list(p1 = result$p1, p2 = result$p2),
      raw_result = result,
      design = "twosample",
      n_exposure = nrow(G1_exposure),
      n_outcome = ny_outcome
    ),
    class = c("fmvmr_separate_twosample", "fmvmr_separate")
  )
  
  if (verbose) cat("=== Complete ===\n\n")
  
  return(out)
}