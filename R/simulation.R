# ============= SIMULATION FUNCTIONS # =============
# Functions to generate genetic instruments, exposures, and outcomes
# for simulation studies and testing

#' Generate multi-exposure data with genetic instruments
#'
#' @param N Sample size
#' @param J Number of genetic instruments (per exposure, if separate_G = TRUE)
#' @param ZXmodel Model type (currently not used)
#' @param nSparse Number of sparse observations per subject
#' @param NT Number of points
#' @param TT Max observation period
#' @param n_exposures Number of exposures to simulate (m)
#' @param shared_effect Whether all exposures share the same time-varying confounding
#' @param separate_G Whether to use separate instruments for each exposure
#' @param shared_G_proportion Proportion of shared instruments (0-1)
#'
#' @return List with per-exposure sparse data and genetic instruments
#' @examples
#' set.seed(1)
#' sim_data <- getX_multi_exposure(N = 50, J = 8, nSparse = 5, n_exposures = 2)
#' length(sim_data$exposures)
#' dim(sim_data$details$G)
#' @export
getX_multi_exposure <- function(N = 10000,
                                J = 30,
                                ZXmodel = 'A',
                                nSparse = 10,
                                NT = 1000,
                                TT = 50,
                                n_exposures = 2,
                                shared_effect = TRUE,
                                separate_G = FALSE,
                                shared_G_proportion = 0.15) {

  m <- n_exposures

  times <- seq(0, TT, len = (NT + 1))
  Times <- times[-1]
  times_squ <- times^2

  # Generate genetic instruments
  if (separate_G) {
    J_vec <- rep(J, m)

    if (shared_G_proportion > 0 && shared_G_proportion <= 1) {
      # Shared + unique instruments
      n_shared <- floor(J * shared_G_proportion)

      G_shared <- matrix(stats::rbinom(n_shared * N, 2, 0.3), N, n_shared)
      G_unique_list <- lapply(seq_len(m), function(k) {
        matrix(stats::rbinom((J_vec[k] - n_shared) * N, 2, 0.3), N, J_vec[k] - n_shared)
      })

      G_list <- lapply(seq_len(m), function(k) cbind(G_shared, G_unique_list[[k]]))
      G <- do.call(cbind, c(list(G_shared), G_unique_list))
    } else {
      # Completely independent instruments
      G_list <- lapply(seq_len(m), function(k) matrix(stats::rbinom(J_vec[k] * N, 2, 0.3), N, J_vec[k]))
      G <- do.call(cbind, G_list)
    }
  } else {
    # Single G matrix shared by all exposures
    G <- matrix(stats::rbinom(J * N, 2, 0.3), N, J)
    G_list <- lapply(seq_len(m), function(k) G)
  }

  # Generate baseline confounder-effect curves
  MGX_list <- lapply(seq_len(m), function(k) {
    a_k <- stats::runif(J, -0.1, 0.1)
    b_k <- stats::runif(J, -0.004, 0.004)
    c_k <- stats::runif(J, -0, 0)
    t(t(matrix(rep(times_squ[-1], J), NT, J) %*% diag(c_k)) +
        t(matrix(rep(times[-1], J), NT, J) %*% diag(b_k)) +
        a_k)
  })

  # Generate time-varying confounding and exposures
  if (shared_effect) {
    UUU <- stats::rnorm(N, 0, 1) + t(apply(matrix(stats::rnorm(NT * N, 0, sqrt(1 / NT)), N, NT), 1, cumsum))
    EX <- apply(matrix(stats::rnorm(NT * N, 0, sqrt(1 / NT)), N, NT), 1, cumsum)

    X_list <- lapply(seq_len(m), function(k) G_list[[k]] %*% t(MGX_list[[k]]) + UUU + t(EX))
  } else {
    X_list <- vector("list", m)
    for (k in seq_len(m)) {
      UUU_k <- stats::rnorm(N, 0, 1) + t(apply(matrix(stats::rnorm(NT * N, 0, sqrt(1 / NT)), N, NT), 1, cumsum))
      EX_k <- apply(matrix(stats::rnorm(NT * N, 0, sqrt(1 / NT)), N, NT), 1, cumsum)
      X_list[[k]] <- G_list[[k]] %*% t(MGX_list[[k]]) + UUU_k + t(EX_k)
    }
  }

  # Generate sparse observations
  Ly_sim_list <- vector("list", m)
  Lt_sim_list <- vector("list", m)

  for (k in seq_len(m)) {
    Ly_sim_k <- list()
    Lt_sim_k <- list()

    for (i in 1:nrow(X_list[[k]])) {
      index_sparse <- (1:length(Times))[sort(sample(1:length(Times), nSparse))]
      time_sparse <- Times[index_sparse]
      Ly_sim_k[[i]] <- X_list[[k]][i, index_sparse]
      Lt_sim_k[[i]] <- time_sparse
    }

    Ly_sim_list[[k]] <- Ly_sim_k
    Lt_sim_list[[k]] <- Lt_sim_k
  }

  # Create output structure
  if (separate_G) {
    G_names_list <- lapply(seq_len(m), function(k) paste0('G', k, '_', 1:ncol(G_list[[k]])))

    DAT_list <- lapply(seq_len(m), function(k) {
      DAT_k <- as.data.frame(cbind(G_list[[k]], X_list[[k]][, (1:TT) * NT / TT]))
      names(DAT_k) <- c(G_names_list[[k]], paste0('X', k, 1:TT))
      DAT_k
    })

    DAT <- as.data.frame(do.call(cbind, c(G_list, lapply(X_list, function(X) X[, (1:TT) * NT / TT]))))
    names(DAT) <- c(unlist(G_names_list), unlist(lapply(seq_len(m), function(k) paste0('X', k, 1:TT))))
  } else {
    G_names <- paste0('G', 1:ncol(G))

    DAT <- as.data.frame(cbind(G, do.call(cbind, lapply(X_list, function(X) X[, (1:TT) * NT / TT]))))
    names(DAT) <- c(G_names, unlist(lapply(seq_len(m), function(k) paste0('X', k, 1:TT))))

    G_list <- vector("list", m)
  }

  exposures <- vector("list", m)
  for (k in seq_len(m)) {
    if (separate_G) {
      exposures[[k]] <- list(DAT = DAT_list[[k]], Ly_sim = Ly_sim_list[[k]], Lt_sim = Lt_sim_list[[k]])
    } else {
      cols <- c(1:J, (J + (k - 1) * TT + 1):(J + k * TT))
      exposures[[k]] <- list(DAT = DAT[, cols], Ly_sim = Ly_sim_list[[k]], Lt_sim = Lt_sim_list[[k]])
    }
  }

  RES <- list(
    exposures = exposures,
    details = list(
      J = J,
      TT = TT,
      NT = NT,
      G = G,
      G_list = G_list,
      times = times,
      X_list = X_list,
      n_exposures = m,
      shared_effect = shared_effect,
      separate_G = separate_G,
      shared_G_proportion = shared_G_proportion
    ),
    DAT = DAT
  )

  return(RES)
}

#' Generate multi-exposure mediation data with genetic instruments
#'
#' @param N Sample size
#' @param J Number of genetic instruments per exposure
#' @param ZXmodel Model type (currently not used, kept for compatibility)
#' @param nSparse Number of sparse observations per subject
#' @param n_exposures Number of exposures to simulate (m)
#' @param mediation_strength m x m numeric matrix of pairwise mediation strengths: entry
#'   \code{[j, k]} (with j < k) is the strength with which exposure j mediates its effect
#'   onto exposure k, generated later in the sequence. Must be strictly upper triangular
#'   (entries with j >= k must be 0). Default: NULL, i.e. no mediation (all-zero matrix).
#' @param separate_G Whether to use separate instruments for each exposure
#' @param shared_G_proportion Proportion of shared instruments (0-1)
#' @param mediation_type Character, or m x m character matrix mirroring mediation_strength:
#'   type of mediation effect for each pair, one of "linear" (default), "nonlinear", or
#'   "time_varying".
#'
#' @return List with same structure as getX_multi_exposure()
#' @examples
#' set.seed(1)
#' # Exposure 1 mediates onto exposure 2 with strength 0.3
#' mediation_strength <- matrix(c(0, 0, 0.3, 0), 2, 2)
#' sim_data <- getX_multi_exposure_mediation(
#'   N = 50, J = 8, nSparse = 5, n_exposures = 2,
#'   mediation_strength = mediation_strength
#' )
#' length(sim_data$exposures)
#' @export
getX_multi_exposure_mediation <- function(N = 10000,
                                          J = 30,
                                          ZXmodel = 'A',
                                          nSparse = 10,
                                          n_exposures = 2,
                                          mediation_strength = NULL,
                                          separate_G = FALSE,
                                          shared_G_proportion = 0,
                                          mediation_type = "linear") {

  m <- n_exposures
  NT <- 1000
  TT <- 50

  if (is.null(mediation_strength)) {
    mediation_strength <- matrix(0, m, m)
  }

  if (!is.matrix(mediation_strength) || !all(dim(mediation_strength) == m)) {
    stop("mediation_strength must be an m x m matrix")
  }

  invalid_mask <- lower.tri(mediation_strength, diag = TRUE)
  if (any(mediation_strength[invalid_mask] != 0)) {
    bad_idx <- which(invalid_mask & (mediation_strength != 0), arr.ind = TRUE)
    stop("mediation_strength[j, k] may only be nonzero for j < k (an earlier exposure ",
         "mediating a later one); invalid entries at: ",
         paste(apply(bad_idx, 1, function(r) paste0("[", r[1], ",", r[2], "]")), collapse = ", "))
  }

  if (is.matrix(mediation_type) && !all(dim(mediation_type) == m)) {
    stop("mediation_type matrix must be m x m")
  }

  get_mediation_type <- function(j, k) {
    if (is.matrix(mediation_type)) mediation_type[j, k] else mediation_type
  }

  times <- seq(0, TT, len = (NT + 1))
  Times <- times[-1]
  times_squ <- times^2

  ## Genetic instruments
  if (separate_G) {
    J_vec <- rep(J, m)

    if (shared_G_proportion > 0 && shared_G_proportion <= 1) {
      n_shared <- floor(J * shared_G_proportion)

      G_shared <- matrix(stats::rbinom(n_shared * N, 2, 0.3), N, n_shared)
      G_unique_list <- lapply(seq_len(m), function(k) {
        matrix(stats::rbinom((J_vec[k] - n_shared) * N, 2, 0.3), N, J_vec[k] - n_shared)
      })

      G_list <- lapply(seq_len(m), function(k) cbind(G_shared, G_unique_list[[k]]))
      G <- do.call(cbind, c(list(G_shared), G_unique_list))
    } else {
      G_list <- lapply(seq_len(m), function(k) matrix(stats::rbinom(J_vec[k] * N, 2, 0.3), N, J_vec[k]))
      G <- do.call(cbind, G_list)
    }
  } else {
    G <- matrix(stats::rbinom(J * N, 2, 0.3), N, J)
    G_list <- lapply(seq_len(m), function(k) G)
  }

  ## Baseline genetic effects
  MGX_list <- lapply(seq_len(m), function(k) {
    noise_k <- stats::rnorm(NT * J, 0, 0.05)
    a_k <- stats::runif(J, -0.1, 0.1)
    b_k <- stats::runif(J, -0.004, 0.004)
    c_k <- stats::runif(J, 0, 0)

    t(t(matrix(rep(times_squ[-1], J), NT, J) %*% diag(c_k)) +
        t(matrix(rep(times[-1], J), NT, J) %*% diag(b_k)) +
        a_k) + noise_k
  })

  ## Time-varying confounding (shared, as in mediation)
  UUU <- stats::rnorm(N, 0, 1) +
    t(apply(matrix(stats::rnorm(NT * N, 0, sqrt(1 / NT)), N, NT), 1, cumsum))

  EX <- apply(matrix(stats::rnorm(NT * N, 0, sqrt(1 / NT)), N, NT), 1, cumsum)

  ## Generate exposures sequentially so each can be mediated by any earlier exposure
  X_list <- vector("list", m)
  mediation_effect_list <- vector("list", m)

  for (k in seq_len(m)) {
    mediation_effect_k <- matrix(0, N, NT)

    for (j in seq_len(k - 1)) {
      s <- mediation_strength[j, k]
      if (s != 0) {
        Xj <- X_list[[j]]
        type_jk <- get_mediation_type(j, k)

        if (type_jk == "linear") {
          mediation_effect_k <- mediation_effect_k + s * Xj
        } else if (type_jk == "nonlinear") {
          Xjs <- scale(Xj)
          mediation_effect_k <- mediation_effect_k + s * (Xjs + 0.1 * Xjs^2)
        } else if (type_jk == "time_varying") {
          w <- sin(times[-1] * pi / max(times[-1])) + 1
          w <- w / mean(w)
          mediation_effect_k <- mediation_effect_k + s * Xj * matrix(rep(w, N), N, NT, byrow = TRUE)
        } else {
          stop("Invalid mediation_type")
        }
      }
    }

    mediation_effect_list[[k]] <- mediation_effect_k
    X_list[[k]] <- G_list[[k]] %*% t(MGX_list[[k]]) + mediation_effect_k + UUU + t(EX)
  }

  ## Sparse observations
  Ly_sim_list <- vector("list", m)
  Lt_sim_list <- vector("list", m)

  for (k in seq_len(m)) {
    Ly_sim_k <- Lt_sim_k <- vector("list", N)

    for (i in 1:N) {
      idx <- sort(sample(seq_along(Times), nSparse))
      Ly_sim_k[[i]] <- X_list[[k]][i, idx]
      Lt_sim_k[[i]] <- Times[idx]
    }

    Ly_sim_list[[k]] <- Ly_sim_k
    Lt_sim_list[[k]] <- Lt_sim_k
  }

  ## Data frames
  if (separate_G) {
    G_names_list <- lapply(seq_len(m), function(k) paste0("G", k, "_", 1:ncol(G_list[[k]])))

    DAT_list <- lapply(seq_len(m), function(k) {
      DAT_k <- as.data.frame(cbind(G_list[[k]], X_list[[k]][, (1:TT) * NT / TT]))
      names(DAT_k) <- c(G_names_list[[k]], paste0("X", k, 1:TT))
      DAT_k
    })

    DAT <- as.data.frame(do.call(cbind, c(G_list, lapply(X_list, function(X) X[, (1:TT) * NT / TT]))))
    names(DAT) <- c(unlist(G_names_list), unlist(lapply(seq_len(m), function(k) paste0("X", k, 1:TT))))
  } else {
    G_names <- paste0("G", 1:ncol(G))
    DAT <- as.data.frame(cbind(G, do.call(cbind, lapply(X_list, function(X) X[, (1:TT) * NT / TT]))))
    names(DAT) <- c(G_names, unlist(lapply(seq_len(m), function(k) paste0("X", k, 1:TT))))
    G_list <- vector("list", m)
  }

  ## Return object
  exposures <- vector("list", m)
  for (k in seq_len(m)) {
    if (separate_G) {
      exposures[[k]] <- list(DAT = DAT_list[[k]], Ly_sim = Ly_sim_list[[k]], Lt_sim = Lt_sim_list[[k]])
    } else {
      cols <- c(1:J, (J + (k - 1) * TT + 1):(J + k * TT))
      exposures[[k]] <- list(DAT = DAT[, cols], Ly_sim = Ly_sim_list[[k]], Lt_sim = Lt_sim_list[[k]])
    }
  }

  RES <- list(
    exposures = exposures,
    details = list(
      J = J,
      TT = TT,
      NT = NT,
      G = G,
      G_list = G_list,
      times = times,
      X_list = X_list,
      n_exposures = m,
      mediation_strength = mediation_strength,
      mediation_type = mediation_type,
      mediation_effect = mediation_effect_list,
      separate_G = separate_G,
      shared_G_proportion = shared_G_proportion
    ),
    DAT = DAT
  )

  return(RES)
}


#' Generate outcome from exposures
#'
#' @param RES Output from getX_multi_exposure() or getX_multi_exposure_mediation()
#' @param XYmodels Length-m vector of effect models per exposure, one of '0'-'9' (default: '1' for all)
#' @param X_effects Length-m logical vector: include each exposure's effect? (default: TRUE for all)
#' @param outcome_type "continuous" or "binary"
#'
#' @return Data frame with outcome Y
#' @examples
#' set.seed(1)
#' sim_data <- getX_multi_exposure(N = 50, J = 8, nSparse = 5, n_exposures = 2)
#' dat <- getY_multi_exposure(sim_data, XYmodels = c("2", "8"), outcome_type = "continuous")
#' head(dat$Y)
#' @export
getY_multi_exposure <- function(RES,
                                XYmodels = NULL,
                                X_effects = NULL,
                                outcome_type = "continuous") {

  m <- RES$details$n_exposures

  XYmodels <- recycle_arg(XYmodels, m, default = '1')
  X_effects <- recycle_arg(X_effects, m, default = TRUE)

  # Effect functions
  get_fun <- function(model_type) {
    switch(model_type,
           '0' = function(t) 0 * (t < Inf),
           '1' = function(t) 0.1 * (t < Inf),
           '2' = function(t) 0.02 * t,
           '3' = function(t) 0.5 - 0.02 * t,
           '4' = function(t) 0.1 * (t < 20),
           '5' = function(t) 0.1 * (t > 30),
           '6' = function(t) 0.05 * (-t + 20) * (t < 20),
           '7' = function(t) 0.05 * (t - 30) * (t > 30),
           '8' = function(t) 0.002 * t^2 - 0.11 * t + 0.5,
           '9' = function(t) -0.00002 * t^3 + 0.004 * t^2 - 0.2 * t + 1
    )
  }

  X_list <- RES$details$X_list
  TT <- RES$details$TT
  NT <- RES$details$NT

  times <- seq(0, TT, len = (NT + 1))

  effect_vec_list <- lapply(seq_len(m), function(k) {
    model_k <- if (is.na(XYmodels[k])) '0' else as.character(XYmodels[k])
    fun_k <- get_fun(model_k)
    effect_vec_k <- as.numeric(fun_k(times[-1]))
    if (!X_effects[k]) effect_vec_k <- rep(0, length(effect_vec_k))
    effect_vec_k
  })

  linear_predictor <- Reduce(`+`, lapply(seq_len(m), function(k) {
    (X_list[[k]] %*% effect_vec_list[[k]]) * TT / NT
  }))

  if (outcome_type == "continuous") {
    Y <- linear_predictor + stats::rnorm(nrow(X_list[[1]]), 0, 1)
  } else if (outcome_type == "binary") {
    prob_Y <- 1 / (1 + exp(-linear_predictor))
    Y <- stats::rbinom(nrow(X_list[[1]]), 1, prob_Y)
  } else {
    stop("outcome_type must be either 'continuous' or 'binary'")
  }

  # Get DAT from RES
  if (!is.null(RES$DAT)) {
    DAT <- RES$DAT
  } else {
    J <- ncol(RES$details$G)
    DAT <- data.frame(
      RES$details$G,
      do.call(cbind, lapply(X_list, function(X) X[, (1:TT) * NT / TT]))
    )
    names(DAT) <- c(paste0('G', 1:J), unlist(lapply(seq_len(m), function(k) paste0('X', k, 1:TT))))
  }

  DAT$Y <- as.numeric(Y)

  attr(DAT, "model_info") <- list(
    XYmodels = XYmodels,
    X_effects = X_effects,
    outcome_type = outcome_type,
    effect_vec_list = effect_vec_list,
    times = times[-1]
  )

  return(DAT)
}
