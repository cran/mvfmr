# ============= SIMULATION FUNCTIONS # =============
# Functions to generate genetic instruments, exposures, and outcomes
# for simulation studies and testing

#' Generate multi-exposure data with genetic instruments
#'
#' @param N Sample size
#' @param J Number of genetic instruments
#' @param ZXmodel Model type (currently not used)
#' @param nSparse Number of sparse observations per subject
#' @param NT Number of points
#' @param TT Max observation period
#' @param shared_effect Whether X1 and X2 share confounding
#' @param separate_G Whether to use separate instruments for each exposure
#' @param shared_G_proportion Proportion of shared instruments (0-1)
#'
#' @return List with X1, X2 sparse data and genetic instruments
#' @export
getX_multi_exposure <- function(N = 10000, 
                                J = 30,
                                ZXmodel = 'A',
                                nSparse = 10,
                                NT = 1000,
                                TT=50,
                                shared_effect = TRUE,
                                separate_G = FALSE,
                                shared_G_proportion = 0.15) {
  
  NT <- 1000  
  TT <- 50 
  
  times <- seq(0, TT, len = (NT + 1))
  Times <- times[-1] 
  times_squ <- times^2 
  
  # Initialize G, G1, G2 as NULL
  G <- NULL
  G1 <- NULL
  G2 <- NULL
  
  # Generate genetic instruments
  if (separate_G) {
    J1 <- J
    J2 <- J
    
    if (shared_G_proportion > 0 && shared_G_proportion <= 1) {
      # Shared + unique instruments
      n_shared <- floor(J * shared_G_proportion)
      
      G_shared <- matrix(stats::rbinom(n_shared * N, 2, 0.3), N, n_shared)
      G1_unique <- matrix(stats::rbinom((J1 - n_shared) * N, 2, 0.3), N, J1 - n_shared)
      G2_unique <- matrix(stats::rbinom((J2 - n_shared) * N, 2, 0.3), N, J2 - n_shared)
      
      G1 <- cbind(G_shared, G1_unique)
      G2 <- cbind(G_shared, G2_unique)
      G <- cbind(G_shared, G1_unique, G2_unique)
    } else {
      # Completely independent
      G1 <- matrix(stats::rbinom(J1 * N, 2, 0.3), N, J1)
      G2 <- matrix(stats::rbinom(J2 * N, 2, 0.3), N, J2)
      G <- cbind(G1, G2)
    }
  } else {
    # Single G matrix for both exposures
    G <- matrix(stats::rbinom(J * N, 2, 0.3), N, J)
    G1 <- G
    G2 <- G
  }
  
  # Generate baseline confounders
  if (shared_effect) {
    a1 <- stats::runif(J, -0.1, 0.1)
    b1 <- stats::runif(J, -0.004, 0.004)
    c1 <- stats::runif(J, -0, 0)
    MGX1 <- t(t(matrix(rep(times_squ[-1], J), NT, J) %*% diag(c1)) +
                t(matrix(rep(times[-1], J), NT, J) %*% diag(b1)) +
                a1)
    
    a2 <- stats::runif(J, -0.1, 0.1)
    b2 <- stats::runif(J, -0.004, 0.004)
    c2 <- stats::runif(J, -0, 0)
    MGX2 <- t(t(matrix(rep(times_squ[-1], J), NT, J) %*% diag(c2)) +
                t(matrix(rep(times[-1], J), NT, J) %*% diag(b2)) +
                a2) 
  } else {
    
    a1 <- stats::runif(J, -0.1, 0.1)
    b1 <- stats::runif(J, -0.004, 0.004)
    c1 <- stats::runif(J, -0, 0)
    MGX1 <- t(t(matrix(rep(times_squ[-1], J), NT, J) %*% diag(c1)) +
                t(matrix(rep(times[-1], J), NT, J) %*% diag(b1)) +
                a1)
    
    a2 <- stats::runif(J, -0.1, 0.1)
    b2 <- stats::runif(J, -0.004, 0.004)
    c2 <- stats::runif(J, -0, 0)
    MGX2 <- t(t(matrix(rep(times_squ[-1], J), NT, J) %*% diag(c2)) +
                t(matrix(rep(times[-1], J), NT, J) %*% diag(b2)) +
                a2)
  }
  
  # Generate time-varying confounding and exposures
  if (shared_effect) {
    UUU <- stats::rnorm(N, 0, 1) + t(apply(matrix(stats::rnorm(NT * N, 0, sqrt(1 / NT)), N, NT), 1, cumsum))
    EX <- apply(matrix(stats::rnorm(NT * N, 0, sqrt(1 / NT)), N, NT), 1, cumsum)
    
    X1 <- G1 %*% t(MGX1) + UUU + t(EX)
    X2 <- G2 %*% t(MGX2) + UUU + t(EX)
  } else {
    UUU1 <- stats::rnorm(N, 0, 1) + t(apply(matrix(stats::rnorm(NT * N, 0, sqrt(1 / NT)), N, NT), 1, cumsum))
    UUU2 <- stats::rnorm(N, 0, 1) + t(apply(matrix(stats::rnorm(NT * N, 0, sqrt(1 / NT)), N, NT), 1, cumsum))
    
    EX1 <- apply(matrix(stats::rnorm(NT * N, 0, sqrt(1 / NT)), N, NT), 1, cumsum)
    EX2 <- apply(matrix(stats::rnorm(NT * N, 0, sqrt(1 / NT)), N, NT), 1, cumsum)
    
    X1 <- G1 %*% t(MGX1) + UUU1 + t(EX1)
    X2 <- G2 %*% t(MGX2) + UUU2 + t(EX2)
  }
  
  # Generate sparse observations
  Ly_sim1 <- list()
  Lt_sim1 <- list()
  Ly_sim2 <- list()
  Lt_sim2 <- list()
  
  for (i in 1:nrow(X1)) {
    index_sparse <- (1:length(Times))[sort(sample(1:length(Times), nSparse))]
    time_sparse <- Times[index_sparse]
    Ly_sim1[[i]] <- X1[i, index_sparse]
    Lt_sim1[[i]] <- time_sparse
  }
  
  for (i in 1:nrow(X2)) {
    index_sparse <- (1:length(Times))[sort(sample(1:length(Times), nSparse))]
    time_sparse <- Times[index_sparse]
    Ly_sim2[[i]] <- X2[i, index_sparse]
    Lt_sim2[[i]] <- time_sparse
  }
  
  # Create output structure
  if (separate_G) {
    G1_names <- paste0('G1_', 1:ncol(G1))
    G2_names <- paste0('G2_', 1:ncol(G2))
    
    DAT1 <- cbind(G1, X1[, (1:50) * NT / 50])
    DAT1 <- as.data.frame(DAT1)
    names(DAT1) <- c(G1_names, paste0('X1', 1:50))
    
    DAT2 <- cbind(G2, X2[, (1:50) * NT / 50])
    DAT2 <- as.data.frame(DAT2)
    names(DAT2) <- c(G2_names, paste0('X2', 1:50))
    
    DAT <- cbind(G1, G2, X1[, (1:50) * NT / 50], X2[, (1:50) * NT / 50])
    DAT <- as.data.frame(DAT)
    names(DAT) <- c(G1_names, G2_names, paste0('X1', 1:50), paste0('X2', 1:50))
  } else {
    G_names <- paste0('G', 1:ncol(G))
    
    DAT <- cbind(G, X1[, (1:50) * NT / 50], X2[, (1:50) * NT / 50])  
    DAT <- as.data.frame(DAT)
    names(DAT) <- c(G_names, paste0('X1', 1:50), paste0('X2', 1:50))
    
    G1 <- NULL
    G2 <- NULL
  }
  
  RES <- list(
    X1 = list(
      DAT = if(separate_G) DAT1 else DAT[, c(1:J, (J+1):(J+TT))],
      Ly_sim = Ly_sim1,
      Lt_sim = Lt_sim1
    ),
    X2 = list(
      DAT = if(separate_G) DAT2 else DAT[, c(1:J, (J+TT+1):(J+TT+TT))],
      Ly_sim = Ly_sim2,
      Lt_sim = Lt_sim2
    ),
    details = list(
      J = J, 
      TT = TT, 
      NT = NT, 
      G = G,
      G1 = G1,
      G2 = G2,
      times = times, 
      X1 = X1,
      X2 = X2,
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
#' @param mediation_strength Strength of mediation X1 -> X2 (default 0.3)
#' @param separate_G Whether to use separate instruments for each exposure
#' @param shared_G_proportion Proportion of shared instruments (0–1)
#' @param mediation_type Character. Type of mediation effect: "linear" (default), "nonlinear", or "time_varying".
#' 
#' @return List with same structure as getX_multi_exposure()
#' @export
getX_multi_exposure_mediation <- function(N = 10000, 
                                          J = 30,
                                          ZXmodel = 'A',
                                          nSparse = 10,
                                          mediation_strength = 0.3,
                                          separate_G = FALSE,
                                          shared_G_proportion = 0,
                                          mediation_type = "linear") {
  
  NT <- 1000
  TT <- 50 
  
  times <- seq(0, TT, len = (NT + 1))
  Times <- times[-1] 
  times_squ <- times^2 
  
  G <- NULL
  G1 <- NULL
  G2 <- NULL
  
  ## Genetic instruments 
  if (separate_G) {
    J1 <- J
    J2 <- J
    
    if (shared_G_proportion > 0 && shared_G_proportion <= 1) {
      n_shared <- floor(J * shared_G_proportion)
      
      G_shared <- matrix(stats::rbinom(n_shared * N, 2, 0.3), N, n_shared)
      G1_unique <- matrix(stats::rbinom((J1 - n_shared) * N, 2, 0.3), N, J1 - n_shared)
      G2_unique <- matrix(stats::rbinom((J2 - n_shared) * N, 2, 0.3), N, J2 - n_shared)
      
      G1 <- cbind(G_shared, G1_unique)
      G2 <- cbind(G_shared, G2_unique)
      G  <- cbind(G_shared, G1_unique, G2_unique)
    } else {
      G1 <- matrix(stats::rbinom(J1 * N, 2, 0.3), N, J1)
      G2 <- matrix(stats::rbinom(J2 * N, 2, 0.3), N, J2)
      G  <- cbind(G1, G2)
    }
  } else {
    G  <- matrix(stats::rbinom(J * N, 2, 0.3), N, J)
    G1 <- G
    G2 <- G
  }
  
  ## Baseline genetic effects
  noise1 <- stats::rnorm(NT * J, 0, 0.05)
  noise2 <- stats::rnorm(NT * J, 0, 0.05)
  
  a1 <- stats::runif(J, -0.1, 0.1)
  b1 <- stats::runif(J, -0.004, 0.004)
  c1 <- stats::runif(J, 0, 0)
  
  a2 <- stats::runif(J, -0.1, 0.1)
  b2 <- stats::runif(J, -0.004, 0.004)
  c2 <- stats::runif(J, 0, 0)
  
  MGX1 <- t(t(matrix(rep(times_squ[-1], J), NT, J) %*% diag(c1)) +
              t(matrix(rep(times[-1], J), NT, J) %*% diag(b1)) +
              a1) + noise1
  
  MGX2 <- t(t(matrix(rep(times_squ[-1], J), NT, J) %*% diag(c2)) +
              t(matrix(rep(times[-1], J), NT, J) %*% diag(b2)) +
              a2) + noise2
  
  ## Time-varying confounding (shared, as in mediation)
  UUU <- stats::rnorm(N, 0, 1) +
    t(apply(matrix(stats::rnorm(NT * N, 0, sqrt(1 / NT)), N, NT), 1, cumsum))
  
  EX <- apply(matrix(stats::rnorm(NT * N, 0, sqrt(1 / NT)), N, NT), 1, cumsum)
  
  ## X1 generation
  X1 <- G1 %*% t(MGX1) + UUU + t(EX)
  
  ## Mediation effect X1 → X2
  if (mediation_type == "linear") {
    mediation_effect <- mediation_strength * X1
  } else if (mediation_type == "nonlinear") {
    X1s <- scale(X1)
    mediation_effect <- mediation_strength * (X1s + 0.1 * X1s^2)
  } else if (mediation_type == "time_varying") {
    w <- sin(times[-1] * pi / max(times[-1])) + 1
    w <- w / mean(w)
    mediation_effect <- mediation_strength * X1 *
      matrix(rep(w, N), N, NT, byrow = TRUE)
  } else {
    stop("Invalid mediation_type")
  }
  
  ## X2 generation
  X2 <- G2 %*% t(MGX2) + mediation_effect + UUU + t(EX)
  
  ## Sparse observations 
  Ly_sim1 <- Lt_sim1 <- Ly_sim2 <- Lt_sim2 <- list()
  
  for (i in 1:N) {
    idx <- sort(sample(seq_along(Times), nSparse))
    Ly_sim1[[i]] <- X1[i, idx]
    Lt_sim1[[i]] <- Times[idx]
    Ly_sim2[[i]] <- X2[i, idx]
    Lt_sim2[[i]] <- Times[idx]
  }
  
  ## Data frames
  if (separate_G) {
    G1_names <- paste0("G1_", 1:ncol(G1))
    G2_names <- paste0("G2_", 1:ncol(G2))
    
    DAT1 <- as.data.frame(cbind(G1, X1[, (1:50) * NT / 50]))
    names(DAT1) <- c(G1_names, paste0("X1", 1:50))
    
    DAT2 <- as.data.frame(cbind(G2, X2[, (1:50) * NT / 50]))
    names(DAT2) <- c(G2_names, paste0("X2", 1:50))
    
    DAT <- as.data.frame(cbind(G1, G2,
                               X1[, (1:50) * NT / 50],
                               X2[, (1:50) * NT / 50]))
    names(DAT) <- c(G1_names, G2_names,
                    paste0("X1", 1:50),
                    paste0("X2", 1:50))
  } else {
    G_names <- paste0("G", 1:ncol(G))
    DAT <- as.data.frame(cbind(G,
                               X1[, (1:50) * NT / 50],
                               X2[, (1:50) * NT / 50]))
    names(DAT) <- c(G_names, paste0("X1", 1:50), paste0("X2", 1:50))
    G1 <- G2 <- NULL
  }
  
  ## Return object 

  RES <- list(
    X1 = list(
      DAT = if (separate_G) DAT1 else DAT[, c(1:J, (J+1):(J+TT))],
      Ly_sim = Ly_sim1,
      Lt_sim = Lt_sim1
    ),
    X2 = list(
      DAT = if (separate_G) DAT2 else DAT[, c(1:J, (J+TT+1):(J+TT+TT))],
      Ly_sim = Ly_sim2,
      Lt_sim = Lt_sim2
    ),
    details = list(
      J = J,
      TT = TT,
      NT = NT,
      G = G,
      G1 = G1,
      G2 = G2,
      times = times,
      X1 = X1,
      X2 = X2,
      mediation_strength = mediation_strength,
      mediation_type = mediation_type,
      mediation_effect = mediation_effect,
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
#' @param X1Ymodel Effect model for X1 (0-9)
#' @param X2Ymodel Effect model for X2 (0-9)  
#' @param X1_effect Include X1 effect?
#' @param X2_effect Include X2 effect?
#' @param outcome_type "continuous" or "binary"
#' 
#' @return Data frame with outcome Y
#' @export
getY_multi_exposure <- function(RES,
                                X1Ymodel = '1',
                                X2Ymodel = '1',
                                X1_effect = TRUE,
                                X2_effect = TRUE,
                                outcome_type = "continuous") {
  
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
  
  fun1 <- get_fun(as.character(X1Ymodel))
  if(!is.na(X2Ymodel)){
    fun2 <- get_fun(as.character(X2Ymodel))
  }else{
    fun2 <- get_fun('0')
  }
  
  X1 <- RES$details$X1
  X2 <- RES$details$X2
  TT <- RES$details$TT
  NT <- RES$details$NT
  
  times <- seq(0, TT, len = (NT + 1))
  effect_vec1 <- as.numeric(fun1(times[-1]))
  effect_vec2 <- as.numeric(fun2(times[-1]))
  
  if (!X1_effect) effect_vec1 <- rep(0, length(effect_vec1))
  if (!X2_effect) effect_vec2 <- rep(0, length(effect_vec2))
  
  linear_predictor <- (X1 %*% effect_vec1) * TT/NT + (X2 %*% effect_vec2) * TT/NT
  
  if (outcome_type == "continuous") {
    Y <- linear_predictor + stats::rnorm(nrow(X1), 0, 1)
  } else if (outcome_type == "binary") {
    prob_Y <- 1 / (1 + exp(-linear_predictor))
    Y <- stats::rbinom(nrow(X1), 1, prob_Y)
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
      X1[, (1:TT) * NT / TT],
      X2[, (1:TT) * NT / TT]
    )
    names(DAT) <- c(paste0('G', 1:J), paste0('X1', 1:TT), paste0('X2', 1:TT))
  }
  
  DAT$Y <- as.numeric(Y)
  
  attr(DAT, "model_info") <- list(
    X1_model = X1Ymodel,
    X2_model = X2Ymodel,
    X1_effect = X1_effect,
    X2_effect = X2_effect,
    outcome_type = outcome_type,
    effect_vec1 = effect_vec1,
    effect_vec2 = effect_vec2,
    times = times[-1]
  )
  
  return(DAT)
}



