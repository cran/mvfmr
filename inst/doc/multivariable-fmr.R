## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## ----install, eval=FALSE------------------------------------------------------
#  # Install from CRAN
#  install.packages("mvfmr")

## ----load---------------------------------------------------------------------
library(mvfmr)
library(fdapace)
library(ggplot2)

## ----simulate_data------------------------------------------------------------
set.seed(473920)

# Generate exposure data
sim_data <- getX_multi_exposure(
  N = 300,              # Sample size
  J = 25,               # Number of genetic instruments
  nSparse = 10,          # Observations per subject
  n_exposures = 2        # Number of exposures (m)
)

# Check dimensions
cat("Sample size:", nrow(sim_data$details$G), "\n")
cat("Number of instruments:", ncol(sim_data$details$G), "\n")

## ----generate_outcome---------------------------------------------------------
outcome_data <- getY_multi_exposure(
  sim_data,
  XYmodels = c("2", "2"),     # Exposure 1/2: linear beta(t) = 0.02*t
  X_effects = c(TRUE, TRUE),
  outcome_type = "continuous"
)

cat("Outcome summary:\n")
summary(outcome_data$Y)

## ----fpca---------------------------------------------------------------------
fpca_results <- lapply(sim_data$exposures, function(exp_k) {
  FPCA(
    exp_k$Ly_sim,
    exp_k$Lt_sim,
    list(dataType = 'Sparse', error = TRUE, verbose = FALSE)
  )
})

cat("FPCA completed:\n")
for (k in seq_along(fpca_results)) {
  cat("  Exposure", k, ":", fpca_results[[k]]$selectK, "components selected\n")
}

## ----joint_estimation---------------------------------------------------------
result_joint <- mvfmr(
  G = sim_data$details$G,
  fpca_results = fpca_results,
  Y = outcome_data$Y,
  outcome_type = "continuous",
  method = "gmm",
  max_nPC = c(4, 4),
  improvement_threshold = 0.001,
  bootstrap = FALSE,
  n_cores = 1,
  true_effects = c("2", "2"),
  X_true = sim_data$details$X_list,
  verbose = FALSE
)

# View results
print(result_joint)

## ----plot_effects, fig.width=10, fig.height=4---------------------------------
# Plot every exposure's effect
plot(result_joint)

## ----coefficients-------------------------------------------------------------
# Estimated beta coefficients for basis functions
coef(result_joint)

# Time-varying effects at each time point (one entry per exposure)
head(result_joint$effects[[1]])
head(result_joint$effects[[2]])

## ----performance--------------------------------------------------------------
cat("Performance Metrics:\n")
for (k in seq_along(result_joint$effects)) {
  cat("\nExposure", k, ":\n")
  cat("  MISE:", round(result_joint$performance$MISE[[k]], 6), "\n")
  cat("  Coverage:", round(result_joint$performance$Coverage[[k]], 3), "\n")
}

## ----separate_estimation------------------------------------------------------
result_separate <- mvfmr_separate(
  G_list = list(sim_data$details$G, sim_data$details$G),
  fpca_results = fpca_results,
  Y = outcome_data$Y,
  outcome_type = "continuous",
  method = "gmm",
  max_nPC = c(4, 4),
  n_cores = 1,
  true_effects = c("2", "2"),
  verbose = FALSE
)

print(result_separate)

## ----comparison_table---------------------------------------------------------
comparison <- data.frame(
  Method = rep(c("Joint (MV-FMR)", "Separate (U-FMR)"), each = 2),
  Exposure = rep(c("X1", "X2"), times = 2),
  MISE = c(
    result_joint$performance$MISE[[1]],
    result_joint$performance$MISE[[2]],
    result_separate$exposures[[1]]$performance$MISE,
    result_separate$exposures[[2]]$performance$MISE
  ),
  Coverage = c(
    result_joint$performance$Coverage[[1]],
    result_joint$performance$Coverage[[2]],
    result_separate$exposures[[1]]$performance$Coverage,
    result_separate$exposures[[2]]$performance$Coverage
  )
)

print(comparison)

## ----diagnostics--------------------------------------------------------------
# Calculate F-statistics
K_total <- sum(result_joint$nPC_used)

PC_stacked <- do.call(cbind, lapply(seq_along(fpca_results), function(k) {
  fpca_results[[k]]$xiEst[, 1:result_joint$nPC_used[k]]
}))

fstats <- IS(
  J = ncol(sim_data$details$G),
  K = K_total,
  PC = 1:K_total,
  datafull = cbind(sim_data$details$G, PC_stacked),
  Y = outcome_data$Y
)

fstats_df <- cbind(
  "Exposure" = unlist(lapply(seq_along(result_joint$nPC_used), function(k) {
    rep(paste0("X", k), result_joint$nPC_used[k])
  })),
  as.data.frame(fstats)
)

print(fstats_df[, c("Exposure", "PC", "cFF")])

## ----binary_outcome, eval=FALSE-----------------------------------------------
#  # Generate binary outcome
#  outcome_binary <- getY_multi_exposure(
#    sim_data,
#    XYmodels = c("2", "2"),
#    outcome_type = "binary"
#  )
#  
#  # Estimate with control function
#  result_binary <- mvfmr(
#    G = sim_data$details$G,
#    fpca_results = fpca_results,
#    Y = outcome_binary$Y,
#    outcome_type = "binary",
#    method = "cf",  # Control function for binary
#    max_nPC = c(3, 3),
#    n_cores = 1,
#    verbose = FALSE
#  )
#  
#  print(result_binary)

## ----m3_example---------------------------------------------------------------
set.seed(163918)#281046
sim_data3 <- getX_multi_exposure(N = 500, J = 50, nSparse = 10, n_exposures = 3)

outcome_data3 <- getY_multi_exposure(
  sim_data3,
  XYmodels = c("2", "2", "2"),
  outcome_type = "continuous"
)

fpca_results3 <- lapply(sim_data3$exposures, function(exp_k) {
  FPCA(exp_k$Ly_sim, exp_k$Lt_sim, list(dataType = 'Sparse', error = TRUE, verbose = FALSE))
})

result_joint3 <- mvfmr(
  G = sim_data3$details$G,
  fpca_results = fpca_results3,
  Y = outcome_data3$Y,
  outcome_type = "continuous",
  method = "gmm",
  max_nPC = c(4, 4, 4),
  n_cores = 1,
  true_effects = c("2", "2", "2"),
  X_true = sim_data3$details$X_list,
  verbose = FALSE
)

print(result_joint3)

## ----m3_plot, fig.width=12, fig.height=4--------------------------------------
plot(result_joint3)

## ----session_info-------------------------------------------------------------
sessionInfo()

