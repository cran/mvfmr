## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## ----install, eval=FALSE------------------------------------------------------
#  # Install from GitHub
#  devtools::install_github("NicoleFontana/mvfmr")

## ----load---------------------------------------------------------------------
library(mvfmr)
library(fdapace)
library(ggplot2)

## ----simulate_data------------------------------------------------------------
set.seed(12345)

# Generate exposure data
sim_data <- getX_multi_exposure(
  N = 300,              # Sample size
  J = 25,               # Number of genetic instruments
  nSparse = 10          # Observations per subject
)

# Check dimensions
cat("Sample size:", nrow(sim_data$details$G), "\n")
cat("Number of instruments:", ncol(sim_data$details$G), "\n")

## ----generate_outcome---------------------------------------------------------
outcome_data <- getY_multi_exposure(
  sim_data,
  X1Ymodel = "2",     # Linear: β₁(t) = 0.02 * t
  X2Ymodel = "5",     # Late-age effect: β₂(t) = 0.1 * (t > 30)
  X1_effect = TRUE,
  X2_effect = TRUE,
  outcome_type = "continuous"
)

cat("Outcome summary:\n")
summary(outcome_data$Y)

## ----fpca---------------------------------------------------------------------
# FPCA for exposure 1
fpca1 <- FPCA(
  sim_data$X1$Ly_sim, 
  sim_data$X1$Lt_sim,
  list(dataType = 'Sparse', error = TRUE, verbose = FALSE)
)

# FPCA for exposure 2
fpca2 <- FPCA(
  sim_data$X2$Ly_sim, 
  sim_data$X2$Lt_sim,
  list(dataType = 'Sparse', error = TRUE, verbose = FALSE)
)

cat("FPCA completed:\n")
cat("  Exposure 1:", fpca1$selectK, "components selected\n")
cat("  Exposure 2:", fpca2$selectK, "components selected\n")

## ----joint_estimation---------------------------------------------------------
result_joint <- mvfmr(
  G = sim_data$details$G,
  fpca_results = list(fpca1, fpca2),
  Y = outcome_data$Y,
  outcome_type = "continuous",
  method = "gmm",
  max_nPC1 = 4,
  max_nPC2 = 4,
  improvement_threshold = 0.001,
  bootstrap = FALSE,
  n_cores = 1,
  true_effects = list(model1 = "2", model2 = "5"),
  X_true = list(X1_true = sim_data$details$X1, X2_true = sim_data$details$X2),
  verbose = FALSE
)

# View results
print(result_joint)

## ----plot_effects, fig.width=10, fig.height=4---------------------------------
# Plot both effects
plot(result_joint)

## ----coefficients-------------------------------------------------------------
# Estimated beta coefficients for basis functions
coef(result_joint)

# Time-varying effects at each time point
head(result_joint$effects$effect1)
head(result_joint$effects$effect2)

## ----performance--------------------------------------------------------------
cat("Performance Metrics:\n")
cat("\nExposure 1 (Linear effect):\n")
cat("  MISE:", round(result_joint$performance$MISE1, 6), "\n")
cat("  Coverage:", round(result_joint$performance$Coverage1, 3), "\n")

cat("\nExposure 2 (Quadratic effect):\n")
cat("  MISE:", round(result_joint$performance$MISE2, 6), "\n")
cat("  Coverage:", round(result_joint$performance$Coverage2, 3), "\n")

## ----separate_estimation------------------------------------------------------
result_separate <- mvfmr_separate(
  G1 = sim_data$details$G,
  G2 = sim_data$details$G,
  fpca_results = list(fpca1, fpca2),
  Y = outcome_data$Y,
  outcome_type = "continuous",
  method = "gmm",
  max_nPC1 = 4,
  max_nPC2 = 4,
  n_cores = 1,
  true_effects = list(model1 = "2", model2 = "5"),
  verbose = FALSE
)

print(result_separate)

## ----comparison_table---------------------------------------------------------
comparison <- data.frame(
  Method = rep(c("Joint (MV-FMR)", "Separate (U-FMR)"), each = 2),
  Exposure = rep(c("X1", "X2"), times = 2),
  MISE = c(
    result_joint$performance$MISE1,
    result_joint$performance$MISE2,
    result_separate$exposure1$performance$MISE,
    result_separate$exposure2$performance$MISE
  ),
  Coverage = c(
    result_joint$performance$Coverage1,
    result_joint$performance$Coverage2,
    result_separate$exposure1$performance$Coverage,
    result_separate$exposure2$performance$Coverage
  )
)

print(comparison)

## ----diagnostics--------------------------------------------------------------
# Calculate F-statistics
K_total <- result_joint$nPC_used$nPC1 + result_joint$nPC_used$nPC2

fstats <- IS(
  J = ncol(sim_data$details$G),
  K = K_total,
  PC = 1:K_total,
  datafull = cbind(
    sim_data$details$G,
    cbind(fpca1$xiEst[, 1:result_joint$nPC_used$nPC1], 
          fpca2$xiEst[, 1:result_joint$nPC_used$nPC2])
  ),
  Y = outcome_data$Y
)

print(fstats)

## ----binary_outcome, eval=FALSE-----------------------------------------------
#  # Generate binary outcome
#  outcome_binary <- getY_multi_exposure(
#    sim_data,
#    X1Ymodel = "2",
#    X2Ymodel = "5",
#    outcome_type = "binary"
#  )
#  
#  # Estimate with control function
#  result_binary <- mvfmr(
#    G = sim_data$details$G,
#    fpca_results = list(fpca1, fpca2),
#    Y = outcome_binary$Y,
#    outcome_type = "binary",
#    method = "cf",  # Control function for binary
#    max_nPC1 = 3,
#    max_nPC2 = 3,
#    n_cores = 1,
#    verbose = FALSE
#  )
#  
#  print(result_binary)

## ----session_info-------------------------------------------------------------
sessionInfo()

