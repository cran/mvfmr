## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## ----install, eval=FALSE------------------------------------------------------
#  devtools::install_github("NicoleFontana/mvfmr")

## ----load---------------------------------------------------------------------
library(mvfmr)
library(fdapace)
library(ggplot2)

## ----simulate-----------------------------------------------------------------
set.seed(12345)

# Generate exposure data (we'll only use X1)
sim_data <- getX_multi_exposure(
  N = 300,
  J = 25,
  nSparse = 10
)

cat("Data simulated:\n")
cat("  Sample size:", nrow(sim_data$details$G), "\n")
cat("  Instruments:", ncol(sim_data$details$G), "\n")

## ----outcome------------------------------------------------------------------
outcome_data <- getY_multi_exposure(
  sim_data,
  X1Ymodel = "2",     # Linear effect: β(t) = 0.02×t
  X2Ymodel = "0",     # No effect from X2
  X1_effect = TRUE,
  X2_effect = FALSE,  # X2 does not affect Y
  outcome_type = "continuous"
)

cat("Outcome summary:\n")
summary(outcome_data$Y)

## ----fpca---------------------------------------------------------------------
# FPCA for exposure 1 only
fpca1 <- FPCA(
  sim_data$X1$Ly_sim, 
  sim_data$X1$Lt_sim,
  list(dataType = 'Sparse', error = TRUE, verbose = FALSE)
)

cat("FPCA completed:\n")
cat("  Components selected:", fpca1$selectK, "\n")
cat("  Variance explained:", 
    round(sum(fpca1$lambda[1:fpca1$selectK]) / sum(fpca1$lambda) * 100, 1), "%\n")

## ----estimation---------------------------------------------------------------
result <- mvfmr_separate(
  G1 = sim_data$details$G,   # Instruments for X1
  G2 = NULL,                  # No second exposure
  fpca_results = list(fpca1),
  Y = outcome_data$Y,
  outcome_type = "continuous",
  method = "gmm",
  max_nPC1 = 4,
  max_nPC2 = 4,  # Not used when G2 = NULL
  n_cores = 1,
  true_effects = list(model1 = "2", model2 = "0"),
  verbose = FALSE
)

print(result)

## ----plot_effect, fig.width=7, fig.height=5-----------------------------------
plot(result)

## ----extract------------------------------------------------------------------
# Coefficients for basis functions
cat("Estimated coefficients:\n")
print(round(coef(result, exposure = 1), 4))

# Time-varying effect curve
cat("\nFirst 10 time points of β(t):\n")
head(result$exposure1$effect, 10)

## ----performance--------------------------------------------------------------
cat("Performance:\n")
cat("  MISE:", round(result$exposure1$performance$MISE, 6), "\n")
cat("  Coverage:", round(result$exposure1$performance$Coverage, 3), "\n")
cat("  Components used:", result$exposure1$nPC_used, "\n")

## ----binary, eval=FALSE-------------------------------------------------------
#  # Generate binary outcome
#  outcome_binary <- getY_multi_exposure(
#    sim_data,
#    X1Ymodel = "2",
#    X2Ymodel = "0",
#    X1_effect = TRUE,
#    X2_effect = FALSE,
#    outcome_type = "binary"
#  )
#  
#  # Estimate with control function
#  result_binary <- mvfmr_separate(
#    G1 = sim_data$details$G,
#    G2 = NULL,
#    fpca_results = list(fpca1),
#    Y = outcome_binary$Y,
#    outcome_type = "binary",
#    method = "cf",      # Control function for binary
#    max_nPC1 = 3,
#    n_cores = 1,
#    verbose = FALSE
#  )
#  
#  print(result_binary)
#  cat("Cases:", sum(outcome_binary$Y == 1), "\n")
#  cat("Controls:", sum(outcome_binary$Y == 0), "\n")

## ----bootstrap, eval=FALSE----------------------------------------------------
#  # Get robust confidence intervals via bootstrap
#  result_boot <- mvfmr_separate(
#    G1 = sim_data$details$G,
#    G2 = NULL,
#    fpca_results = list(fpca1),
#    Y = outcome_data$Y,
#    outcome_type = "continuous",
#    bootstrap = TRUE,
#    n_bootstrap = 100,
#    max_nPC1 = 4,
#    verbose = FALSE
#  )
#  
#  # Bootstrap CIs are stored in result_boot$exposure1$...

## ----twosample, eval=FALSE----------------------------------------------------
#  # Simulate GWAS summary statistics
#  by_outcome <- rnorm(25, 0.02, 0.01)
#  sy_outcome <- runif(25, 0.005, 0.015)
#  
#  result_2sample <- fmvmr_separate_twosample(
#    G1_exposure = sim_data$details$G,
#    G2_exposure = NULL,
#    fpca_results = list(fpca1),
#    by_outcome1 = by_outcome,
#    by_outcome2 = NULL,
#    sy_outcome1 = sy_outcome,
#    sy_outcome2 = NULL,
#    ny_outcome = 50000,
#    max_nPC1 = 3,
#    verbose = FALSE
#  )
#  
#  print(result_2sample)

## ----session------------------------------------------------------------------
sessionInfo()

