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
#  
#  # Or, for the development version:
#  # devtools::install_github("NicoleFontana/mvfmr")

## ----load---------------------------------------------------------------------
library(mvfmr)
library(fdapace)
library(ggplot2)

## ----simulate-----------------------------------------------------------------
set.seed(12345)

# Generate exposure data with a single exposure (m = 1)
sim_data <- getX_multi_exposure(
  N = 300,
  J = 25,
  nSparse = 10,
  n_exposures = 1
)

cat("Data simulated:\n")
cat("  Sample size:", nrow(sim_data$details$G), "\n")
cat("  Instruments:", ncol(sim_data$details$G), "\n")

## ----outcome------------------------------------------------------------------
outcome_data <- getY_multi_exposure(
  sim_data,
  XYmodels = "2",     # Linear effect: beta(t) = 0.02*t
  X_effects = TRUE,
  outcome_type = "continuous"
)

cat("Outcome summary:\n")
summary(outcome_data$Y)

## ----fpca---------------------------------------------------------------------
fpca1 <- FPCA(
  sim_data$exposures[[1]]$Ly_sim,
  sim_data$exposures[[1]]$Lt_sim,
  list(dataType = 'Sparse', error = TRUE, verbose = FALSE)
)

cat("FPCA completed:\n")
cat("  Components selected:", fpca1$selectK, "\n")
cat("  Variance explained:",
    round(sum(fpca1$lambda[1:fpca1$selectK]) / sum(fpca1$lambda) * 100, 1), "%\n")

## ----estimation---------------------------------------------------------------
result <- mvfmr_separate(
  G_list = list(sim_data$details$G),
  fpca_results = list(fpca1),
  Y = outcome_data$Y,
  outcome_type = "continuous",
  method = "gmm",
  max_nPC = 4,
  n_cores = 1,
  true_effects = "2",
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
cat("\nFirst 10 time points of beta(t):\n")
head(result$exposures[[1]]$effect, 10)

## ----performance--------------------------------------------------------------
cat("Performance:\n")
cat("  MISE:", round(result$exposures[[1]]$performance$MISE, 6), "\n")
cat("  Coverage:", round(result$exposures[[1]]$performance$Coverage, 3), "\n")
cat("  Components used:", result$exposures[[1]]$nPC_used, "\n")

## ----binary, eval=FALSE-------------------------------------------------------
#  # Generate binary outcome
#  outcome_binary <- getY_multi_exposure(
#    sim_data,
#    XYmodels = "2",
#    X_effects = TRUE,
#    outcome_type = "binary"
#  )
#  
#  # Estimate with control function
#  result_binary <- mvfmr_separate(
#    G_list = list(sim_data$details$G),
#    fpca_results = list(fpca1),
#    Y = outcome_binary$Y,
#    outcome_type = "binary",
#    method = "cf",      # Control function for binary
#    max_nPC = 3,
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
#    G_list = list(sim_data$details$G),
#    fpca_results = list(fpca1),
#    Y = outcome_data$Y,
#    outcome_type = "continuous",
#    bootstrap = TRUE,
#    n_bootstrap = 100,
#    max_nPC = 4,
#    verbose = FALSE
#  )
#  
#  # Bootstrap CIs are stored in result_boot$exposures[[1]]$...

## ----twosample, eval=FALSE----------------------------------------------------
#  # Simulate GWAS summary statistics
#  by_outcome <- rnorm(25, 0.02, 0.01)
#  sy_outcome <- runif(25, 0.005, 0.015)
#  
#  result_2sample <- fmvmr_separate_twosample(
#    G_list = list(sim_data$details$G),
#    fpca_results = list(fpca1),
#    by_outcome_list = list(by_outcome),
#    sy_outcome_list = list(sy_outcome),
#    ny_outcome = 50000,
#    max_nPC = 3,
#    verbose = FALSE
#  )
#  
#  print(result_2sample)

## ----session------------------------------------------------------------------
sessionInfo()

