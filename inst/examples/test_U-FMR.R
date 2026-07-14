#!/usr/bin/env Rscript
# =============================================================================
# UNIVARIABLE FUNCTIONAL MR (U-FMR) TEST - SINGLE EXPOSURE
# =============================================================================
# This script tests the fmvmr package for SINGLE EXPOSURE analysis
# Demonstrates how to use mvfmr_separate() for univariable functional MR

cat("\n========================================\n")
cat("   TESTING U-FMR - SINGLE EXPOSURE\n")
cat("========================================\n\n")

# Load packages
suppressPackageStartupMessages({
  library(mvfmr)
  library(fdapace)
  library(ggplot2)
  library(dplyr)
})

cat("Packages loaded\n\n")

# =============================================================================
# TEST 1: SIMULATE DATA WITH SINGLE EXPOSURE
# =============================================================================

cat("TEST 1: Data Simulation (Single Exposure)\n")
cat("------------------------------------------\n")

set.seed(12345)

# Simulate a single exposure (n_exposures = 1)
RES <- getX_multi_exposure(
  N = 500,           # Sample size
  J = 30,            # Number of genetic instruments
  nSparse = 10,      # Sparse observations per subject
  n_exposures = 1,   # Single exposure
  shared_effect = FALSE  # No shared confounding (single exposure analysis)
)

cat("Exposure data simulated:", nrow(RES$details$G), "subjects\n")
cat("  - Genetic instruments (J):", ncol(RES$details$G), "\n")
cat("  - Observations per subject:", 10, "\n\n")

# =============================================================================
# TEST 2: FUNCTIONAL PCA FOR EXPOSURE 1
# =============================================================================

cat("TEST 2: Functional PCA for Exposure 1\n")
cat("--------------------------------------\n")

# Perform FPCA on the (only) exposure
res1 <- FPCA(
  RES$exposures[[1]]$Ly_sim,
  RES$exposures[[1]]$Lt_sim,
  list(dataType = 'Sparse', error = TRUE, verbose = FALSE)
)

cat("FPCA completed for Exposure 1\n")
cat("  - Automatic components selected:", res1$selectK, "\n")
cat("  - Variance explained:", round(sum(res1$lambda[1:res1$selectK]) / sum(res1$lambda) * 100, 2), "%\n\n")

# =============================================================================
# TEST 3: GENERATE OUTCOME (ONLY EXPOSURE 1 HAS EFFECT)
# =============================================================================

cat("TEST 3: Outcome Simulation\n")
cat("--------------------------\n")

# Generate outcome where the exposure has a causal effect
DAT <- getY_multi_exposure(
  RES,
  XYmodels = "2",     # Linear time-varying effect
  X_effects = TRUE,   # Include the effect
  outcome_type = "continuous"
)

cat("Outcome generated (continuous)\n")
cat("  - Effect model for X1: Linear (model '2')\n")
cat("  - Outcome mean:", round(mean(DAT$Y), 3), "\n")
cat("  - Outcome SD:", round(sd(DAT$Y), 3), "\n\n")

# =============================================================================
# TEST 4: UNIVARIABLE FUNCTIONAL MR ESTIMATION
# =============================================================================

cat("TEST 4: Univariable Functional MR Estimation\n")
cat("---------------------------------------------\n")

# Estimate causal effect of EXPOSURE 1 ONLY using mvfmr_separate()
# with a G_list of length 1
result <- mvfmr_separate(
  G_list = list(RES$details$G),  # Genetic instruments for the exposure
  fpca_results = list(res1),
  Y = DAT$Y,
  outcome_type = "continuous",
  method = "gmm",               # Generalized Method of Moments
  max_nPC = 5,                  # Maximum components to consider
  improvement_threshold = 0.001,
  bootstrap = FALSE,            # Set TRUE for bootstrap inference
  n_cores = 2,
  true_effects = "2",           # True effect for validation
  verbose = FALSE
)

cat("Estimation completed!\n")
cat("  - Components selected for X1:", result$exposures[[1]]$nPC_used, "\n")
cat("  - Method: GMM (Generalized Method of Moments)\n\n")

# =============================================================================
# TEST 5: DISPLAY RESULTS
# =============================================================================

cat("TEST 5: Results Summary\n")
cat("-----------------------\n\n")

# Print result object
print(result)

cat("\n")

# Display coefficients
cat("Estimated Coefficients (Basis Functions):\n")
print(round(coef(result, exposure = 1), 4))

cat("\n")

# =============================================================================
# TEST 6: PERFORMANCE METRICS
# =============================================================================

cat("TEST 6: Performance Metrics\n")
cat("---------------------------\n")

cat("\nExposure 1 Performance:\n")
cat("  - MISE (Mean Integrated Squared Error):", round(result$exposures[[1]]$performance$MISE, 6), "\n")
cat("  - Coverage Rate:", round(result$exposures[[1]]$performance$Coverage, 3), "\n")


# =============================================================================
# TEST 7: INSTRUMENT STRENGTH DIAGNOSTICS
# =============================================================================

cat("TEST 7: Instrument Strength (F-statistics)\n")
cat("-------------------------------------------\n")

# Calculate F-statistics for exposure 1
fstats <- IS(
  J = ncol(RES$details$G),
  K = res1$selectK,
  PC = 1:res1$selectK,
  datafull = cbind(RES$details$G, res1$xiEst[, 1:res1$selectK]),
  Y = DAT$Y
)

fstats_df <- as.data.frame(fstats) %>%
  mutate(Component = paste0("PC", PC)) %>%
  select(Component, RR, FF, cFF)

cat("\nF-statistics for Exposure 1:\n")
print(fstats_df)

# =============================================================================
# TEST 8: VISUALIZE TIME-VARYING EFFECT
# =============================================================================

cat("TEST 8: Visualization\n")
cat("---------------------\n")

# Plot the estimated time-varying effect
plot(result)

cat("Effect plot displayed\n")
cat("  - Black line: Estimated effect\n")
cat("  - Dashed lines: 95% confidence interval\n")
cat("  - Blue line: True effect (if provided)\n\n")

# Create custom plot with more details
ggdata <- result$raw_result$ggdata[[1]]

custom_plot <- ggplot(ggdata, aes(x = time)) +
  geom_line(aes(y = effect), linewidth = 1, color = "black") +
  geom_ribbon(aes(ymin = effect_low, ymax = effect_up),
              alpha = 0.2, fill = "blue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  labs(
    title = "Time-Varying Causal Effect of Exposure 1",
    subtitle = paste0("U-FMR with ", result$exposures[[1]]$nPC_used, " components"),
    x = "Time / Age",
    y = "Causal Effect β(t)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11),
    panel.grid.minor = element_blank()
  )

if (!is.na(ggdata$true_shape[1])) {
  custom_plot <- custom_plot +
    geom_line(aes(y = true_shape), color = "red", linewidth = 1, linetype = "dashed")
}

print(custom_plot)

cat("\nCustom plot created\n\n")

# =============================================================================
# TEST 9: QUADRATIC EFFECT MODEL
# =============================================================================

cat("TEST 9: Quadratic Effect Model\n")
cat("-------------------------------\n")

# Generate new outcome with quadratic effect
DAT_quad <- getY_multi_exposure(
  RES,
  XYmodels = "8",     # Quadratic effect (inverted U-shape)
  X_effects = TRUE,
  outcome_type = "continuous"
)

# Estimate
result_quad <- mvfmr_separate(
  G_list = list(RES$details$G),
  fpca_results = list(res1),
  Y = DAT_quad$Y,
  outcome_type = "continuous",
  method = "gmm",
  max_nPC = 5,
  true_effects = "8",
  verbose = FALSE
)

cat("Quadratic model estimated\n")
cat("  - Components selected:", result_quad$exposures[[1]]$nPC_used, "\n")
cat("  - MISE:", round(result_quad$exposures[[1]]$performance$MISE, 6), "\n")
cat("  - Coverage:", round(result_quad$exposures[[1]]$performance$Coverage, 3), "\n\n")

# =============================================================================
# TEST 10: BINARY OUTCOME
# =============================================================================

cat("TEST 10: Binary Outcome\n")
cat("-----------------------\n")

# Generate binary outcome
DAT_binary <- getY_multi_exposure(
  RES,
  XYmodels = "2",     # Linear effect
  X_effects = TRUE,
  outcome_type = "binary"
)

cat("Binary outcome generated\n")
cat("  - Cases (Y=1):", sum(DAT_binary$Y == 1), "\n")
cat("  - Controls (Y=0):", sum(DAT_binary$Y == 0), "\n")
cat("  - Prevalence:", round(mean(DAT_binary$Y), 3), "\n\n")

# Estimate with control function method
result_binary <- mvfmr_separate(
  G_list = list(RES$details$G),
  fpca_results = list(res1),
  Y = DAT_binary$Y,
  outcome_type = "binary",
  method = "cf",      # Control function for binary outcomes
  max_nPC = 3,
  verbose = FALSE
)

cat("Binary outcome estimation completed\n")
cat("  - Method: Control Function (2SRI)\n")
cat("  - Components selected:", result_binary$exposures[[1]]$nPC_used, "\n\n")

# =============================================================================
# TEST 11: BOOTSTRAP INFERENCE
# =============================================================================

cat("TEST 11: Bootstrap Inference\n")
cat("----------------------------\n")
cat("Note: Bootstrap is computationally intensive\n")
cat("Running with small number of replicates (10) for demonstration\n\n")

result_bootstrap <- mvfmr_separate(
  G_list = list(RES$details$G),
  fpca_results = list(res1),
  Y = DAT$Y,
  outcome_type = "continuous",
  method = "gmm",
  max_nPC = 5,
  bootstrap = TRUE,
  n_bootstrap = 10,  # Use 100-200 for real analysis
  n_cores = 2,
  true_effects = "2",
  verbose = FALSE
)

cat("Bootstrap inference completed\n")
cat("  - Bootstrap replicates: 10\n")
cat("  - Bootstrap confidence intervals computed\n\n")

# =============================================================================
# TEST 12: DIFFERENT EFFECT MODELS COMPARISON
# =============================================================================

cat("TEST 12: Comparing Different Effect Models\n")
cat("-------------------------------------------\n")

effect_models <- c("1", "2", "8")  # Constant, Linear up, Quadratic
model_names <- c("Constant", "Linear Increasing", "Linear Decreasing", "Quadratic")

comparison_results <- data.frame(
  Model = character(),
  Name = character(),
  nPC = integer(),
  MISE = numeric(),
  Coverage = numeric(),
  stringsAsFactors = FALSE
)

for (i in seq_along(effect_models)) {
  cat("  Testing model", effect_models[i], ":", model_names[i], "...\n")

  # Generate outcome
  outcome_temp <- getY_multi_exposure(
    RES,
    XYmodels = effect_models[i],
    X_effects = TRUE,
    outcome_type = "continuous"
  )

  # Estimate
  result_temp <- mvfmr_separate(
    G_list = list(RES$details$G),
    fpca_results = list(res1),
    Y = outcome_temp$Y,
    outcome_type = "continuous",
    method = "gmm",
    max_nPC = 5,
    true_effects = effect_models[i],
    verbose = FALSE
  )

  # Store results
  comparison_results <- rbind(comparison_results, data.frame(
    Model = effect_models[i],
    Name = model_names[i],
    nPC = result_temp$exposures[[1]]$nPC_used,
    MISE = result_temp$exposures[[1]]$performance$MISE,
    Coverage = result_temp$exposures[[1]]$performance$Coverage
  ))
}

cat("\n")
cat("Comparison of Effect Models:\n")
print(comparison_results)

# =============================================================================
# TEST 13: EXTRACT AND EXPORT RESULTS
# =============================================================================

cat("TEST 13: Extract and Export Results\n")
cat("------------------------------------\n")

output_dir <- tempdir()

# Extract time-varying effect curve
effect_curve <- data.frame(
  time = result$raw_result$ggdata[[1]]$time,
  effect = result$raw_result$ggdata[[1]]$effect,
  lower_ci = result$raw_result$ggdata[[1]]$effect_low,
  upper_ci = result$raw_result$ggdata[[1]]$effect_up
)

# Save results to CSV
write.csv(effect_curve, file.path(output_dir, "ufmr_effect_curve.csv"), row.names = FALSE)
cat("Results saved to:", file.path(output_dir, "ufmr_effect_curve.csv"), "\n")

# Save coefficients
coef_df <- data.frame(
  Component = paste0("Beta_", 1:length(result$exposures[[1]]$coefficients)),
  Coefficient = result$exposures[[1]]$coefficients,
  SE = sqrt(diag(result$exposures[[1]]$vcov))
)

write.csv(coef_df, file.path(output_dir, "ufmr_coefficients.csv"), row.names = FALSE)
cat("Coefficients saved to:", file.path(output_dir, "ufmr_coefficients.csv"), "\n\n")
