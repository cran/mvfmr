#!/usr/bin/env Rscript
# ============= MULTIVARIABLE FUNCTIONAL MR (MV-FMR) TEST - TWO EXPOSURES # =============
# This script tests the fmvmr package for MULTIVARIABLE analysis
# Demonstrates joint estimation of two correlated time-varying exposures

cat("\n========================================\n")
cat("   TESTING FMVMR PACKAGE - COMPLETE\n")
cat("========================================\n\n")


# Load packages
suppressPackageStartupMessages({
  library(mvfmr)
  library(fdapace)
  library(ggplot2)
  library(gridExtra)
  library(dplyr)
})

cat("Packages loaded\n\n")


# ============= TEST 1: SIMULATE DATA WITH TWO EXPOSURES # =============

cat("TEST 1: Data Simulation (Two Exposures)\n")
cat("----------------------------------------\n")

set.seed(12345)
RES <- getX_multi_exposure(
  N = 500,           # Sample size
  J = 30,            # Number of genetic instruments
  ZXmodel = 'B1',
  nSparse = 10,      # Sparse observations per subject
  n_exposures = 2,   # Number of exposures (m)
  shared_effect = TRUE  # Shared confounding between exposures
)

cat("Exposure data simulated:", nrow(RES$details$G), "subjects\n")
cat("  - Genetic instruments (J):", ncol(RES$details$G), "\n")
cat("  - Observations per subject:", 10, "\n\n")

# ============= TEST 2: FUNCTIONAL PCA # =============

cat("TEST 2: Functional PCA\n")
cat("----------------------\n")

res1 <- FPCA(RES$exposures[[1]]$Ly_sim, RES$exposures[[1]]$Lt_sim, list(dataType = 'Sparse', error = TRUE, verbose = FALSE))
res2 <- FPCA(RES$exposures[[2]]$Ly_sim, RES$exposures[[2]]$Lt_sim, list(dataType = 'Sparse', error = TRUE, verbose = FALSE))

cat("FPCA completed\n")
cat("  - Exposure 1:", res1$selectK, "components\n")
cat("  - Exposure 2:", res2$selectK, "components\n\n")

# ============= TEST 3: OUTCOME SIMULATION # =============

cat("TEST 3: Outcome Simulation\n")
cat("--------------------------\n")

DAT <- getY_multi_exposure(
  RES,
  XYmodels = c("2", "8"),  # Linear effect for exposure 1, quadratic for exposure 2
  X_effects = c(TRUE, TRUE),
  outcome_type = "continuous"
)

cat("Outcome generated (continuous)\n")
cat("  - Effect model X1: Linear (model '2')\n")
cat("  - Effect model X2: Quadratic (model '8')\n")
cat("  - Outcome mean:", round(mean(DAT$Y), 3), "\n")
cat("  - Outcome SD:", round(sd(DAT$Y), 3), "\n\n")

# ============= TEST 4: JOINT MULTIVARIABLE ESTIMATION (MV-FMR) # =============

cat("TEST 4: Joint Multivariable Estimation (mvfmr)\n")
cat("-----------------------------------------------\n")

result_joint <- mvfmr(
  G = RES$details$G,
  fpca_results = list(res1, res2),
  Y = DAT$Y,
  outcome_type = "continuous",
  method = "gmm",
  max_nPC = c(5, 5),
  improvement_threshold = 0.001,
  bootstrap = FALSE,
  n_cores = 2,
  true_effects = c("2", "8"),
  X_true = RES$details$X_list,
  verbose = FALSE
)

cat("Joint estimation completed!\n")
cat("  - Components selected: nPC1 =", result_joint$nPC_used[1], ", nPC2 =", result_joint$nPC_used[2], "\n\n")

# ============= TEST 5: SEPARATE UNIVARIABLE ESTIMATION (U-FMR) # =============

cat("TEST 5: Separate Univariable Estimation (mvfmr_separate)\n")
cat("---------------------------------------------------------\n")

result_separate <- mvfmr_separate(
  G_list = list(RES$details$G, RES$details$G),
  fpca_results = list(res1, res2),
  Y = DAT$Y,
  outcome_type = "continuous",
  method = "gmm",
  max_nPC = c(5, 5),
  improvement_threshold = 0.001,
  bootstrap = FALSE,
  n_cores = 2,
  true_effects = c("2", "8"),
  verbose = FALSE
)

cat("Separate estimation completed!\n")
cat("  - Components X1:", result_separate$exposures[[1]]$nPC_used, "\n")
cat("  - Components X2:", result_separate$exposures[[2]]$nPC_used, "\n\n")

# ============= TEST 6: INSTRUMENT STRENGTH DIAGNOSTICS # =============

cat("TEST 6: Conditional F-statistics\n")
cat("---------------------------------\n")

# Calculate F-statistics for joint estimation
fstats <- IS(
  J = 30,
  K = (res1$selectK + res2$selectK),
  PC = 1:(res1$selectK + res2$selectK),
  datafull = cbind(
    RES$details$G,
    cbind(res1$xiEst[, 1:res1$selectK], res2$xiEst[, 1:res2$selectK])),
  Y = DAT$Y)

fstats_df <- cbind(
  "Exposure" = c(rep("X1", res1$selectK), rep("X2", res2$selectK)),
  as.data.frame(fstats)) %>% dplyr::select(Exposure, PC, cFF)

cat("\nConditional F-statistics:\n")
print(fstats_df)

# ============= TEST 7: PERFORMANCE COMPARISON # =============
cat("TEST 7: Performance Comparison\n")
cat("------------------------------\n")

cat("\nJoint Estimation (MV-FMR):\n")
cat("  Exposure 1 - MISE:", round(result_joint$performance$MISE[[1]], 6), "\n")
cat("  Exposure 1 - Coverage:", round(result_joint$performance$Coverage[[1]], 3), "\n")
cat("  Exposure 2 - MISE:", round(result_joint$performance$MISE[[2]], 6), "\n")
cat("  Exposure 2 - Coverage:", round(result_joint$performance$Coverage[[2]], 3), "\n")

cat("\nSeparate Estimation (U-FMR):\n")
cat("  Exposure 1 - MISE:", round(result_separate$exposures[[1]]$performance$MISE, 6), "\n")
cat("  Exposure 1 - Coverage:", round(result_separate$exposures[[1]]$performance$Coverage, 3), "\n")
cat("  Exposure 2 - MISE:", round(result_separate$exposures[[2]]$performance$MISE, 6), "\n")
cat("  Exposure 2 - Coverage:", round(result_separate$exposures[[2]]$performance$Coverage, 3), "\n")

# ============= TEST 8: VISUALIZATION - BUILT-IN PLOTS # =============

cat("TEST 8: Visualization - Built-in Plots\n")
cat("---------------------------------------\n")

# Use the built-in plot method for joint estimation
cat("\nDisplaying built-in plots for joint estimation:\n")
plot(result_joint)

cat("Displaying built-in plots for separate estimation:\n")
plot(result_separate)

cat("\n")

# ============= TEST 9: CUSTOM VISUALIZATION - JOINT ESTIMATION # =============

cat("TEST 9: Custom Visualization - Joint Estimation\n")
cat("------------------------------------------------\n")

# Extract data for plotting (one data frame per exposure)
ggdata_joint_1 <- result_joint$raw_result$ggdata[[1]]
ggdata_joint_2 <- result_joint$raw_result$ggdata[[2]]

# Create custom plots for both exposures
p_joint_1 <- ggplot(ggdata_joint_1, aes(x = time)) +
  geom_line(aes(y = effect), linewidth = 1, color = "darkblue") +
  geom_ribbon(aes(ymin = effect_low, ymax = effect_up),
              alpha = 0.2, fill = "darkblue") +
  geom_line(aes(y = true_shape), linewidth = 1,
            color = "red", linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
  labs(
    title = "Exposure 1: Joint Estimation (MV-FMR)",
    subtitle = paste0("Linear effect - ", result_joint$nPC_used[1], " components"),
    x = "Time / Age",
    y = "Causal Effect beta1(t)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10),
    panel.grid.minor = element_blank()
  )

p_joint_2 <- ggplot(ggdata_joint_2, aes(x = time)) +
  geom_line(aes(y = effect), linewidth = 1, color = "darkgreen") +
  geom_ribbon(aes(ymin = effect_low, ymax = effect_up),
              alpha = 0.2, fill = "darkgreen") +
  geom_line(aes(y = true_shape), linewidth = 1,
            color = "red", linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
  labs(
    title = "Exposure 2: Joint Estimation (MV-FMR)",
    subtitle = paste0("Quadratic effect - ", result_joint$nPC_used[2], " components"),
    x = "Time / Age",
    y = "Causal Effect beta2(t)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10),
    panel.grid.minor = element_blank()
  )

# Display joint plots
grid.arrange(p_joint_1, p_joint_2, ncol = 2, top = "Joint Multivariable Estimation (MV-FMR)")

cat("Custom joint estimation plots displayed\n\n")

# ============= TEST 10: CUSTOM VISUALIZATION - SEPARATE ESTIMATION # =============

cat("TEST 10: Custom Visualization - Separate Estimation\n")
cat("----------------------------------------------------\n")

# Extract data for separate estimation
ggdata_sep1 <- result_separate$raw_result$ggdata[[1]]
ggdata_sep2 <- result_separate$raw_result$ggdata[[2]]

# Create plots for separate estimation
p_sep_1 <- ggplot(ggdata_sep1, aes(x = time)) +
  geom_line(aes(y = effect), linewidth = 1, color = "purple") +
  geom_ribbon(aes(ymin = effect_low, ymax = effect_up),
              alpha = 0.2, fill = "purple") +
  geom_line(aes(y = true_shape), linewidth = 1,
            color = "red", linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
  labs(
    title = "Exposure 1: Separate Estimation (U-FMR)",
    subtitle = paste0("Linear effect - ", result_separate$exposures[[1]]$nPC_used, " components"),
    x = "Time / Age",
    y = "Causal Effect β1(t)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10),
    panel.grid.minor = element_blank()
  )

p_sep_2 <- ggplot(ggdata_sep2, aes(x = time)) +
  geom_line(aes(y = effect), linewidth = 1, color = "orange") +
  geom_ribbon(aes(ymin = effect_low, ymax = effect_up),
              alpha = 0.2, fill = "orange") +
  geom_line(aes(y = true_shape), linewidth = 1,
            color = "red", linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
  labs(
    title = "Exposure 2: Separate Estimation (U-FMR)",
    subtitle = paste0("Quadratic effect - ", result_separate$exposures[[2]]$nPC_used, " components"),
    x = "Time / Age",
    y = "Causal Effect β2(t)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10),
    panel.grid.minor = element_blank()
  )

# Display separate plots
grid.arrange(p_sep_1, p_sep_2, ncol = 2, top = "Separate Univariable Estimation (U-FMR)")

cat("Custom separate estimation plots displayed\n\n")

# ============= TEST 11: COMPARISON PLOT - JOINT VS SEPARATE # =============

cat("TEST 11: Comparison Plot - Joint vs Separate\n")
cat("---------------------------------------------\n")

# Prepare comparison data for Exposure 1
comp_data_1 <- data.frame(
  time = ggdata_joint_1$time,
  joint = ggdata_joint_1$effect,
  separate = ggdata_sep1$effect,
  true = ggdata_joint_1$true_shape
)

comp_data_1_long <- tidyr::pivot_longer(
  comp_data_1,
  cols = c(joint, separate, true),
  names_to = "Method",
  values_to = "Effect"
)

# Comparison plot for Exposure 1
p_comp_1 <- ggplot(comp_data_1_long, aes(x = time, y = Effect, color = Method)) +
  geom_line(linewidth = 1) +
  scale_color_manual(
    values = c("joint" = "darkgreen", "separate" = "darkblue", "true" = "red"),
    labels = c("joint" = "Joint (MV-FMR)", "separate" = "Separate (U-FMR)",
               "true" = "True Effect")
  ) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
  labs(
    title = "Exposure 1: Comparison of Joint vs Separate Estimation",
    x = "Time / Age",
    y = "Causal Effect beta1(t)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

# Prepare comparison data for Exposure 2
comp_data_2 <- data.frame(
  time = ggdata_joint_2$time,
  joint = ggdata_joint_2$effect,
  separate = ggdata_sep2$effect,
  true = ggdata_joint_2$true_shape
)

comp_data_2_long <- tidyr::pivot_longer(
  comp_data_2,
  cols = c(joint, separate, true),
  names_to = "Method",
  values_to = "Effect"
)

# Comparison plot for Exposure 2
p_comp_2 <- ggplot(comp_data_2_long, aes(x = time, y = Effect, color = Method)) +
  geom_line(linewidth = 1) +
  scale_color_manual(
    values = c("joint" = "darkgreen", "separate" = "darkblue", "true" = "red"),
    labels = c("joint" = "Joint (MV-FMR)", "separate" = "Separate (U-FMR)",
               "true" = "True Effect")
  ) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
  labs(
    title = "Exposure 2: Comparison of Joint vs Separate Estimation",
    x = "Time / Age",
    y = "Causal Effect beta2(t)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

# Display comparison plots
grid.arrange(p_comp_1, p_comp_2, ncol = 2,
             top = "Joint vs Separate Estimation Comparison")

cat("Comparison plots displayed\n\n")

# ============= TEST 12: BINARY OUTCOME # =============

cat("TEST 12: Binary Outcome\n")
cat("-----------------------\n")

DAT_binary <- getY_multi_exposure(
  RES,
  XYmodels = c("2", "8"),
  X_effects = c(TRUE, TRUE),
  outcome_type = "binary"
)

cat("Binary outcome generated\n")
cat("  - Cases (Y=1):", sum(DAT_binary$Y == 1), "\n")
cat("  - Controls (Y=0):", sum(DAT_binary$Y == 0), "\n\n")

result_binary <- mvfmr(
  G = RES$details$G,
  fpca_results = list(res1, res2),
  Y = DAT_binary$Y,
  outcome_type = "binary",
  method = "cf",
  max_nPC = c(3, 3),
  n_cores = 2,
  verbose = FALSE
)

cat("Binary estimation completed!\n")
cat("  - Components: nPC1 =", result_binary$nPC_used[1], ", nPC2 =", result_binary$nPC_used[2], "\n\n")

# ============= TEST 13: EXPORT RESULTS # =============

cat("TEST 13: Export Results\n")
cat("-----------------------\n")

output_dir <- tempdir()

# Export joint estimation effects
joint_effects <- data.frame(
  time = ggdata_joint_1$time,
  effect1 = ggdata_joint_1$effect,
  effect1_lower = ggdata_joint_1$effect_low,
  effect1_upper = ggdata_joint_1$effect_up,
  effect2 = ggdata_joint_2$effect,
  effect2_lower = ggdata_joint_2$effect_low,
  effect2_upper = ggdata_joint_2$effect_up
)

write.csv(joint_effects, file.path(output_dir, "mvfmr_joint_effects.csv"), row.names = FALSE)
cat("Joint effects saved to:", file.path(output_dir, "mvfmr_joint_effects.csv"), "\n")

# Export coefficients
coef_data_joint <- data.frame(coefficient = coef(result_joint))
write.csv(coef_data_joint, file.path(output_dir, "mvfmr_coefficients.csv"), row.names = FALSE)
cat("Coefficients saved to:", file.path(output_dir, "mvfmr_coefficients.csv"), "\n")

# Export performance comparison
performance_comparison <- data.frame(
  Method = rep(c("Joint (MV-FMR)", "Separate (U-FMR)"), each = 2),
  Exposure = rep(c("Exposure 1", "Exposure 2"), times = 2),
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

write.csv(performance_comparison, file.path(output_dir, "mvfmr_performance.csv"), row.names = FALSE)
cat("Performance comparison saved to:", file.path(output_dir, "mvfmr_performance.csv"), "\n\n")
