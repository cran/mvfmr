# ============= TEST FOR THE 3 PAPER SCENARIOS # =============
# Tests the 3 main scenarios presented in the paper (parameters simplified for computing them faster):
# - Scenario 1: Shared confounding with G1, G2, G12 (15% overlap)
# - Scenario 2: No X2 effect with shared confounding
# - Scenario 3: Mediation scenario
#
# Compares Multi vs Separate approaches using MISE and Coverage

library(mvfmr)
library(fdapace)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)

# ============= SETUP # =============

set.seed(12345)

# Simulation parameters
N <- 500  # Smaller N for testing (paper uses 5000)
J <- 30   # Smaller J for testing (paper uses 100)
nSparse <- 10
n_sim <- 1  # Smaller n_sim for testing (paper uses 500)

# Results storage
results_all <- data.frame() # For performance
isres_all <- data.frame() # For cFF/FF

# ============= SCENARIO 1: PLEIOTROPY MODEL # =============
# Paper Section 4a
# Genetic variants G12 influence both the exposures X1 and X2
# shared_G_proportion = 0.15 (15% of instruments are shared)
# Both X1 and X2 have causal effects on Y

cat("\n========================================\n")
cat("SCENARIO 1: Shared Confounding (G1 & G2 & G12)\n")
cat("========================================\n\n")

X1Ymodel_vec <- c("2", "8")  # Linear and Quadratic
X2Ymodel_vec <- c("2", "8")

XY_grid <- expand.grid(X1 = X1Ymodel_vec, X2 = X2Ymodel_vec)

for (i in 1:nrow(XY_grid)) {
  for (sim in 1:n_sim) {
    cat(paste0("Scenario 1 - Grid ", i, "/", nrow(XY_grid), " - Sim ", sim, "/", n_sim, "\n"))

    set.seed(sim * 12345)

    # Generate data with shared confounding
    RES <- getX_multi_exposure(
      N = N,
      J = J,
      nSparse = nSparse,
      n_exposures = 2,
      shared_effect = TRUE,
      separate_G = TRUE,
      shared_G_proportion = 0.15
    )

    # FPCA
    res1 <- FPCA(RES$exposures[[1]]$Ly_sim, RES$exposures[[1]]$Lt_sim, list(dataType = 'Sparse', error = TRUE, verbose = FALSE))
    res2 <- FPCA(RES$exposures[[2]]$Ly_sim, RES$exposures[[2]]$Lt_sim, list(dataType = 'Sparse', error = TRUE, verbose = FALSE))

    # Generate outcome
    DAT <- getY_multi_exposure(
      RES,
      XYmodels = c(XY_grid[i, "X1"], XY_grid[i, "X2"]),
      X_effects = c(TRUE, TRUE),
      outcome_type = "continuous"
    )

    # Calculate Instrument Strength (IS)

    # For MULTIVARIABLE approach: cFF (conditional F-statistic)
    K_total <- res1$selectK + res2$selectK
    IS_multi <- IS(
      J = ncol(RES$details$G),
      K = K_total,
      PC = 1:K_total,
      datafull = cbind(
        RES$details$G,
        cbind(res1$xiEst[, 1:res1$selectK], res2$xiEst[, 1:res2$selectK])
      )
    )

    IS_multi_df <- as.data.frame(IS_multi) %>%
      mutate(
        Exposure = c(rep("X1", res1$selectK), rep("X2", res2$selectK)),
        nPC = c(1:res1$selectK, 1:res2$selectK),
        Scenario = "1_Shared_G1_G2_G12",
        Simulation = sim,
        X1Ymodel = XY_grid[i, "X1"],
        X2Ymodel = XY_grid[i, "X2"],
        Approach = "Multivariable"
      )

    # For SEPARATE approach: FF (standard F-statistic)
    # Uses separate G matrices for each exposure
    IS_sep1 <- IS(
      J = ncol(RES$details$G_list[[1]]),
      K = res1$selectK,
      PC = 1:res1$selectK,
      datafull = cbind(RES$details$G_list[[1]], res1$xiEst[, 1:res1$selectK])
    )

    IS_sep2 <- IS(
      J = ncol(RES$details$G_list[[2]]),
      K = res2$selectK,
      PC = 1:res2$selectK,
      datafull = cbind(RES$details$G_list[[2]], res2$xiEst[, 1:res2$selectK])
    )

    IS_sep_df <- rbind(
      as.data.frame(IS_sep1) %>%
        mutate(
          Exposure = "X1",
          nPC = 1:res1$selectK,
          Scenario = "1_Shared_G1_G2_G12",
          Simulation = sim,
          X1Ymodel = XY_grid[i, "X1"],
          X2Ymodel = XY_grid[i, "X2"],
          Approach = "Separate"
        ),
      as.data.frame(IS_sep2) %>%
        mutate(
          Exposure = "X2",
          nPC = 1:res2$selectK,
          Scenario = "1_Shared_G1_G2_G12",
          Simulation = sim,
          X1Ymodel = XY_grid[i, "X1"],
          X2Ymodel = XY_grid[i, "X2"],
          Approach = "Separate"
        )
    )

    # Combine IS results
    isres_all <- rbind(isres_all, IS_multi_df, IS_sep_df)

    # Multi approach (joint estimation)
    result_multi <- mvfmr(
      G = RES$details$G,
      fpca_results = list(res1, res2),
      Y = DAT$Y,
      outcome_type = "continuous",
      method = "gmm",
      true_effects = c(XY_grid[i, "X1"], XY_grid[i, "X2"]),
      verbose = FALSE
    )

    # Separate approach (univariable estimation)
    result_separate <- mvfmr_separate(
      G_list = list(RES$details$G_list[[1]], RES$details$G_list[[2]]),
      fpca_results = list(res1, res2),
      Y = DAT$Y,
      outcome_type = "continuous",
      method = "gmm",
      true_effects = c(XY_grid[i, "X1"], XY_grid[i, "X2"]),
      verbose = FALSE
    )

    # Store results
    temp_res <- data.frame(
      Scenario = "Scenario 1: Pleiotropy model",
      Simulation = sim,
      X1Ymodel = XY_grid[i, "X1"],
      X2Ymodel = XY_grid[i, "X2"],
      Multi_nPC1 = result_multi$nPC_used[1],
      Multi_nPC2 = result_multi$nPC_used[2],
      Multi_MISE1 = result_multi$performance$MISE[[1]],
      Multi_Coverage1 = result_multi$performance$Coverage[[1]],
      Multi_MISE2 = result_multi$performance$MISE[[2]],
      Multi_Coverage2 = result_multi$performance$Coverage[[2]],
      Separate_nPC1 = result_separate$exposures[[1]]$nPC_used,
      Separate_nPC2 = result_separate$exposures[[2]]$nPC_used,
      Separate_MISE1 = result_separate$exposures[[1]]$performance$MISE,
      Separate_Coverage1 = result_separate$exposures[[1]]$performance$Coverage,
      Separate_MISE2 = result_separate$exposures[[2]]$performance$MISE,
      Separate_Coverage2 = result_separate$exposures[[2]]$performance$Coverage
    )

    results_all <- rbind(results_all, temp_res)
  }
}

# ============= SCENARIO 2: NULL EFFECT CONTROL (beta_2 = 0) # =============
# Paper Section 5
# X1 has effect, X2 has NO effect on Y
# Tests robustness when one exposure is null

cat("\n========================================\n")
cat("SCENARIO 2: No X2 Effect (beta_2 = 0)\n")
cat("========================================\n\n")

X1Ymodel_vec <- c("2", "8")
X2Ymodel_vec <- c("0")  # Null effect
XY_grid <- expand.grid(X1 = X1Ymodel_vec, X2 = X2Ymodel_vec)

for (i in 1:nrow(XY_grid)) {
  for (sim in 1:n_sim) {
    cat(paste0("Scenario 2 - Grid ", i, "/", nrow(XY_grid), " - Sim ", sim, "/", n_sim, "\n"))

    set.seed(sim * 23456)

    # Generate data
    RES <- getX_multi_exposure(
      N = N,
      J = J,
      nSparse = nSparse,
      n_exposures = 2,
      shared_effect = TRUE,
      separate_G = TRUE,
      shared_G_proportion = 0.15
    )

    # FPCA
    res1 <- FPCA(RES$exposures[[1]]$Ly_sim, RES$exposures[[1]]$Lt_sim, list(dataType = 'Sparse', error = TRUE, verbose = FALSE))
    res2 <- FPCA(RES$exposures[[2]]$Ly_sim, RES$exposures[[2]]$Lt_sim, list(dataType = 'Sparse', error = TRUE, verbose = FALSE))

    # Generate outcome - X2 has NO effect
    DAT <- getY_multi_exposure(
      RES,
      XYmodels = c(XY_grid[i, "X1"], XY_grid[i, "X2"]),  # X2 model "0" = null
      X_effects = c(TRUE, FALSE),  # X2 has no effect
      outcome_type = "continuous"
    )

    # Calculate IS statistics
    K_total <- res1$selectK + res2$selectK
    IS_multi <- IS(
      J = ncol(RES$details$G),
      K = K_total,
      PC = 1:K_total,
      datafull = cbind(
        RES$details$G,
        cbind(res1$xiEst[, 1:res1$selectK], res2$xiEst[, 1:res2$selectK])
      )
    )

    IS_multi_df <- as.data.frame(IS_multi) %>%
      mutate(
        Exposure = c(rep("X1", res1$selectK), rep("X2", res2$selectK)),
        nPC = c(1:res1$selectK, 1:res2$selectK),
        Scenario = "2_No_X2_Effect",
        Simulation = sim,
        X1Ymodel = XY_grid[i, "X1"],
        X2Ymodel = XY_grid[i, "X2"],
        Approach = "Multivariable"
      )

    IS_sep1 <- IS(
      J = ncol(RES$details$G_list[[1]]),
      K = res1$selectK,
      PC = 1:res1$selectK,
      datafull = cbind(RES$details$G_list[[1]], res1$xiEst[, 1:res1$selectK])
    )

    IS_sep2 <- IS(
      J = ncol(RES$details$G_list[[2]]),
      K = res2$selectK,
      PC = 1:res2$selectK,
      datafull = cbind(RES$details$G_list[[2]], res2$xiEst[, 1:res2$selectK])
    )

    IS_sep_df <- rbind(
      as.data.frame(IS_sep1) %>%
        mutate(
          Exposure = "X1",
          nPC = 1:res1$selectK,
          Scenario = "2_No_X2_Effect",
          Simulation = sim,
          X1Ymodel = XY_grid[i, "X1"],
          X2Ymodel = XY_grid[i, "X2"],
          Approach = "Separate"
        ),
      as.data.frame(IS_sep2) %>%
        mutate(
          Exposure = "X2",
          nPC = 1:res2$selectK,
          Scenario = "2_No_X2_Effect",
          Simulation = sim,
          X1Ymodel = XY_grid[i, "X1"],
          X2Ymodel = XY_grid[i, "X2"],
          Approach = "Separate"
        )
    )

    isres_all <- rbind(isres_all, IS_multi_df, IS_sep_df)

    # Multi approach
    result_multi <- mvfmr(
      G = RES$details$G,
      fpca_results = list(res1, res2),
      Y = DAT$Y,
      outcome_type = "continuous",
      method = "gmm",
      true_effects = c(XY_grid[i, "X1"], XY_grid[i, "X2"])
    )

    # Separate approach
    result_separate <- mvfmr_separate(
      G_list = list(RES$details$G_list[[1]], RES$details$G_list[[2]]),
      fpca_results = list(res1, res2),
      Y = DAT$Y,
      outcome_type = "continuous",
      method = "gmm",
      true_effects = c(XY_grid[i, "X1"], XY_grid[i, "X2"])
    )

    # Store results

    temp_res <- data.frame(
      Scenario = "Scenario 2: Null effect control",
      Simulation = sim,
      X1Ymodel = XY_grid[i, "X1"],
      X2Ymodel = XY_grid[i, "X2"],
      Multi_nPC1 = result_multi$nPC_used[1],
      Multi_nPC2 = result_multi$nPC_used[2],
      Multi_MISE1 = result_multi$performance$MISE[[1]],
      Multi_Coverage1 = result_multi$performance$Coverage[[1]],
      Multi_MISE2 = result_multi$performance$MISE[[2]],
      Multi_Coverage2 = result_multi$performance$Coverage[[2]],
      Separate_nPC1 = result_separate$exposures[[1]]$nPC_used,
      Separate_nPC2 = result_separate$exposures[[2]]$nPC_used,
      Separate_MISE1 = result_separate$exposures[[1]]$performance$MISE,
      Separate_Coverage1 = result_separate$exposures[[1]]$performance$Coverage,
      Separate_MISE2 = result_separate$exposures[[2]]$performance$MISE,
      Separate_Coverage2 = result_separate$exposures[[2]]$performance$Coverage
    )

    results_all <- rbind(results_all, temp_res)
  }
}

# ============= SCENARIO 3: MEDIATION MODEL # =============
# Paper Section 6
# X1 affects X2 (mediation), both affect Y
# Tests ability to handle correlated exposures with causal pathway

cat("\n========================================\n")
cat("SCENARIO 3: Mediation (X1 -> X2 -> Y)\n")
cat("========================================\n\n")


X1Ymodel_vec <- c("2", "8")
X2Ymodel_vec <- c("2", "8")

XY_grid <- expand.grid(X1 = X1Ymodel_vec, X2 = X2Ymodel_vec)

# Mediation: exposure 1 mediates onto exposure 2 with strength 0.3.
# mediation_strength[j, k] is the strength with which exposure j mediates
# its effect onto exposure k (only entries with j < k may be nonzero).
mediation_strength <- matrix(0, 2, 2)
mediation_strength[1, 2] <- 0.3

for (i in 1:nrow(XY_grid)) {
  for (sim in 1:n_sim) {
    cat(paste0("Scenario 3 - Grid ", i, "/", nrow(XY_grid),
               " - Sim ", sim, "/", n_sim, "\n"))

    set.seed(sim * 34567)

    # Generate data with mediation
    RES <- getX_multi_exposure_mediation(
      N = N,
      J = J,
      nSparse = nSparse,
      n_exposures = 2,
      mediation_strength = mediation_strength,
      separate_G = TRUE,
      shared_G_proportion = 0.15
    )

    # FPCA
    res1 <- FPCA(RES$exposures[[1]]$Ly_sim, RES$exposures[[1]]$Lt_sim, list(dataType = 'Sparse', error = TRUE, verbose = FALSE))
    res2 <- FPCA(RES$exposures[[2]]$Ly_sim, RES$exposures[[2]]$Lt_sim, list(dataType = 'Sparse', error = TRUE, verbose = FALSE))

    # Generate outcome
    DAT <- getY_multi_exposure(
      RES,
      XYmodels = c(XY_grid[i, "X1"], XY_grid[i, "X2"]),
      X_effects = c(TRUE, TRUE),
      outcome_type = "continuous"
    )

    # Calculate IS statistics
    K_total <- res1$selectK + res2$selectK
    IS_multi <- IS(
      J = ncol(RES$details$G),
      K = K_total,
      PC = 1:K_total,
      datafull = cbind(
        RES$details$G,
        cbind(res1$xiEst[, 1:res1$selectK], res2$xiEst[, 1:res2$selectK])
      )
    )

    IS_multi_df <- as.data.frame(IS_multi) %>%
      mutate(
        Exposure = c(rep("X1", res1$selectK), rep("X2", res2$selectK)),
        nPC = c(1:res1$selectK, 1:res2$selectK),
        Scenario = "3_Mediation",
        Simulation = sim,
        X1Ymodel = XY_grid[i, "X1"],
        X2Ymodel = XY_grid[i, "X2"],
        Approach = "Multivariable"
      )

    IS_sep1 <- IS(
      J = ncol(RES$details$G_list[[1]]),
      K = res1$selectK,
      PC = 1:res1$selectK,
      datafull = cbind(RES$details$G_list[[1]], res1$xiEst[, 1:res1$selectK])
    )

    IS_sep2 <- IS(
      J = ncol(RES$details$G_list[[2]]),
      K = res2$selectK,
      PC = 1:res2$selectK,
      datafull = cbind(RES$details$G_list[[2]], res2$xiEst[, 1:res2$selectK])
    )

    IS_sep_df <- rbind(
      as.data.frame(IS_sep1) %>%
        mutate(
          Exposure = "X1",
          nPC = 1:res1$selectK,
          Scenario = "3_Mediation",
          Simulation = sim,
          X1Ymodel = XY_grid[i, "X1"],
          X2Ymodel = XY_grid[i, "X2"],
          Approach = "Separate"
        ),
      as.data.frame(IS_sep2) %>%
        mutate(
          Exposure = "X2",
          nPC = 1:res2$selectK,
          Scenario = "3_Mediation",
          Simulation = sim,
          X1Ymodel = XY_grid[i, "X1"],
          X2Ymodel = XY_grid[i, "X2"],
          Approach = "Separate"
        )
    )

    isres_all <- rbind(isres_all, IS_multi_df, IS_sep_df)

    # Multi approach
    result_multi <- mvfmr(
      G = RES$details$G,
      fpca_results = list(res1, res2),
      Y = DAT$Y,
      outcome_type = "continuous",
      method = "gmm",
      true_effects = c(XY_grid[i, "X1"], XY_grid[i, "X2"])
    )

    # Separate approach
    result_separate <- mvfmr_separate(
      G_list = list(RES$details$G_list[[1]], RES$details$G_list[[2]]),
      fpca_results = list(res1, res2),
      Y = DAT$Y,
      outcome_type = "continuous",
      method = "gmm",
      true_effects = c(XY_grid[i, "X1"], XY_grid[i, "X2"])
    )

    # Store results
    temp_res <- data.frame(
      Scenario = "Scenario 3: Mediation model",
      Simulation = sim,
      X1Ymodel = XY_grid[i, "X1"],
      X2Ymodel = XY_grid[i, "X2"],
      Multi_nPC1 = result_multi$nPC_used[1],
      Multi_nPC2 = result_multi$nPC_used[2],
      Multi_MISE1 = result_multi$performance$MISE[[1]],
      Multi_Coverage1 = result_multi$performance$Coverage[[1]],
      Multi_MISE2 = result_multi$performance$MISE[[2]],
      Multi_Coverage2 = result_multi$performance$Coverage[[2]],
      Separate_nPC1 = result_separate$exposures[[1]]$nPC_used,
      Separate_nPC2 = result_separate$exposures[[2]]$nPC_used,
      Separate_MISE1 = result_separate$exposures[[1]]$performance$MISE,
      Separate_Coverage1 = result_separate$exposures[[1]]$performance$Coverage,
      Separate_MISE2 = result_separate$exposures[[2]]$performance$MISE,
      Separate_Coverage2 = result_separate$exposures[[2]]$performance$Coverage
    )

    results_all <- rbind(results_all, temp_res)
  }
}

# ============= SAVE RESULTS # =============

# Use tempdir()
output_dir <- tempdir()
write.csv(results_all, file.path(output_dir, "paper_scenarios_results.csv"), row.names = FALSE)
write.csv(isres_all, file.path(output_dir, "paper_scenarios_IS_statistics.csv"), row.names = FALSE)

cat("\nResults saved to:\n")
cat("  -", file.path(output_dir, "paper_scenarios_results.csv"), "\n")
cat("  -", file.path(output_dir, "paper_scenarios_IS_statistics.csv"), "\n")

# ============= SUMMARY STATISTICS # =============

cat("\n========================================\n")
cat("SUMMARY STATISTICS\n")
cat("========================================\n\n")

summary_stats <- results_all %>%
  group_by(Scenario, X1Ymodel, X2Ymodel) %>%
  summarize(
    # Multi (MV-FMR)
    Multi_MISE1_mean = mean(Multi_MISE1, na.rm = TRUE),
    Multi_MISE1_sd = sd(Multi_MISE1, na.rm = TRUE),
    Multi_Coverage1_mean = mean(Multi_Coverage1, na.rm = TRUE),
    Multi_MISE2_mean = mean(Multi_MISE2, na.rm = TRUE),
    Multi_MISE2_sd = sd(Multi_MISE2, na.rm = TRUE),
    Multi_Coverage2_mean = mean(Multi_Coverage2, na.rm = TRUE),

    # Separate (U-FMR)
    Separate_MISE1_mean = mean(Separate_MISE1, na.rm = TRUE),
    Separate_MISE1_sd = sd(Separate_MISE1, na.rm = TRUE),
    Separate_Coverage1_mean = mean(Separate_Coverage1, na.rm = TRUE),
    Separate_MISE2_mean = mean(Separate_MISE2, na.rm = TRUE),
    Separate_MISE2_sd = sd(Separate_MISE2, na.rm = TRUE),
    Separate_Coverage2_mean = mean(Separate_Coverage2, na.rm = TRUE),

    .groups = 'drop'
  )

print(summary_stats)

# IS Statistics summary
cat("\n========================================\n")
cat("INSTRUMENT STRENGTH (IS) SUMMARY\n")
cat("========================================\n\n")

is_summary <- isres_all %>%
  group_by(Scenario, Approach, Exposure, nPC) %>%
  summarize(
    FF_mean = mean(FF, na.rm = TRUE),
    FF_sd = sd(FF, na.rm = TRUE),
    cFF_mean = mean(cFF, na.rm = TRUE),
    cFF_sd = sd(cFF, na.rm = TRUE),
    .groups = 'drop'
  )

print(is_summary)

# ============= VISUALIZATION - MISE COMPARISON # =============

cat("\n========================================\n")
cat("CREATING PLOTS\n")
cat("========================================\n\n")

# Helper function to label models
model_to_label <- function(model_num) {
  case_when(
    model_num == "2" ~ "Linear",
    model_num == "8" ~ "Quadratic",
    model_num == "0" ~ "Null",
    TRUE ~ model_num
  )
}

results_all$X1_label <- model_to_label(results_all$X1Ymodel)
results_all$X2_label <- model_to_label(results_all$X2Ymodel)
results_all$ModelLabel <- paste0(results_all$X1_label, " vs ", results_all$X2_label)

# Reshape for plotting - MISE
mise_long <- results_all %>%
  select(Scenario, Simulation, ModelLabel,
         Multi_MISE1, Multi_MISE2, Separate_MISE1, Separate_MISE2) %>%
  pivot_longer(
    cols = c(Multi_MISE1, Multi_MISE2, Separate_MISE1, Separate_MISE2),
    names_to = "Variable",
    values_to = "MISE"
  ) %>%
  mutate(
    Exposure = ifelse(grepl("MISE1", Variable), "X1", "X2"),
    Method = ifelse(grepl("Multi", Variable), "MV-FMR", "U-FMR")
  )

# MISE plot
p_mise <- ggplot(mise_long, aes(x = ModelLabel, y = MISE, fill = Method)) +
  geom_boxplot() +
  facet_grid(Exposure ~ Scenario, scales = "free_y") +
  theme_bw() +
  labs(
    title = "Integrated Squared Error: MV-FMR vs U-FMR",
    subtitle = paste0("N=", N, ", J=", J, ", n_sim=", n_sim),
    y = "MISE",
    x = "Effect Models"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    strip.text = element_text(size = 9),
    legend.position = "bottom"
  ) +
  scale_fill_manual(values = c("MV-FMR" = "#E41A1C", "U-FMR" = "#377EB8"))

# Reshape for plotting - Coverage
coverage_long <- results_all %>%
  select(Scenario, Simulation, ModelLabel,
         Multi_Coverage1, Multi_Coverage2,
         Separate_Coverage1, Separate_Coverage2) %>%
  pivot_longer(
    cols = c(Multi_Coverage1, Multi_Coverage2,
             Separate_Coverage1, Separate_Coverage2),
    names_to = "Variable",
    values_to = "Coverage"
  ) %>%
  mutate(
    Exposure = ifelse(grepl("Coverage1", Variable), "X1", "X2"),
    Method = ifelse(grepl("Multi", Variable), "MV-FMR", "U-FMR")
  )

# Coverage plot
p_coverage <- ggplot(coverage_long, aes(x = ModelLabel, y = Coverage, fill = Method)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  facet_grid(Exposure ~ Scenario, scales = "free_y") +
  theme_bw() +
  labs(
    title = "Coverage Rate: MV-FMR vs U-FMR",
    subtitle = paste0("N=", N, ", J=", J, ", n_sim=", n_sim),
    y = "Coverage Rate",
    x = "Effect Models"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    strip.text = element_text(size = 9),
    legend.position = "bottom"
  ) +
  scale_fill_manual(values = c("MV-FMR" = "#E41A1C", "U-FMR" = "#377EB8"))

# Save plots
ggsave(file.path(output_dir, "paper_scenarios_MISE.png"), p_mise, width = 14, height = 8)
ggsave(file.path(output_dir, "paper_scenarios_Coverage.png"), p_coverage, width = 14, height = 8)


cat("\nPlots saved:\n")
cat("  -", file.path(output_dir, "paper_scenarios_MISE.png"), "\n")
cat("  -", file.path(output_dir, "paper_scenarios_Coverage.png"), "\n")

# ============= COMPARISON TABLE (LIKE PAPER TABLE 1) # =============

# Format function
format_mean_sd <- function(mean_val, sd_val, digits = 3) {
  paste0(format(round(mean_val, digits), nsmall = digits),
         " (", format(round(sd_val, digits), nsmall = digits), ")")
}

# Create comparison table
comparison_table <- summary_stats %>%
  mutate(
    Multi_MISE1 = format_mean_sd(Multi_MISE1_mean, Multi_MISE1_sd),
    Multi_Coverage1 = format(round(Multi_Coverage1_mean, 3), nsmall = 3),
    Multi_MISE2 = format_mean_sd(Multi_MISE2_mean, Multi_MISE2_sd),
    Multi_Coverage2 = format(round(Multi_Coverage2_mean, 3), nsmall = 3),
    Separate_MISE1 = format_mean_sd(Separate_MISE1_mean, Separate_MISE1_sd),
    Separate_Coverage1 = format(round(Separate_Coverage1_mean, 3), nsmall = 3),
    Separate_MISE2 = format_mean_sd(Separate_MISE2_mean, Separate_MISE2_sd),
    Separate_Coverage2 = format(round(Separate_Coverage2_mean, 3), nsmall = 3)
  ) %>%
  select(Scenario, X1Ymodel, X2Ymodel,
         Multi_MISE1, Multi_Coverage1, Multi_MISE2, Multi_Coverage2,
         Separate_MISE1, Separate_Coverage1, Separate_MISE2, Separate_Coverage2)

cat("\n========================================\n")
cat("COMPARISON TABLE (Paper Style)\n")
cat("========================================\n\n")
print(comparison_table)

# Save table
write.csv(comparison_table, file.path(output_dir, "paper_scenarios_table.csv"), row.names = FALSE)
cat("\nTable saved to:", file.path(output_dir, "paper_scenarios_table.csv"), "\n")
