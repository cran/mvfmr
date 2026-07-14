# mvfmr: Multivariable Functional Mendelian Randomization

[![arXiv](https://img.shields.io/badge/arXiv-2512.19064-b31b1b.svg)](https://arxiv.org/abs/2512.19064)

## Citation

If you use this package, please cite:

```
Fontana, N., Ieva, F., Zuccolo, L., Di Angelantonio, E., & Secchi, P. (2025). Unraveling time-varying causal effects of multiple exposures: Integrating functional data analysis with multivariable Mendelian randomization. arXiv. https://arxiv.org/abs/2512.19064
```

```bibtex
@misc{fontana2025unravelingtimevaryingcausaleffects,
      title={Unraveling time-varying causal effects of multiple exposures: integrating Functional Data Analysis with Multivariable Mendelian Randomization}, 
      author={Nicole Fontana and Francesca Ieva and Luisa Zuccolo and Emanuele Di Angelantonio and Piercesare Secchi},
      year={2025},
      eprint={2512.19064},
      archivePrefix={arXiv},
      primaryClass={stat.AP},
      url={https://arxiv.org/abs/2512.19064}, 
}
```

## Overview

The `mvfmr` package implements Multivariable Functional Mendelian Randomization methods to estimate time-varying causal effects of longitudinal exposures on health outcomes. The package supports:

- **Multivariable Functional MR (MV-FMR)**: Joint estimation of an arbitrary number of correlated time-varying exposures
- **Univariable Functional MR (U-FMR)**: Separate estimation of each exposure independently (including single-exposure analysis)
- **Both continuous and binary outcomes**
- **One-sample and two-sample MR designs** (the latter using outcome GWAS summary statistics)
- **Automatic component selection** via cross-validation
- **Mediation pathways** between exposures
- **Bootstrap inference** for confidence intervals
- **Instrument strength diagnostics**

## Installation

Install the released version from CRAN:

```r
install.packages("mvfmr")
```

Or install the development version from GitHub:

```r
install.packages("devtools") # Install devtools if not already installed
devtools::install_github("NicoleFontana/mvfmr")
```

## Test Scripts and Simulations

The package ships with test scripts demonstrating different use cases. Once installed, they live inside the package itself, so you can run or open them with `demo()` / `system.file()` — no download needed:

### 1. **Manuscript Simulations** (`demo/tests_manuscript.R`)
Reproduces the main simulation scenarios from the manuscript:
- **Scenario 1-3**: pleiotropy, null effects, mediation
- **Exposure effects**: linear and quadratic
- **Performance comparison**: MV-FMR vs U-FMR across scenarios
- **Evaluation**: MISE, coverage rates

```r
# Run the manuscript simulations
demo("tests_manuscript", package = "mvfmr")
```

### 2. **Multivariable FMR Tutorial** (`inst/examples/test_MV-FMR.R`)
Complete tutorial for using **joint multivariable** estimation:
- Data simulation with two exposures
- FPCA and component selection
- Joint estimation with `mvfmr()`
- Instrument diagnostics
- Performance metrics and visualization
- Binary outcome analysis
- Comparison with univariable estimation

```r
# Run the MV-FMR tutorial
source(system.file("examples", "test_MV-FMR.R", package = "mvfmr"))

# Or just locate it to open/copy/adapt it:
system.file("examples", "test_MV-FMR.R", package = "mvfmr")
```

### 3. **Univariable FMR Tutorial** (`inst/examples/test_U-FMR.R`)
Complete tutorial for **single exposure** analysis:
- Single exposure simulation
- Univariable estimation with `mvfmr_separate()`
- Instrument diagnostics
- Performance metrics and visualization
- Comparison of different exposure effects
- Binary outcome analysis

```r
# Run the U-FMR tutorial
source(system.file("examples", "test_U-FMR.R", package = "mvfmr"))
```

**Note**: These scripts serve as templates for your own analyses. Open them (via `system.file()`), modify the parameters, and adapt to your specific research questions.

## Quick Start

### Example 1: Joint Multivariable FMR (two exposures)

Exposure-related arguments (`fpca_results`, `max_nPC`, `true_effects`, `X_true`, ...) are lists or vectors of length `m`, one entry per exposure.

```r
library(mvfmr)
library(fdapace)

# Step 1: Simulate exposure data (m = 2 exposures)
set.seed(12345)
sim_data <- getX_multi_exposure(
  N = 1000,           # Sample size
  J = 50,             # Number of genetic instruments
  nSparse = 10,        # Sparse observations per subject
  n_exposures = 2      # Number of exposures (m)
)

# Step 2: Generate outcome
outcome_data <- getY_multi_exposure(
  sim_data,
  XYmodels = c("2", "8"),     # Linear effect for exposure 1, quadratic for exposure 2
  X_effects = c(TRUE, TRUE),
  outcome_type = "continuous"
)

# Step 3: Functional PCA for each exposure
fpca_results <- lapply(sim_data$exposures, function(exp_k) {
  FPCA(
    exp_k$Ly_sim,
    exp_k$Lt_sim,
    list(dataType = 'Sparse', error = TRUE, verbose = FALSE)
  )
})

# Step 4: Joint estimation with MV-FMR
result <- mvfmr(
  G = sim_data$details$G,
  fpca_results = fpca_results,
  Y = outcome_data$Y,
  outcome_type = "continuous",
  method = "gmm",
  max_nPC = c(10, 10),
  bootstrap = TRUE,
  n_bootstrap = 100
)

# View results
print(result)
summary(result)
plot(result)  # Visualize time-varying effects for every exposure

# Extract coefficients and effects
coef(result)
result$effects[[1]]  # Time-varying effect for exposure 1
result$effects[[2]]  # Time-varying effect for exposure 2
```

### Example 2: Separate estimation of multiple exposures

Continuing from Example 1: compare joint vs. separate (univariable) estimation
for the same two exposures. For `mvfmr_separate()`, instruments are passed as
`G_list`, a list of length `m` (here the same shared instrument matrix is
reused for both exposures):

```r
# Separate estimation (U-FMR for each exposure independently), reusing
# sim_data / outcome_data / fpca_results simulated in Example 1
result_separate <- mvfmr_separate(
  G_list = list(sim_data$details$G, sim_data$details$G),
  fpca_results = fpca_results,
  Y = outcome_data$Y
)

# Compare performance: joint (`result`, from Example 1) vs. separate
result$performance
result_separate$exposures[[1]]$performance
result_separate$exposures[[2]]$performance
```

### Example 3: Univariable Functional MR (one exposure)

The package can also be used for single exposure analysis (U-FMR), by simulating and passing a single exposure (`n_exposures = 1`, `G_list` of length 1):

```r
library(mvfmr)
library(fdapace)

# Step 1: Simulate a single exposure
set.seed(12345)
sim_data <- getX_multi_exposure(
  N = 1000,
  J = 50,
  nSparse = 10,
  n_exposures = 1
)

# Step 2: Generate outcome
outcome_data <- getY_multi_exposure(
  sim_data,
  XYmodels = "2",     # Linear effect
  X_effects = TRUE,
  outcome_type = "continuous"
)

# Step 3: FPCA for the (only) exposure
fpca1 <- FPCA(
  sim_data$exposures[[1]]$Ly_sim,
  sim_data$exposures[[1]]$Lt_sim,
  list(dataType = 'Sparse', error = TRUE, verbose = FALSE)
)

# Step 4: Univariable estimation
result <- mvfmr_separate(
  G_list = list(sim_data$details$G),  # A list of length 1
  fpca_results = list(fpca1),
  Y = outcome_data$Y,
  outcome_type = "continuous",
  method = "gmm",
  max_nPC = 10
)

# View results
print(result)
coef(result, exposure = 1)
result$exposures[[1]]$effect  # Time-varying effect
```

### Example 4: Extending to more than two exposures (m = 3)

Nothing changes in the API when moving from 2 to `m` exposures: `fpca_results`, `max_nPC`, `true_effects` and `X_true` simply grow to length `m`.

```r
set.seed(2026)
sim_data3 <- getX_multi_exposure(N = 1000, J = 50, nSparse = 10, n_exposures = 3)

outcome_data3 <- getY_multi_exposure(
  sim_data3,
  XYmodels = c("2", "5", "8"),
  outcome_type = "continuous"
)

fpca_results3 <- lapply(sim_data3$exposures, function(exp_k) {
  FPCA(exp_k$Ly_sim, exp_k$Lt_sim, list(dataType = 'Sparse', error = TRUE, verbose = FALSE))
})

# Joint estimation across all 3 exposures
result_joint3 <- mvfmr(
  G = sim_data3$details$G,
  fpca_results = fpca_results3,
  Y = outcome_data3$Y,
  outcome_type = "continuous",
  method = "gmm",
  max_nPC = c(10, 10, 10),
  true_effects = c("2", "5", "8"),
  X_true = sim_data3$details$X_list
)

print(result_joint3)
plot(result_joint3)  # One panel per exposure

# Separate estimation across all 3 exposures
result_separate3 <- mvfmr_separate(
  G_list = list(sim_data3$details$G, sim_data3$details$G, sim_data3$details$G),
  fpca_results = fpca_results3,
  Y = outcome_data3$Y,
  max_nPC = c(10, 10, 10),
  true_effects = c("2", "5", "8")
)

# Access any exposure by index (1..m)
result_joint3$effects[[3]]
coef(result_separate3, exposure = 3)
```

### Example 5: Multivariable two-sample Functional MR

Use outcome GWAS summary statistics instead of individual-level outcome data.

```r
library(mvfmr)
library(fdapace)

# Step 1: Simulate exposure data (individual-level)
set.seed(12345)
sim_data <- getX_multi_exposure(
  N = 5000,           # Exposure sample size
  J = 30,             # Number of genetic instruments (SNPs)
  nSparse = 10,
  n_exposures = 2
)

# Perform FPCA on longitudinal exposures
fpca_results <- lapply(sim_data$exposures, function(exp_k) {
  FPCA(exp_k$Ly_sim, exp_k$Lt_sim, list(dataType = 'Sparse', error = TRUE, verbose = FALSE))
})

# Step 2: Get outcome GWAS summary statistics (from a separate study)
# Simulate obtaining summary statistics from a separate GWAS
# (this mimics what you'd get from a published GWAS)
by_outcome <- rnorm(30, mean = 0.02, sd = 0.01)  # SNP-outcome associations
sy_outcome <- runif(30, 0.005, 0.015)            # Standard errors
ny_outcome <- 100000                             # GWAS sample size

# Step 3: Two-sample MV-FMR estimation
result_twosample <- fmvmr_twosample(
  G_exposure = sim_data$details$G,   # Genotypes from the exposure sample
  fpca_results = fpca_results,       # FPCA from the exposures
  by_outcome = by_outcome,           # GWAS betas (from the outcome study)
  sy_outcome = sy_outcome,           # GWAS standard errors
  ny_outcome = ny_outcome,           # GWAS sample size
  max_nPC = c(3, 3),
  verbose = TRUE
)

# Step 4: View results
print(result_twosample)

# Extract time-varying effects
result_twosample$effects[[1]]
result_twosample$effects[[2]]
```

## Main Functions

### Data Simulation

**`getX_multi_exposure()`** - Generate genetic instruments and exposure data for `m` exposures
```r
getX_multi_exposure(
  N = 1000,                  # Sample size
  J = 50,                    # Number of genetic instruments
  nSparse = 10,               # Observations per subject
  n_exposures = 2,            # Number of exposures (m)
  shared_effect = TRUE,       # Whether all exposures share the same time-varying confounding
  separate_G = FALSE,         # Whether to use separate instruments per exposure
  shared_G_proportion = 0.15  # Proportion of shared instruments (0-1, if separate_G = TRUE)
)
```

**`getX_multi_exposure_mediation()`** - Generate data with mediation pathways between exposures
```r
getX_multi_exposure_mediation(
  N = 1000,                  # Sample size
  J = 50,                    # Number of genetic instruments
  nSparse = 10,               # Observations per subject
  n_exposures = 2,            # Number of exposures (m)
  mediation_strength = NULL,  # m x m matrix: entry [j, k] (j < k) is the strength
                              # with which exposure j mediates its effect onto
                              # exposure k. Default: NULL = no mediation.
  mediation_type = "linear"   # "linear", "nonlinear", "time_varying" (scalar or
                              # m x m matrix mirroring mediation_strength)
)
```

**`getY_multi_exposure()`** - Generate outcome with time-varying effects
```r
getY_multi_exposure(
  RES,                         # Output from getX_multi_exposure() or getX_multi_exposure_mediation()
  XYmodels = NULL,             # Length-m vector of effect models, one per exposure (see below); default '1' for all
  X_effects = NULL,            # Length-m logical vector: include each exposure's effect?; default TRUE for all
  outcome_type = "continuous"  # "continuous" or "binary"
)
```

**Available effect models:**
- `"0"` - No effect (null)
- `"1"` - Constant effect (β = 0.1)
- `"2"` - Linear increasing (β(t) = 0.02×t)
- `"3"` - Linear decreasing (β(t) = 0.5 - 0.02×t)
- `"4"` - Early life effect (β(t) = 0.1 for t < 20)
- `"5"` - Late life effect (β(t) = 0.1 for t > 30)
- `"6"` - Early decreasing (β(t) = 0.05×(20-t) for t < 20)
- `"7"` - Late increasing (β(t) = 0.05×(t-30) for t > 30)
- `"8"` - Quadratic (β(t) = 0.002×t² - 0.11×t + 0.5)
- `"9"` - Cubic (β(t) = -0.00002×t³ + 0.004×t² - 0.2×t + 1)

### Estimation Functions

**`mvfmr()`** - Joint multivariable estimation
```r
mvfmr(
  G,                                     # Genetic instrument matrix (N x J)
  fpca_results,                          # List of length m of FPCA objects, one per exposure
  Y,                                     # Outcome vector
  outcome_type = "continuous",           # "continuous" or "binary"
  method = "gmm",                        # "gmm", "cf" (control function), or "cf-lasso"
  nPC = NA,                              # Fixed number of components per exposure (length 1 or m; NA = select automatically)
  max_nPC = NA,                          # Maximum number of components per exposure (length 1 or m)
  improvement_threshold = 0.001,         # Minimum CV improvement required to add a component
  bootstrap = FALSE,                     # Whether to compute bootstrap confidence intervals
  n_bootstrap = 100,                     # Number of bootstrap replicates
  n_cores = parallel::detectCores() - 1, # Number of CPU cores for parallel computations
  true_effects = NULL,                   # Length-m vector of true effect model codes (simulation only)
  X_true = NULL,                         # Length-m list of true X curves (simulation only)
  verbose = FALSE                        # Print progress and diagnostic messages
)
```

**`mvfmr_separate()`** - Separate univariable estimation
```r
mvfmr_separate(
  G_list,                # List of length m of genetic instrument matrices, one per exposure
                         # (use a list of length 1 to analyze a single exposure)
  fpca_results,          # List of length m of FPCA objects, same length as G_list
  Y,                     # Outcome vector
  outcome_type = "continuous",
  method = "gmm",
  nPC = NA,
  max_nPC = NA,
  improvement_threshold = 0.001,
  bootstrap = FALSE,
  n_bootstrap = 100,
  n_cores = parallel::detectCores() - 1,
  true_effects = NULL,
  X_true = NULL,
  verbose = FALSE
)
```

**`fmvmr_twosample()`** - Two-sample joint multivariable estimation
```r
fmvmr_twosample(
  G_exposure,            # Genetic instrument matrix from the exposure sample (N x J)
  fpca_results,          # List of length m of FPCA objects
  by_outcome,            # Vector of SNP-outcome betas from the outcome GWAS, length J
  sy_outcome,            # Vector of standard errors for SNP-outcome effects, length J
  ny_outcome,            # Sample size of the outcome GWAS
  max_nPC = NA,          # Maximum number of components per exposure (length 1 or m)
  true_effects = NULL,   # Length-m vector of true effect model codes (simulation only)
  verbose = TRUE
)
```

**`fmvmr_separate_twosample()`** - Two-sample separate univariable estimation
```r
fmvmr_separate_twosample(
  G_list,                # List of length m of genetic instrument matrices
  fpca_results,          # List of length m of FPCA objects
  by_outcome_list,       # List of length m of SNP-outcome beta vectors
  sy_outcome_list,       # List of length m of SNP-outcome standard error vectors
  ny_outcome,            # Outcome GWAS sample size
  max_nPC = NA,
  true_effects = NULL,
  verbose = TRUE
)
```

### Utility Functions

**`IS()`** - Calculate instrument strength (F-statistics)
```r
IS(
  J,                     # Number of genetic instruments
  K,                     # Number of exposures/components
  PC,                    # Vector of indices indicating which columns in datafull correspond to the principal components
  datafull,              # Data frame containing instruments (first J columns) and principal components (subsequent columns) [G, X]
  Y                      # Optional outcome vector; if provided, Q-statistic for overidentification is calculated
)
```

## Methods

The package supports three estimation methods:

1. **GMM (Generalized Method of Moments)** - For continuous outcomes
   - Efficient two-step GMM estimation
   - Optimal weighting matrix

2. **Control Function (CF)** - For binary outcomes
   - Two-stage residual inclusion (2SRI)
   - Logistic regression second stage

3. **Control Function (CF)-LASSO** - Control function with LASSO regularization
   - Cross-validated penalty selection

## Output Objects

### `mvfmr` object (from `mvfmr()`)

```r
result <- mvfmr(...)
names(result)
```

Components:
- `coefficients` - Estimated β coefficients for basis functions (stacked across all exposures)
- `vcov` - Variance-covariance matrix
- `effects` - List of length m, one time-varying effect curve per exposure
- `confidence_intervals` - `lower`/`upper`, each a list of length m
- `nPC_used` - Vector of length m: components selected per exposure
- `performance` - `MISE` and `Coverage` (lists of length m), only for simulations
- `plots` - `effects` (list of m ggplot2 objects) and `plot_beta` (combined coefficient plot)

Methods:
- `print()`, `summary()` - Display results
- `plot()` - Visualize time-varying effects for every exposure
- `coef()` - Extract coefficients
- `vcov()` - Extract variance-covariance matrix

### `mvfmr_separate` object (from `mvfmr_separate()`)

```r
result <- mvfmr_separate(...)
names(result)
```

Components:
- `exposures` - List of length m; each entry has `coefficients`, `vcov`, `effect`, `nPC_used`, `performance`
- `plots` - `effects`, a list of m ggplot2 objects

Methods:
- `coef(result, exposure = k)` - Extract coefficients for exposure `k` (1..m)
- `vcov(result, exposure = k)` - Extract variance-covariance matrix for exposure `k`

## Binary Outcomes

For binary outcomes, use `method = "cf"` or `method = "cf-lasso"`:

```r
# Generate binary outcome
outcome_binary <- getY_multi_exposure(
  sim_data,
  XYmodels = c("2", "8"),
  outcome_type = "binary"
)

# Estimate with control function
result <- mvfmr(
  G = sim_data$details$G,
  fpca_results = list(fpca1, fpca2),
  Y = outcome_binary$Y,
  outcome_type = "binary",
  method = "cf"
)
```

## Advanced Features

### Component Selection

Automatic selection via cross-validation:
```r
result <- mvfmr(
  G = G,
  fpca_results = list(fpca1, fpca2),
  Y = Y,
  max_nPC = c(10, 10),          # Search up to 10 components per exposure
  improvement_threshold = 0.01  # Stop if improvement < 1%
)

# View selected components
result$nPC_used
```

### Bootstrap Inference

```r
result <- mvfmr(
  G = G,
  fpca_results = list(fpca1, fpca2),
  Y = Y,
  bootstrap = TRUE,
  n_bootstrap = 200  # Number of bootstrap replicates
)

# Bootstrap confidence intervals available in:
result$confidence_intervals
```

### Parallel Processing

```r
result <- mvfmr(
  G = G,
  fpca_results = list(fpca1, fpca2),
  Y = Y,
  n_cores = 4  # Use 4 cores for cross-validation
)
```

### Mediation Analysis

`mediation_strength` is an m x m matrix: entry `[j, k]` (with `j < k`) is the
strength with which exposure `j` mediates its effect onto exposure `k`. Any
exposure can mediate onto any later one, each with its own strength, so
mediation chains with more than two exposures (e.g. X1 -> X2, X1 -> X3, X2 -> X3)
are supported directly.

```r
# Generate data where exposure 1 mediates onto exposure 2
mediation_strength <- matrix(0, 2, 2)
mediation_strength[1, 2] <- 0.5

sim_mediation <- getX_multi_exposure_mediation(
  N = 1000,
  J = 50,
  n_exposures = 2,
  mediation_strength = mediation_strength,
  mediation_type = "linear"
)

outcome <- getY_multi_exposure(
  sim_mediation,
  XYmodels = c("2", "1"),  # Direct effect of exposure 1; effect of exposure 2 (mediator)
  outcome_type = "continuous"
)

fpca_results <- lapply(sim_mediation$exposures, function(exp_k) {
  FPCA(exp_k$Ly_sim, exp_k$Lt_sim, list(dataType = 'Sparse', error = TRUE, verbose = FALSE))
})

# Estimate with MV-FMR to capture mediation
result <- mvfmr(
  G = sim_mediation$details$G,
  fpca_results = fpca_results,
  Y = outcome$Y
)
```

## Instrument Strength Diagnostics

Check instrument strength with F-statistics (`IS()` is generic in the number of exposures/components `K`):

```r
# After FPCA
K_total <- sum(sapply(fpca_results, function(f) f$selectK))

PC_stacked <- do.call(cbind, lapply(fpca_results, function(f) f$xiEst[, 1:f$selectK]))

fstats <- IS(
  J = ncol(G),
  K = K_total,
  PC = 1:K_total,
  datafull = cbind(G, PC_stacked)
)

# View conditional F-statistics (cFF)
print(fstats)
```

## Performance Metrics

When true effects are provided (simulations):

- **MISE (Mean Integrated Squared Error)**: Average squared difference between estimated and true effect curves
- **Coverage**: Proportion of time points where true effect falls within 95% CI


## Acknowledgments and Related Work

This package extends the **univariable functional Mendelian Randomization** framework to the multivariable setting. Key related work:

### TVMR Package (Univariable Functional MR)
Our implementation builds upon and extends the TVMR package by Tian et al.:

> **Tian, H., Mason, A. M., Liu, C., & Burgess, S.** (2024). Estimating time‐varying exposure effects through continuous‐time modelling in Mendelian randomization. *Statistics in Medicine*, 43(26), 5006-5024. https://doi.org/10.1002/sim.10222

**GitHub**: https://github.com/HDTian/TVMR

## Author

Nicole Fontana

## License

MIT — see the [LICENSE](LICENSE) file for details.

## Getting Help

For questions and issues:
- Open an issue on GitHub
- Email: nicole.fontana@polimi.it
