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

- **Multivariable Functional MR (MV-FMR)**: Joint estimation of multiple time-varying exposures
- **Univariable Functional MR (U-FMR)**: Estimation of single time-varying exposures
- **Both continuous and binary outcomes**
- **Automatic component selection** via cross-validation
- **Bootstrap inference** for confidence intervals
- **Instrument strength diagnostics**

## Installation

Install the package:

```r
install.packages("devtools") # Install devtools if not already installed
devtools::install_github("NicoleFontana/mvfmr") # Install mvfmr from GitHub
```

Required dependencies:
```r
install.packages(c("fdapace", "ggplot2", "glmnet", "pROC", "parallel", "doParallel", "foreach", "progress", "dplyr", "gridExtra"))
```

## Test Scripts and Simulations

The package includes comprehensive test scripts demonstrating different use cases:

### 1. **Manuscript Simulations** (`tests_manuscript.R`)
Reproduces main simulation scenarios from the manuscript:
- **Scenario 1-3**: pleiotropy, null effects, mediation
- **Exposure effects**: linear and quadratic
- **Performance comparison**: MV-FMR vs U-FMR across scenarios
- **Evaluation**: MISE, coverage rates

```r
# Download and run manuscript simulations
source("https://raw.githubusercontent.com/NicoleFontana/mvfmr/demo/tests_manuscript.R")
```

### 2. **Multivariable FMR Tutorial** (`test_MV-FMR.R`)
Complete tutorial for using **joint multivariable** estimation:
- Data simulation with two exposures
- FPCA and component selection
- Joint estimation with `fmvmr()`
- Instrument diagnostics
- Performance metrics and visualization
- Binary outcome analysis
- Comparison with univariable estimation

```r
# Download and explore MV-FMR tutorial
source("https://raw.githubusercontent.com/NicoleFontana/mvfmr/demo/test_MV-FMR.R")
```

### 3. **Univariable FMR Tutorial** (`test_U-FMR.R`)
Complete tutorial for **single exposure** analysis:
- Single exposure simulation
- Univariable estimation with `fmvmr_separate()`
- Instrument diagnostics
- Performance metrics and visualization
- Comparison of different exposure effects
- Binary outcome analysis

```r
# Download and explore U-FMR tutorial
source("https://raw.githubusercontent.com/NicoleFontana/mvfmr/demo/test_U-FMR.R")
```

**Note**: These scripts serve as templates for your own analyses. Download them, modify the parameters, and adapt to your specific research questions.


## Quick Start

### Example 1: Multivariable Functional MR (two exposures)

```r
library(fdapace)

# Step 1: Simulate exposures data
set.seed(12345)
sim_data <- getX_multi_exposure(
  N = 1000,           # Sample size
  J = 50,             # Number of genetic instruments
  nSparse = 10       # Sparse observations per subject
)

# Step 2: Generate outcome 
outcome_data <- getY_multi_exposure(
  sim_data,
  X1Ymodel = "2",     # Linear effect for exposure 1
  X2Ymodel = "8",     # Quadratic effect for exposure 2
  X1_effect = TRUE,
  X2_effect = TRUE,
  outcome_type = "continuous"
)

# Step 3: Functional PCA for both the exposures
fpca1 <- FPCA(
  sim_data$X1$Ly_sim, 
  sim_data$X1$Lt_sim,
  list(dataType = 'Sparse', error = TRUE, verbose = FALSE)
)

fpca2 <- FPCA(
  sim_data$X2$Ly_sim, 
  sim_data$X2$Lt_sim,
  list(dataType = 'Sparse', error = TRUE, verbose = FALSE)
)

# Step 4: Joint estimation with MV-FMR
result <- fmvmr(
  G = sim_data$details$G,
  fpca_results = list(fpca1, fpca2),
  Y = outcome_data$Y,
  outcome_type = "continuous",
  method = "gmm",
  max_nPC1 = 10,
  max_nPC2 = 10,
  bootstrap = TRUE,
  n_bootstrap = 100
)

# View results
print(result)
summary(result)
plot(result)  # Visualize time-varying effects

# Extract coefficients and effects
coef(result)
result$effects$effect1  # Time-varying effect for exposure 1
result$effects$effect2  # Time-varying effect for exposure 2
```

### Example 2: Univariable Functional MR (one exposure)

The package can also be used for single exposure analysis (U-FMR):

```r
library(fdapace)

# Step 1: Simulate exposure data
set.seed(12345)
sim_data <- getX_multi_exposure(
  N = 1000,
  J = 50,
  nSparse = 10,
  shared_effect = FALSE
)

# Step 2: Generate outcome (only exposure 1 has effect)
outcome_data <- getY_multi_exposure(
  sim_data,
  X1Ymodel = "2",     # Linear effect for exposure 1
  X2Ymodel = "0",     # No effect for exposure 2
  X1_effect = TRUE,
  X2_effect = FALSE,
  outcome_type = "continuous"
)

# Step 3: FPCA for exposure 1 only
fpca1 <- FPCA(
  sim_data$X1$Ly_sim, 
  sim_data$X1$Lt_sim,
  list(dataType = 'Sparse', error = TRUE, verbose = FALSE)
)

# Step 4: Univariable estimation for exposure 1 only
result <- fmvmr_separate(
  G1 = sim_data$details$G,  # Instruments for exposure 1
  G2 = NULL,                # No instruments for exposure 2
  fpca_results = list(fpca1),
  Y = outcome_data$Y,
  outcome_type = "continuous",
  method = "gmm",
  max_nPC1 = 10,
  max_nPC2 = 10
)

# View results for exposure 1
print(result)
coef(result, exposure = 1)
result$exposure1$effect  # Time-varying effect
```

### Example 3: Univariable separate estimation of multiple exposures

Compare joint vs. univariable separate estimation:

```r
# Joint estimation (MV-FMR)
result_joint <- fmvmr(
  G = sim_data$details$G,
  fpca_results = list(fpca1, fpca2),
  Y = outcome_data$Y
)

# Separate estimation (U-FMR for each exposure independently)
result_separate <- fmvmr_separate(
  G1 = sim_data$details$G,
  G2 = sim_data$details$G,
  fpca_results = list(fpca1, fpca2),
  Y = outcome_data$Y
)

# Compare performance
result_joint$performance
result_separate$exposure1$performance
result_separate$exposure2$performance
```

### Example 4: Multvariable two-sample Functional MR

Use outcome GWAS summary statistics instead of individual-level outcome data. 

```r
library(fdapace)

# Step 1: Simulate exposure data (individual-level)
set.seed(12345)
sim_data <- getX_multi_exposure(
  N = 5000,           # Exposure sample size
  J = 30,             # Number of genetic instruments (SNPs)
  nSparse = 10
)

# Perform FPCA on longitudinal exposures
fpca1 <- FPCA(
  sim_data$X1$Ly_sim, 
  sim_data$X1$Lt_sim,
  list(dataType = 'Sparse', error = TRUE, verbose = FALSE)
)

fpca2 <- FPCA(
  sim_data$X2$Ly_sim, 
  sim_data$X2$Lt_sim,
  list(dataType = 'Sparse', error = TRUE, verbose = FALSE)
)

# Step 2: Get outcome GWAS summary statistics (from separate study)
# Load by_outcome, sy_outcome, ny_outcome from a public GWAS

# Simulate obtaining summary statistics from a separate GWAS
# (This mimics what you'd get from a published GWAS)
by_outcome <- rnorm(30, mean = 0.02, sd = 0.01)  # SNP-outcome associations
sy_outcome <- runif(30, 0.005, 0.015)            # Standard errors
ny_outcome <- 100000                             # GWAS sample size

# Step 3: Two-sample MV-FMR estimation
result_twosample <- fmvmr_twosample(
  G_exposure = sim_data$details$G,        # Genotypes from exposure sample
  fpca_results = list(fpca1, fpca2),      # FPCA from exposures
  by_outcome = by_outcome,                # GWAS betas (from outcome study)
  sy_outcome = sy_outcome,                # GWAS standard errors
  ny_outcome = ny_outcome,                # GWAS sample size
  max_nPC1 = 3,
  max_nPC2 = 3,
  verbose = TRUE
)

# Step 4: View results
print(result_twosample)

# Extract time-varying effects
result_twosample$effects$effect1 
result_twosample$effects$effect2  

```


## Main Functions

### Data Simulation

**`getX_multi_exposure()`** - Generate genetic instruments and exposure data
```r
getX_multi_exposure(
  N = 1000,                  # Sample size
  J = 50,                    # Number of genetic instruments
  nSparse = 10,              # Observations per subject
  shared_G_proportion = 0.15 # Proportion of shared instruments (0-1)
)
```

**`getX_multi_exposure_mediation()`** - Generate data with mediation
```r
getX_multi_exposure_mediation(
  N = 1000,                  # Sample size
  J = 50,                    # Number of genetic instruments
  nSparse = 10,              # Observations per subject
  mediation_strength = 0.3,  # Strength of X1 → X2
  mediation_type = "linear"  # "linear", "nonlinear", "time_varying"
)
```

**`getY_multi_exposure()`** - Generate outcome with time-varying effects
```r
getY_multi_exposure(
  RES,                         # Output from getX_multi_exposure()
  X1Ymodel = "2",              # Effect model for X1 (see below)
  X2Ymodel = "8",              # Effect model for X2
  X1_effect = TRUE,            # X1 has an effect on Y
  X2_effect = TRUE,            # X2 has an effect on Y
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

**`fmvmr()`** - Joint multivariable estimation
```r
fmvmr(
  G,                                     # Genetic instrument matrix (N × J)
  fpca_results,                          # List of 2 FPCA objects
  Y,                                     # Outcome vector
  outcome_type = "continuous",           # Type of outcome: "continuous" for numeric outcomes, "binary" for 0/1 outcomes
  method = "gmm",                        # Estimation method: "gmm" (Generalized Method of Moments), "cf" (control function), or "cf-lasso" (control function with Lasso)
  nPC1_selected = NA,                    # Fixed number of principal components to retain for exposure 1 (NA = select automatically)
  max_nPC1 = NA,                         # Maximum number of principal components to retain for exposure 1 (NA = automatically determined)
  nPC2_selected = NA,                    # Fixed number of principal components to retain for exposure 2 (NA = select automatically)
  max_nPC2 = NA,                         # Maximum number of principal components to retain for exposure 2 (NA = automatically determined)
  improvement_threshold = 0.001,         # Minimum cross-validation improvement required to add an additional principal component
  bootstrap = FALSE,                     # Whether to compute confidence intervals using bootstrap resampling
  n_bootstrap = 100,                     # Number of bootstrap replicates (only used if bootstrap = TRUE)
  n_cores = parallel::detectCores() - 1, # Number of CPU cores to use for parallel computations
  verbose = TRUE                         # Print progress and diagnostic messages during computation
)
```

**`fmvmr_separate()`** - Separate univariable estimation
```r
fmvmr_separate(
  G1,                    # Genetic instrument matrix for exposure 1
  G2,                    # Genetic instrument matrix for exposure 2, or NULL if only a single exposure is analyzed
  fpca_results,          # List of FPCA objects
  Y,                     # Outcome vector
  outcome_type = "continuous", # Type of outcome: "continuous" for numeric outcomes, "binary" for 0/1 outcomes
  method = "gmm",        # Estimation method: "gmm" (Generalized Method of Moments), "cf" (control function), or "cf-lasso" (control function with Lasso)
  nPC1_selected = NA,    # Fixed number of principal components to retain for exposure 1 (NA = select automatically)
  max_nPC1 = NA,         # Maximum number of principal components to retain for exposure 1 (NA = automatically determined)
  nPC2_selected = NA,    # Fixed number of principal components to retain for exposure 2 (NA = select automatically)
  max_nPC2 = NA,         # Maximum number of principal components to retain for exposure 2 (NA = automatically determined; ignored if G2 is NULL)
  improvement_threshold = 0.001,         # Minimum cross-validation improvement required to add an additional principal component
  bootstrap = FALSE,     # Whether to compute confidence intervals using bootstrap resampling
  n_bootstrap = 100,     # Number of bootstrap replicates (only used if bootstrap = TRUE)
  n_cores = parallel::detectCores() - 1, # Number of CPU cores to use for parallel computations
  verbose = TRUE                         # Print progress and diagnostic messages during computation
)

```

**`fmvmr_twosample()`** - Two-sample joint multivariable estimation
```r
rfmvmr_twosample(
  G_exposure,            # Genetic instrument matrix from the exposure sample (N × J)
  fpca_results,          # List of 2 FPCA objects corresponding to the two exposures from the exposure data
  by_outcome,            # Vector of SNP-outcome effect estimates (betas) from the outcome GWAS, length J
  sy_outcome,            # Vector of standard errors for SNP-outcome effects, length J
  ny_outcome,            # Sample size of the outcome GWAS
  max_nPC1 = NA,         # Maximum number of principal components to retain for exposure 1 (NA = default)
  max_nPC2 = NA,         # Maximum number of principal components to retain for exposure 2 (NA = default)
  true_effects = NULL,   # For simulation studies: list containing true effects for exposure 1 and exposure 2 (e.g., list(model1, model2))
  verbose = TRUE         # Print progress messages and diagnostics during computation
)
```


**`fmvmr_separate_twosample()`** - Two-sample separate univariable estimation 
```r
rfmvmr_separate_twosample(
  G1_exposure,           # Genetic instrument matrix from exposure 1 (N × J1)
  G2_exposure = NULL,    # Genetic instrument matrix from exposure 2 (N × J2) or NULL for single exposure
  fpca_results,          # List of 2 FPCA objects
  by_outcome1,           # SNP-outcome betas for exposure 1 instruments
  by_outcome2 = NULL,    # SNP-outcome betas for exposure 2 or NULL
  sy_outcome1,           # Standard errors for exposure 1
  sy_outcome2 = NULL,    # Standard errors for exposure 2 or NULL
  ny_outcome,            # Outcome GWAS sample size
  max_nPC1 = NA,         # Maximum number of principal components to retain for exposure 1 (NA = automatically determined)
  max_nPC2 = NA,         # Maximum number of principal components to retain for exposure 2 (NA = automatically determined)
  true_effects = NULL,   # List containing true effects for exposure 1 and exposure 2 (simulation only)
  verbose = TRUE         # Print progress messages and diagnostics during computation
)
```

### Utility Functions

**`IS()`** - Calculate instrument strength (F-statistics)
```r
IS(
  J,                     # Number of genetic instruments
  K,                     # Number of exposures
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

### `fmvmr` object (from `fmvmr()`)

```r
result <- fmvmr(...)
names(result)
```

Components:
- `coefficients` - Estimated β coefficients for basis functions
- `vcov` - Variance-covariance matrix
- `effects` - List with `effect1`, `effect2`, `time_grid`
- `confidence_intervals` - Upper and lower bounds
- `nPC_used` - Selected components (nPC1, nPC2)
- `performance` - MISE and coverage (only for simulations)
- `plots` - ggplot2 objects for visualization

Methods:
- `print()`, `summary()` - Display results
- `plot()` - Visualize time-varying effects
- `coef()` - Extract coefficients
- `vcov()` - Extract variance-covariance matrix

### `fmvmr_separate` object (from `fmvmr_separate()`)

```r
result <- fmvmr_separate(...)
names(result)
```

Components:
- `exposure1` - Results for exposure 1
  - `coefficients`, `vcov`, `effect`, `nPC_used`, `performance`
- `exposure2` - Results for exposure 2 (if provided)
  - `coefficients`, `vcov`, `effect`, `nPC_used`, `performance`
- `plots` - Visualization objects

Methods:
- `coef(result, exposure = 1)` - Extract coefficients for specific exposure
- `vcov(result, exposure = 1)` - Extract variance-covariance matrix

## Binary Outcomes

For binary outcomes, use `method = "cf"` or `method = "cf-lasso"`:

```r
# Generate binary outcome
outcome_binary <- getY_multi_exposure(
  sim_data,
  X1Ymodel = "2",
  X2Ymodel = "8",
  outcome_type = "binary"
)

# Estimate with control function
result <- fmvmr(
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
result <- fmvmr(
  G = G,
  fpca_results = list(fpca1, fpca2),
  Y = Y,
  max_nPC1 = 10,              # Search up to 10 components
  max_nPC2 = 10,
  improvement_threshold = 0.01 # Stop if improvement < 1%
)

# View selected components
result$nPC_used
```

### Bootstrap Inference

```r
result <- fmvmr(
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
result <- fmvmr(
  G = G,
  fpca_results = list(fpca1, fpca2),
  Y = Y,
  n_cores = 4  # Use 4 cores for cross-validation
)
```

### Mediation Analysis

```r
# Generate data with X1 → X2 mediation
sim_mediation <- getX_multi_exposure_mediation(
  N = 1000,
  J = 50,
  mediation_strength = 0.5,
  mediation_type = "linear"
)

outcome <- getY_multi_exposure(
  sim_mediation,
  X1Ymodel = "2",  # Direct effect of X1
  X2Ymodel = "1",  # Effect of X2 (mediator)
  outcome_type = "continuous"
)

# Estimate with MV-FMR to capture mediation
result <- fmvmr(
  G = sim_mediation$details$G,
  fpca_results = list(fpca1, fpca2),
  Y = outcome$Y
)
```

## Instrument Strength Diagnostics

Check instrument strength with F-statistics:

```r
# After FPCA
K_total <- fpca1$selectK + fpca2$selectK

fstats <- IS(
  J = ncol(G),
  K = K_total,
  PC = 1:K_total,
  datafull = cbind(
    G,
    cbind(fpca1$xiEst[, 1:fpca1$selectK], 
          fpca2$xiEst[, 1:fpca2$selectK])
  )
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

[License information]

## Getting Help

For questions and issues:
- Open an issue on GitHub
- Email: nicole.fontana@polimi.it
