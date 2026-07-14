# Shared small fixtures for m-exposure smoke tests.
# Kept intentionally small (N, J, max_nPC) so the CV component-selection
# loop in AUTOMATIC_Multi_MVFMR does minimal work and the suite stays fast.

make_fixture <- function(m, N = 60, J = 8, nSparse = 5, seed = 1) {
  set.seed(seed)
  sim <- getX_multi_exposure(N = N, J = J, nSparse = nSparse, n_exposures = m)

  fpca_results <- lapply(sim$exposures, function(exp_k) {
    fdapace::FPCA(
      exp_k$Ly_sim, exp_k$Lt_sim,
      list(dataType = 'Sparse', error = TRUE, verbose = FALSE)
    )
  })

  models <- rep(c("2", "8", "4", "5"), length.out = m)
  dat <- getY_multi_exposure(sim, XYmodels = models, outcome_type = "continuous")

  list(sim = sim, fpca_results = fpca_results, dat = dat, models = models, N = N, J = J)
}
