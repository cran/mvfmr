make_twosample_fixture <- function(m, N = 200, J = 10, nSparse = 6, seed = 1) {
  fx <- make_fixture(m, N = N, J = J, nSparse = nSparse, seed = seed)
  by_outcome <- sapply(seq_len(J), function(j) {
    stats::coef(stats::lm(fx$dat$Y ~ fx$sim$details$G[, j]))[2]
  })
  sy_outcome <- sapply(seq_len(J), function(j) {
    summary(stats::lm(fx$dat$Y ~ fx$sim$details$G[, j]))$coefficients[2, 2]
  })
  c(fx, list(by_outcome = by_outcome, sy_outcome = sy_outcome))
}

test_that("fmvmr_twosample() works for m = 2 and m = 3", {
  for (m in 2:3) {
    fx <- make_twosample_fixture(m, seed = 300 + m)

    result <- fmvmr_twosample(
      G_exposure = fx$sim$details$G,
      fpca_results = fx$fpca_results,
      by_outcome = fx$by_outcome,
      sy_outcome = fx$sy_outcome,
      ny_outcome = fx$N,
      max_nPC = rep(2, m),
      true_effects = fx$models,
      verbose = FALSE
    )

    expect_s3_class(result, "fmvmr_twosample")
    expect_s3_class(result, "fmvmr")
    expect_equal(result$n_exposures, m)
    expect_equal(length(result$effects), m)
    expect_equal(length(result$nPC_used), m)
    expect_true(is.numeric(result$Q_stat))
  }
})

test_that("fmvmr_separate_twosample() works for m = 1..3", {
  for (m in 1:3) {
    fx <- make_twosample_fixture(m, seed = 400 + m)

    result <- fmvmr_separate_twosample(
      G_list = lapply(seq_len(m), function(k) fx$sim$details$G),
      fpca_results = fx$fpca_results,
      by_outcome_list = lapply(seq_len(m), function(k) fx$by_outcome),
      sy_outcome_list = lapply(seq_len(m), function(k) fx$sy_outcome),
      ny_outcome = fx$N,
      max_nPC = rep(2, m),
      true_effects = fx$models,
      verbose = FALSE
    )

    expect_s3_class(result, "fmvmr_separate_twosample")
    expect_s3_class(result, "fmvmr_separate")
    expect_equal(result$n_exposures, m)
    expect_equal(length(result$exposures), m)
  }
})
