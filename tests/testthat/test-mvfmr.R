check_mvfmr_joint <- function(m) {
  fx <- make_fixture(m, seed = 100 + m)

  result <- mvfmr(
    G = fx$sim$details$G,
    fpca_results = fx$fpca_results,
    Y = fx$dat$Y,
    outcome_type = "continuous",
    method = "gmm",
    max_nPC = rep(2, m),
    n_cores = 1,
    true_effects = fx$models,
    X_true = fx$sim$details$X_list,
    verbose = FALSE
  )

  expect_s3_class(result, "mvfmr")
  expect_equal(result$n_exposures, m)
  expect_equal(length(result$effects), m)
  expect_equal(length(result$nPC_used), m)
  expect_equal(length(result$plots$effects), m)
  expect_false(is.null(result$performance))
  expect_equal(length(result$performance$MISE), m)
  expect_equal(length(result$coefficients), sum(result$nPC_used))
  expect_equal(dim(result$vcov), c(sum(result$nPC_used), sum(result$nPC_used)))
}

test_that("mvfmr() works for m = 1..2 exposures", {
  for (m in 1:2) check_mvfmr_joint(m)
})

test_that("mvfmr() works for m = 3..4 exposures", {
  skip_on_cran()
  for (m in 3:4) check_mvfmr_joint(m)
})

test_that("mvfmr() supports binary outcome with control function methods", {
  fx <- make_fixture(2, N = 100, J = 10, nSparse = 6, seed = 5)
  dat_bin <- getY_multi_exposure(fx$sim, XYmodels = fx$models, outcome_type = "binary")

  for (method in c("cf", "cf-lasso")) {
    result <- mvfmr(
      G = fx$sim$details$G,
      fpca_results = fx$fpca_results,
      Y = dat_bin$Y,
      outcome_type = "binary",
      method = method,
      max_nPC = c(2, 2),
      n_cores = 1,
      verbose = FALSE
    )
    expect_s3_class(result, "mvfmr")
    expect_equal(length(result$effects), 2)
  }
})

test_that("mvfmr() bootstrap produces confidence intervals for m = 2", {
  fx <- make_fixture(2, N = 60, J = 8, nSparse = 5, seed = 9)

  result <- mvfmr(
    G = fx$sim$details$G,
    fpca_results = fx$fpca_results,
    Y = fx$dat$Y,
    outcome_type = "continuous",
    method = "gmm",
    max_nPC = c(2, 2),
    bootstrap = TRUE,
    n_bootstrap = 3,
    n_cores = 1,
    verbose = FALSE
  )

  expect_equal(length(result$raw_result$CI_beta_k), 2)
  expect_equal(length(result$raw_result$CI_beta_t), 2)
})

check_mvfmr_separate <- function(m) {
  fx <- make_fixture(m, seed = 200 + m)

  result <- mvfmr_separate(
    G_list = lapply(seq_len(m), function(k) fx$sim$details$G),
    fpca_results = fx$fpca_results,
    Y = fx$dat$Y,
    outcome_type = "continuous",
    method = "gmm",
    max_nPC = rep(2, m),
    n_cores = 1,
    true_effects = fx$models,
    verbose = FALSE
  )

  expect_s3_class(result, "mvfmr_separate")
  expect_equal(result$n_exposures, m)
  expect_equal(length(result$exposures), m)
  expect_equal(length(result$plots$effects), m)

  for (k in seq_len(m)) {
    expect_false(is.null(result$exposures[[k]]$coefficients))
    expect_false(is.null(result$exposures[[k]]$performance))
  }
}

test_that("mvfmr_separate() works for m = 1..2 exposures", {
  for (m in 1:2) check_mvfmr_separate(m)
})

test_that("mvfmr_separate() works for m = 3..4 exposures", {
  skip_on_cran()
  for (m in 3:4) check_mvfmr_separate(m)
})

test_that("coef.mvfmr_separate / vcov.mvfmr_separate validate the exposure index", {
  fx <- make_fixture(3, N = 60, J = 8, nSparse = 5, seed = 11)

  result <- mvfmr_separate(
    G_list = lapply(seq_len(3), function(k) fx$sim$details$G),
    fpca_results = fx$fpca_results,
    Y = fx$dat$Y,
    outcome_type = "continuous",
    method = "gmm",
    max_nPC = rep(2, 3),
    n_cores = 1,
    verbose = FALSE
  )

  expect_equal(coef(result, exposure = 1), result$exposures[[1]]$coefficients)
  expect_equal(coef(result, exposure = 3), result$exposures[[3]]$coefficients)
  expect_equal(vcov(result, exposure = 2), result$exposures[[2]]$vcov)
  expect_error(coef(result, exposure = 4))
  expect_error(vcov(result, exposure = 0))
})

test_that("print/summary/plot S3 methods run without error for m = 3", {
  fx <- make_fixture(3, N = 60, J = 8, nSparse = 5, seed = 13)

  result_joint <- mvfmr(
    G = fx$sim$details$G,
    fpca_results = fx$fpca_results,
    Y = fx$dat$Y,
    outcome_type = "continuous",
    method = "gmm",
    max_nPC = rep(2, 3),
    n_cores = 1,
    verbose = FALSE
  )

  result_sep <- mvfmr_separate(
    G_list = lapply(seq_len(3), function(k) fx$sim$details$G),
    fpca_results = fx$fpca_results,
    Y = fx$dat$Y,
    outcome_type = "continuous",
    method = "gmm",
    max_nPC = rep(2, 3),
    n_cores = 1,
    verbose = FALSE
  )

  expect_output(print(result_joint))
  expect_output(print(result_sep))
  expect_output(summary(result_joint))
  expect_output(summary(result_sep))

  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_silent(plot(result_joint))
  expect_silent(plot(result_sep))
})
