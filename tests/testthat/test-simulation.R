test_that("getX_multi_exposure generates m exposures", {
  for (m in 1:4) {
    set.seed(m)
    sim <- getX_multi_exposure(N = 50, J = 8, nSparse = 5, n_exposures = m)

    expect_equal(length(sim$exposures), m)
    expect_equal(sim$details$n_exposures, m)
    expect_equal(length(sim$details$X_list), m)
    expect_equal(nrow(sim$details$G), 50)

    for (k in seq_len(m)) {
      expect_equal(length(sim$exposures[[k]]$Ly_sim), 50)
      expect_equal(nrow(sim$details$X_list[[k]]), 50)
    }
  }
})

test_that("getX_multi_exposure supports separate_G with shared proportion", {
  sim <- getX_multi_exposure(N = 40, J = 8, nSparse = 5, n_exposures = 3,
                              separate_G = TRUE, shared_G_proportion = 0.25)

  expect_equal(length(sim$details$G_list), 3)
  for (k in 1:3) {
    expect_equal(nrow(sim$details$G_list[[k]]), 40)
  }
})

test_that("getY_multi_exposure produces an outcome vector of correct length for m exposures", {
  for (m in 1:4) {
    fx <- make_fixture(m, N = 40, J = 8, nSparse = 5, seed = m)
    expect_equal(length(fx$dat$Y), 40)
    expect_true(is.numeric(fx$dat$Y))
  }
})

test_that("getY_multi_exposure supports binary outcomes", {
  fx <- make_fixture(2, N = 60, J = 8, nSparse = 5, seed = 42)
  dat_bin <- getY_multi_exposure(fx$sim, XYmodels = fx$models, outcome_type = "binary")
  expect_true(all(dat_bin$Y %in% c(0, 1)))
})

test_that("getX_multi_exposure_mediation accepts a valid upper-triangular mediation matrix", {
  m <- 3
  mediation_strength <- matrix(0, m, m)
  mediation_strength[1, 2] <- 0.2
  mediation_strength[1, 3] <- 0.1
  mediation_strength[2, 3] <- 0.3

  set.seed(7)
  sim <- getX_multi_exposure_mediation(N = 40, J = 8, nSparse = 5, n_exposures = m,
                                        mediation_strength = mediation_strength)

  expect_equal(length(sim$exposures), m)
  expect_equal(sim$details$n_exposures, m)
})

test_that("getX_multi_exposure_mediation rejects invalid (non-upper-triangular) mediation matrices", {
  m <- 2
  bad_strength <- matrix(0, m, m)
  bad_strength[2, 1] <- 0.3  # invalid: exposure 2 cannot mediate onto exposure 1

  expect_error(
    getX_multi_exposure_mediation(N = 20, J = 5, nSparse = 4, n_exposures = m,
                                   mediation_strength = bad_strength)
  )
})

test_that("getX_multi_exposure_mediation defaults to no mediation", {
  set.seed(3)
  sim <- getX_multi_exposure_mediation(N = 30, J = 6, nSparse = 4, n_exposures = 2)
  expect_true(all(sim$details$mediation_strength == 0))
})

test_that("getX_multi_exposure_mediation correctly combines multiple mediating parents (m = 3)", {
  # Exposure 3 is mediated by BOTH exposure 1 and exposure 2 with different
  # strengths; the resulting mediation_effect for exposure 3 must equal the
  # exact linear combination of the two parents' contributions.
  m <- 3
  strength <- matrix(0, m, m)
  strength[1, 3] <- 0.05
  strength[2, 3] <- 0.9

  set.seed(123)
  sim <- getX_multi_exposure_mediation(N = 50, J = 6, nSparse = 5, n_exposures = m,
                                        mediation_strength = strength)

  X1 <- sim$details$X_list[[1]]
  X2 <- sim$details$X_list[[2]]
  expected_me3 <- strength[1, 3] * X1 + strength[2, 3] * X2

  expect_equal(sim$details$mediation_effect[[3]], expected_me3)

  # exposure 2 has no mediating parents here, so its mediation effect is zero
  expect_true(all(sim$details$mediation_effect[[2]] == 0))
})

test_that("mvfmr() and mvfmr_separate() run end-to-end on m = 3 mediated data", {
  m <- 3
  strength <- matrix(0, m, m)
  strength[1, 2] <- 0.2
  strength[1, 3] <- 0.1
  strength[2, 3] <- 0.3

  set.seed(321)
  sim <- getX_multi_exposure_mediation(N = 60, J = 8, nSparse = 6, n_exposures = m,
                                        mediation_strength = strength)
  fpca_results <- lapply(sim$exposures, function(exp_k) {
    fdapace::FPCA(exp_k$Ly_sim, exp_k$Lt_sim,
                  list(dataType = 'Sparse', error = TRUE, verbose = FALSE))
  })
  dat <- getY_multi_exposure(sim, XYmodels = c("2", "8", "4"), outcome_type = "continuous")

  result_joint <- mvfmr(
    G = sim$details$G, fpca_results = fpca_results, Y = dat$Y,
    max_nPC = rep(2, m), n_cores = 1, true_effects = c("2", "8", "4"),
    X_true = sim$details$X_list, verbose = FALSE
  )
  expect_s3_class(result_joint, "mvfmr")
  expect_equal(length(result_joint$effects), m)

  result_sep <- mvfmr_separate(
    G_list = lapply(seq_len(m), function(k) sim$details$G),
    fpca_results = fpca_results, Y = dat$Y,
    max_nPC = rep(2, m), n_cores = 1, true_effects = c("2", "8", "4"),
    verbose = FALSE
  )
  expect_s3_class(result_sep, "mvfmr_separate")
  expect_equal(length(result_sep$exposures), m)
})

