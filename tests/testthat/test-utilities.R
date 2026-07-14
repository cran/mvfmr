test_that("compute_offsets and block_idx produce correct blocks", {
  nPC_vec <- c(3, 4, 2)
  offsets <- compute_offsets(nPC_vec)

  expect_equal(offsets, c(0, 3, 7, 9))
  expect_equal(block_idx(offsets, 1), 1:3)
  expect_equal(block_idx(offsets, 2), 4:7)
  expect_equal(block_idx(offsets, 3), 8:9)
})

test_that("recycle_arg recycles scalars and validates length", {
  expect_equal(recycle_arg(NULL, 3, default = NA), rep(NA, 3))
  expect_equal(recycle_arg(5, 3, default = NA), rep(5, 3))
  expect_equal(recycle_arg(c(1, 2, 3), 3, default = NA), c(1, 2, 3))
  expect_error(recycle_arg(c(1, 2), 3, default = NA))
})
