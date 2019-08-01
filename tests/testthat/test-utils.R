context("test-utils")

test_that("seq_log() works", {
  expect_equal(seq_log(1, 1000, 4), 10^(0:3))
  expect_equal(seq_log(1, 100,  5), 10^(0:4 / 2))
  expect_equal(seq_log(1000, 1, 4), rev(seq_log(1, 1000, 4)))
  expect_equal(seq_log(100,  1, 5), rev(seq_log(1, 100,  5)))
  expect_equal(seq_log(1, 1, 1), 1)
  expect_equal(seq_log(1, 1, 5), rep(1, 5))
  expect_error(seq_log(1, 1000, -4), "'length.out' must be a non-negative number")
})
