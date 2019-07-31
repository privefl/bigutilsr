context("test-rollmean")

test_that("rollmean() works", {

  x <- rnorm(100)
  expect_equal(rollmean(x, 0), x)

  x.roll <- sapply(c(0, 1, 3, 5, 10, 20), function(size) rollmean(x, size))
  expect_false(is.unsorted(rev(apply(x.roll, 2, var))))

  expect_error(rollmean(x, 50), "Parameter 'size' is too large.")
  expect_error(rollmean(x, -1), "Parameter 'size' must be positive.")

  expect_equal(rollmean(rep(1, 100), 20), rep(1, 100))
  expect_equal(rollmean(rep(c(0, 1), 100), 50), rep(0.5, 200), tolerance = 1e-2)
})
