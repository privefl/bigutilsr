context("test-outlier-hist")

test_that("hist_out() works with size 200", {

  set.seed(1)
  x <- rnorm(200)
  expect_equal(hist_out(x)$lim, c(-Inf, Inf))

  # Easy to separate
  x2 <- c(x, rnorm(20, mean = 6))
  hist(x2, breaks = nclass.scottRob)
  expect_equal(hist_out(x2)$lim, c(-Inf, 2.75))

  # More difficult to separate
  x3 <- c(x, rnorm(30, mean = 5, sd = 1))
  hist(x3, breaks = nclass.scottRob)
  expect_equal(hist_out(x3)$lim, c(-Inf, Inf))
  expect_equal(hist_out(x3, nboot = 999)$lim, c(-Inf, 3.75))
})

test_that("hist_out() works with size 1000", {

  set.seed(1)
  x <- rnorm(1000)
  expect_equal(hist_out(x)$lim, c(-Inf, Inf))

  # Easy to separate
  x2 <- c(x, rnorm(50, mean = 7))
  # hist(x2, breaks = nclass.scottRob)
  expect_equal(hist_out(x2)$lim, c(-Inf, 4.25))

  # More difficult to separate
  x3 <- c(x, rnorm(50, mean = 6))
  # hist(x3, breaks = nclass.scottRob)
  expect_equal(hist_out(x3)$lim, c(-Inf, Inf))
  expect_equal(hist_out(x3, nboot = 999)$lim, c(-Inf, 4.75))
})

test_that("hist_out() works with size 100,000", {

  set.seed(1)
  x <- rnorm(100e3)
  expect_equal(hist_out(x)$lim, c(-3.95, 4.05))

  # Easy to separate
  x2 <- c(x, rnorm(200, mean = 6))
  hist(x2, breaks = nclass.scottRob)
  expect_equal(hist_out(x2)$lim, c(-3.95, 4.05))
})
