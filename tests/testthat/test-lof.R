context("test-lof")

test_that("to_maha() works", {

  X <- readRDS(system.file("testdata", "three-pops.rds", package = "bigutilsr"))
  svd <- svds(scale(X), k = 5)
  U <- svd$u
  maha <- covRob(U, estim = "pairwiseGK")
  U.maha <- to_maha(U)
  expect_equal(rowSums(U.maha^2), maha$dist)

  mat <- matrix(rnorm(500), 100, 5)
  mat.maha <- attr(U.maha, "trans")(mat)
  expect_equal(rowSums(mat.maha^2),
               stats::mahalanobis(mat, maha$center, maha$cov))
})

test_that("LOF() works", {

  X <- readRDS(system.file("testdata", "three-pops.rds", package = "bigutilsr"))
  svd <- svds(scale(X), k = 10)
  # plot(head(svd$d, -1), log = "xy")
  llof <- LOF(svd$u)
  expect_identical(LOF(svd$u, ncores = 2), llof)
  # hist(llof, breaks = 20)
  lof <- LOF(svd$u, log = FALSE)
  expect_equal(log(lof), llof)
  expect_equal(sum(llof > tukey_mc_up(llof)), 0)
  llof_maha <- LOF(svd$u, robMaha = TRUE)
  # hist(llof_maha)
  expect_equal(sum(llof_maha > tukey_mc_up(llof_maha)), 0)

  X2 <- rbind(X[1, ], X)  # duplicate first row
  svd2 <- svds(scale(X2), k = 10)
  llof2 <- LOF(svd2$u)
  expect_equal(which(llof2 > tukey_mc_up(llof2)), 1:2)
  llof2_maha <- LOF(svd2$u, robMaha = TRUE)
  expect_equal(which(llof2_maha > tukey_mc_up(llof2_maha)), 1:2)
})
