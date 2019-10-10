context("test-lof")

test_that("LOF() works", {

  X <- readRDS(system.file("testdata", "three-pops.rds", package = "bigutilsr"))
  svd <- svd(scale(X))
  # plot(head(svd$d, -1), log = "xy")
  llof <- LOF(svd$u[, 1:10])
  expect_identical(LOF(svd$u[, 1:10], ncores = 2), llof)
  # hist(llof, breaks = 20)
  lof <- LOF(svd$u[, 1:10], log = FALSE)
  expect_equal(log(lof), llof)
  expect_equal(sum(llof > tukey_mc_up(llof)), 0)
  llof_maha <- LOF(svd$u[, 1:10], robMaha = TRUE)
  # hist(llof_maha)
  expect_equal(sum(llof_maha > tukey_mc_up(llof_maha)), 0)

  X2 <- rbind(X[1, ], X)  # duplicate first row
  svd2 <- svd(scale(X2))
  llof2 <- LOF(svd2$u[, 1:10])
  expect_equal(which(llof2 > tukey_mc_up(llof2)), 1:2)
  llof2_maha <- LOF(svd2$u[, 1:10], robMaha = TRUE)
  expect_equal(which(llof2_maha > tukey_mc_up(llof2_maha)), 1:2)
})
