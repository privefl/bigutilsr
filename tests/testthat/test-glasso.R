################################################################################

context("GEO_MEDIAN")

################################################################################

test_that("regul_glasso() works", {

  cov <- cov(iris[1:4])
  lambda <- 1 / sqrt(nrow(iris))
  cov_regul <- regul_glasso(cov, lambda)
  expect_equal(diag(cov_regul), diag(cov))
  expect_lte(max(abs(cov_regul - cov)), lambda * 1.0001)

  cov_regul2 <- regul_glasso(cov * 10, lambda * 10, tol = 1e-3)
  expect_equal(cov_regul2, structure(cov_regul * 10, lambda = lambda * 10),
               tolerance = 1e-4)

  cov_regul3 <- regul_glasso(cov / 10, lambda / 10, tol = 1e-5)
  expect_equal(cov_regul3, structure(cov_regul / 10, lambda = lambda / 10),
               tolerance = 1e-4)

  expect_error(regul_glasso(as.data.frame(cov), lambda),
               "Not compatible with requested type")
})

################################################################################
