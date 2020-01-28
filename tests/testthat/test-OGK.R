################################################################################

context("OGK")

################################################################################

expect_error(covrob_ogk(iris), "'U' is not of class 'matrix'.")
mat <- as_model_matrix(iris)
expect_error(covrob_ogk(mat), "Problem with columns [5, 6, 7] of 'U'.", fixed = TRUE)
mat[1, 1] <- NA
expect_error(covrob_ogk(mat), "You can't have missing values in 'U'.")
mat2 <- as.matrix(iris[1:4])
expect_type(covrob_ogk(mat2), "list")

test_that("covrob_ogk() is the same as rrcov::CovOgk()", {

  skip_if_not_installed("mvtnorm")
  skip_if_not_installed("rrcov")

  replicate(100, {
    n <- sample(200:2000, 1)
    p <- sample(3:10, 1)
    sigma1 <- crossprod(matrix(rnorm(n * p),                       nrow = n))
    sigma2 <- crossprod(matrix(rnorm(n * p, sd = runif(1, 1, 10)), nrow = n))
    X <- rbind(
      mvtnorm::rmvnorm(n,                 mean = rep(0, p), sigma = sigma1),
      mvtnorm::rmvnorm(sample(10:200, 1), mean = rnorm(p),  sigma = sigma2)
    )
    beta <- runif(1, 0.7, 0.95)
    niter <- sample(0:5, 1)
    test <- covrob_ogk(X, niter = niter, beta = beta)
    true <- rrcov::CovOgk(X, niter = niter, beta = beta)
    expect_equal(test$cov, true@cov)
    expect_equal(test$center, true@center)
  })
})

################################################################################
