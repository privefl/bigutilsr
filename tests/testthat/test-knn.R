################################################################################

context("test-knn")

options(bigstatsr.check.parallel.blas = FALSE)

################################################################################

test_that("knn_parallel() works", {

  N <- 100
  M <- 50
  K <- 10
  mat1 <- matrix(rnorm(N * K), N, K)
  mat2 <- matrix(rnorm(M * K), M, K)

  # sequential
  expect_identical(nabor::knn(mat1, k = 8),
                   knn_parallel(mat1, k = 8, ncores = 1))
  expect_identical(nabor::knn(mat1, mat2, k = 8),
                   knn_parallel(mat1, mat2, k = 8, ncores = 1))
  # parallel
  expect_identical(nabor::knn(mat1, k = 8),
                   knn_parallel(mat1, k = 8, ncores = 2))
  expect_identical(nabor::knn(mat1, mat2, k = 8),
                   knn_parallel(mat1, mat2, k = 8, ncores = 2))
  # with extra params
  expect_failure(expect_identical(nabor::knn(mat1, k = 8, eps = 1),
                                  knn_parallel(mat1, k = 8, ncores = 2)))
  expect_identical(nabor::knn(mat1, k = 8, eps = 1),
                   knn_parallel(mat1, k = 8, eps = 1, ncores = 2))
})

################################################################################
