################################################################################

context("UTILS")

################################################################################

test_that("as_model_matrix() works", {

  my_iris <- datasets::iris
  mat0 <- as_model_matrix(my_iris)
  expect_equal(dim(mat0), c(150, 4 + 3))

  # expect_identical(as_model_matrix(tibble    ::as_tibble    (my_iris)), mat0)
  # expect_identical(as_model_matrix(data.table::as.data.table(my_iris)), mat0)

  expect_error(as_model_matrix(as.matrix(my_iris)), "not implemented")
  expect_error(as_model_matrix(   unlist(my_iris)), "not implemented")

  expect_equal(as_model_matrix(setNames(my_iris, NULL)), mat0,
               check.attributes = FALSE)

  my_iris$Species_chr <- as.character(my_iris$Species)
  mat2 <- as_model_matrix(my_iris)
  expect_equal(dim(mat2), c(150, 4 + 3 + 2))

  mat3 <- as_model_matrix(my_iris, intercept = TRUE)
  expect_equal(dim(mat2), c(150, 1 + 4 + 2 + 2))
})

################################################################################
