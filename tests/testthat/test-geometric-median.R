################################################################################

context("GEO_MEDIAN")

################################################################################

test_that("geometric_median() works", {

  X <- readRDS(system.file("testdata", "three-pops.rds", package = "bigutilsr"))
  svd <- svds(scale(X), k = 5)

  U <- sweep(svd$u, 2, svd$d, '*')
  med_all <- geometric_median(U)
  expect_equal(med_all, drop(Gmedian::Weiszfeld(U)$median))

  pop <- rep(1:3, c(143, 167, 207))
  med_pop <- geometric_median(U, by_grp = pop)
  expect_equal(dim(med_pop), c(3, 5))

  med_pop2 <- do.call("rbind", lapply(split(seq_along(pop), pop), function(ind) {
    Gmedian::Weiszfeld(U[ind, ])$median
  }))
  expect_equal(med_pop, med_pop2, check.attributes = FALSE)
})

################################################################################
