################################################################################

context("VARIMAX")

################################################################################

test_that("varimax2() works", {

  replicate(10, {

    for (scale in c(TRUE, FALSE)) {
      for (ord in c(TRUE, FALSE)) {
        X <- matrix(rnorm(1000 * 5), ncol = 5)
        X_rot <- varimax2(X, normalize = scale, reorder = ord)
        rot   <- varimax2(X, normalize = scale, reorder = ord, rotmat = TRUE)
        expect_equal(X_rot, X %*% rot)
        if (ord) {
          X_rot2 <- varimax(X, normalize = scale)$loadings[]
          expect_equal(sum(cor(X_rot2, X_rot) > 0.9999), 5)
        } else {
          expect_equal(X_rot, varimax(X, normalize = scale)$loadings[])
          expect_equal(rot,   varimax(X, normalize = scale)$rotmat)
        }
      }
    }
  })

})

################################################################################
