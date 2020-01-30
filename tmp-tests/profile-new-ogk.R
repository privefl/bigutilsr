mat <- do.call("cbind", iris)
mat2 <- mat[rep(1:150, 5000), ]
object.size(mat2) / 1024^2  # 29 MB

profvis::profvis({test2 <- bigutilsr::covrob_ogk(mat2)})  # 2 sec / 1.3 GB

# profvis::profvis({test <- covRob_rcpp(mat2)})  # 2.7 sec / 1.6 GB
# profvis::profvis({true <- robust::covRob(mat2, estim = "pairwiseGK")})  # 33 sec / 0.4 GB
# profvis::profvis({test3 <- robustbase::covOGK(mat2, sigmamu = robustbase::scaleTau2)}) # 4 sec and 4GB
# profvis::profvis({test4 <- rrcov::CovOgk(mat2)}) # 2 sec and 300 MB
# all.equal(true[c("cov", "center", "dist")], test[c("cov", "center", "dist")],
#           check.attributes = FALSE)
# all.equal(true[c("cov", "center", "dist")], test2[c("cov", "center", "dist")],
#           check.attributes = FALSE)
# all.equal(true[c("cov", "center", "dist")], test3[c("cov", "center", "distances")],
#           check.attributes = FALSE)
