X <- readRDS(system.file("testdata", "three-pops.rds", package = "bigutilsr"))
pca <- prcomp(X, scale. = TRUE, rank. = 10)
U <- pca$x

library(bigutilsr)
# dist <- sqrt(covRob(U, estim = "MCD")$dist)
dist <- sqrt(covRob(U, estim = "pairwiseGK")$dist)
library(ggplot2)
theme_set(bigstatsr::theme_bigstatsr(0.8))
qplot(U[, 1], U[, 2], color = dist) + coord_equal() +
  scale_color_viridis_c()


U[1, 1] <- 30
U[1, 1] <- -18.9001

tmp <- performance::check_outliers(as.data.frame(U), method = "optics")
dist2 <- sqrt(attr(tmp, "data")$Distance_OPTICS)

qplot(U[, 1], U[, 2], color = dist2) + coord_equal() +
  scale_color_viridis_c()
qplot(U[, 3], U[, 4], color = dist2) + coord_equal() +
  scale_color_viridis_c()
qplot(U[, 5], U[, 6], color = dist2) + coord_equal() +
  scale_color_viridis_c()
qplot(U[, 7], U[, 8], color = dist2) + coord_equal() +
  scale_color_viridis_c()
qplot(U[, 9], U[, 10], color = dist2) + coord_equal() +
  scale_color_viridis_c()

qplot(dist, dist2)

dist3 <- LOF(U, robMaha = TRUE)
qplot(dist2, dist3)

all_dist <- do.call("cbind", lapply(2:10, function(k) {

  # U.sub <- U[, k + 0:1]

  # tmp <- performance::check_outliers(as.data.frame(U.sub), method = "optics")
  # dist2 <- sqrt(attr(tmp, "data")$Distance_OPTICS)

  # LOF(U.sub)

  # cbind(LOF(U[, 1:k], log = FALSE),
  #       sqrt(covRob(U[, k:10], estim = "pairwiseGK")$dist))

  # sqrt(LOF(U[, 1:k], log = FALSE))
  LOF(U[, 1:k])
}))
round(100 * cor(all_dist), 1)
round(100 * cor(all_dist, method = "kendall"), 1)
plot(all_dist[, c(6, 9)])
dist4 <- sqrt(covRob(all_dist, estim = "pairwiseGK")$dist)
hist(dist4, nclass.scottRob)
dist4 <- apply(all_dist, 1, max)
hist(dist4, nclass.scottRob)
abline(v = tukey_mc_up(dist4), col = "red")
abline(v = hist_out(dist4)$lim[2], col = "blue")

all_dist2 <- matrixStats::rowCummaxs(all_dist)
round(100 * cor(all_dist2), 1)
round(100 * cor(all_dist2, method = "kendall"), 1)
plot(all_dist2[, c(1, 9)])

qplot(U[, 1], U[, 2], color = dist4) + coord_equal() +
  scale_color_viridis_c()
qplot(U[, 1], U[, 2], color = dist4 > tukey_mc_up(dist4)) + coord_equal() +
  labs(color = NULL)
qplot(U[, 1], U[, 2], color = dist4 > hist_out(dist4)$lim[2]) + coord_equal() +
  labs(color = NULL)
qplot(U[, 3], U[, 4], color = dist4) + coord_equal() +
  scale_color_viridis_c()
qplot(U[, 5], U[, 6], color = dist4) + coord_equal() +
  scale_color_viridis_c()
qplot(U[, 7], U[, 8], color = dist4) + coord_equal() +
  scale_color_viridis_c()
qplot(U[, 9], U[, 10], color = dist4) + coord_equal() +
  scale_color_viridis_c()

