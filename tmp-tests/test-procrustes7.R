library(bigsnpr)
bedfile <- download_1000G("tmp-data")
(obj.bed <- bed(bedfile))

N <- 1800; M <- ncol(obj.bed)
set.seed(1); ind <- sample(nrow(obj.bed), N)
K <- 10
obj.svd <- bed_autoSVD2(obj.bed, ind.row = ind, k = K, ncores = nb_cores())
PC.ref <- predict(obj.svd)

ind2 <- rows_along(obj.bed)[-ind]
keep <- attr(obj.svd, "subset.col")
system.time(
  X2 <- bigsnpr:::read_bed_scaled(obj.bed, ind2, keep,
                                  center = obj.svd$center,
                                  scale = obj.svd$scale)
)

system.time(
  all_proj <- bigutilsr::pca_OADP_proj(X2, obj.svd$v, obj.svd$d)
)
U3 <- all_proj$simple_proj
proj <- all_proj$OADP_proj

plot(proj[, 5:10], U3[, 5:10], col = col(proj[, 5:10]))

col <- 5:6 + 2
plot(PC.ref[, col], xlab = paste0("PC", col[1]), ylab = paste0("PC", col[2]))
points(U3[, col], col = "red", pch = 20)
points(proj[, col], col = "blue", pch = 20)
