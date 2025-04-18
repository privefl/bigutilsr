% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/varimax.R
\name{varimax2}
\alias{varimax2}
\title{Varimax rotation}
\usage{
varimax2(X, normalize = FALSE, reorder = TRUE, rotmat = FALSE)
}
\arguments{
\item{X}{A matrix with more rows than columns.}

\item{normalize}{Whether to apply Kaiser normalization? See \link[stats:varimax]{stats::varimax}.
Default is \code{FALSE}.}

\item{reorder}{Whether to permute rotation vectors to maximize the conservation
of the order of the initial columns of \code{X}. Default is \code{TRUE}.}

\item{rotmat}{Whether to return the rotation matrix \code{rot}, or the rotated
matrix \code{X \%*\% rot} (the default, \code{FALSE}).}
}
\value{
Either the rotation matrix \code{rot}, or the rotated matrix \code{X \%*\% rot},
depending on \code{rotmat}.
}
\description{
Varimax rotation
}
\examples{
X <- as.matrix(iris[1:4])
X_rot <- varimax2(X)
X_rot2 <- varimax(X, normalize = FALSE)$loadings[]
all.equal(X_rot2, X_rot[, c(3, 2, 1, 4)], check.attributes = FALSE)
varimax2(X, rotmat = TRUE)

X2 <- prcomp(X)$x
X2_rot <- varimax2(X2)
X2_rot2 <- varimax(X2, normalize = FALSE)$loadings[]
all.equal(X2_rot, X2_rot2, check.attributes = FALSE)
}
