#include <RcppArmadillo.h>
#include <RcppEigen.h>
using namespace Eigen;

//' Solve (A + lam I) x = b
//'
//' @param A A *symmetric* square matrix.
//' @param b A vector.
//' @param add_to_diag One value to add to the diagonal of A (lam). Default is 0.
//'
//' @return The best solution `x` of this linear system.
//'
//' @examples
//' A <- matrix(rnorm(4), 2); A[1, 2] <- A[2, 1]  # should be symmetric
//' x <- rnorm(2)
//' b <- A %*% x
//' x2 <- drop(solve(A, b))
//' x3 <- solve_linear_system(A, b)
//' rbind(x, x2, x3)
//'
//' @export
// [[Rcpp::export]]
Eigen::VectorXd solve_linear_system(Eigen::MatrixXd& A,
                                    const Eigen::VectorXd& b,
                                    double add_to_diag = 0) {
  int K = A.cols();
  for (int k = 0; k < K; k++) A(k, k) += add_to_diag;
  Eigen::MINRES<Eigen::MatrixXd> minres(A);
  return minres.solve(b);
}
