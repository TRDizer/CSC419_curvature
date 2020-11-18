#include "../include/principal_curvatures.h"
#include <igl/adjacency_matrix.h>
#include <set>
#include <igl/pinv.h>
#include <cmath>

#define sqr(x)    std::pow(x, 2)

void principal_curvatures(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & D1,
  Eigen::MatrixXd & D2,
  Eigen::VectorXd & K1,
  Eigen::VectorXd & K2)
{
  // Replace with your code
  K1 = Eigen::VectorXd::Zero(V.rows());
  K2 = Eigen::VectorXd::Zero(V.rows());
  D1 = Eigen::MatrixXd::Zero(V.rows(),3);
  D2 = Eigen::MatrixXd::Zero(V.rows(),3);

  // Gather sampling points in the vicinity by combining primary adjacency (directly adjacent) and secondary adjacency (two-ring)
  Eigen::SparseMatrix<int> adjacency_matrix;
  igl::adjacency_matrix(F, adjacency_matrix);

  // Ref (god-sent): https://math.stackexchange.com/questions/207826/algorithm-to-find-all-vertices-exactly-k-steps-away-in-an-undirected-graph
  // samples = direct adjacency + secondary adjacency
  Eigen::SparseMatrix<int> samples_mat = adjacency_matrix + adjacency_matrix * adjacency_matrix;

  for (int i = 0; i < V.rows(); i++) {
    // Assemble P
    std::set<int> sample_verteices;
    for (Eigen::SparseMatrix<int>::InnerIterator iter(samples_mat, i); iter; ++iter) {
      sample_verteices.insert(iter.row());
    }

    Eigen::MatrixXd P(sample_verteices.size(), 3);
    int P_row_idx = 0;
    for (int index : sample_verteices) {
      P.row(P_row_idx) = V.row(index) - V.row(i);
      P_row_idx++;
    }

    // PCA of P^T * P and assembly S and B
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> pca(P.transpose() * P);
    Eigen::Vector3d pc1 = pca.eigenvectors().col(2);
    Eigen::Vector3d pc2 = pca.eigenvectors().col(1);
    Eigen::Vector3d height = pca.eigenvectors().col(0);

    Eigen::MatrixXd S(P.rows(), 2);
    S.col(0) = P * pc1;
    S.col(1) = P * pc2;

    Eigen::VectorXd B = P * height;

    Eigen::MatrixXd least_sqr_goal(P.rows(), 5);
    // w = a1 * u + a2 * v + a3 * u^2 + a4 * u * w + a5 * v^2
    least_sqr_goal << S.col(0), S.col(1), (S.col(0).array().square()).matrix(), S.col(0).cwiseProduct(S.col(1)), (S.col(1).array().square()).matrix();
  
    Eigen::MatrixXd least_sqr_inv;
    igl::pinv(least_sqr_goal, least_sqr_inv);

    Eigen::VectorXd a = least_sqr_inv * B;

    //       | e   f | | E    F |^-1
    // S = - |       | |        |
    //       | f   g | | F    G |

    double E = 1 + sqr(a(0));
    double F_term = a(0) * a(1);
    double G = 1 + sqr(a(1));
    double denom = std::sqrt(sqr(a(0)) + 1 + sqr(a(1)));
    double e = 2.0 * a(2) / denom;
    double f = a(3) / denom;
    double g = 2.0 * a(4) / denom;

    Eigen::Matrix2d temp_mat1, temp_mat2, shape_operator;
    temp_mat1 << e, f,
                 f, g;
    temp_mat2 << E,      F_term,
                 F_term, G;
    shape_operator = -temp_mat1 * temp_mat2.inverse();

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> decomp(shape_operator);
    K1(i) = decomp.eigenvalues()(0);
    K2(i) = decomp.eigenvalues()(1);
    D1.row(i) = decomp.eigenvectors()(1,1) * pc1 + decomp.eigenvectors()(1,0) * pc2;
    D2.row(i) = decomp.eigenvectors()(0,1) * pc1 + decomp.eigenvectors()(0,0) * pc2;
  }
}
