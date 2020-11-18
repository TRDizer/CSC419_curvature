#include "../include/mean_curvature.h"
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <igl/per_vertex_normals.h>

void mean_curvature(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & H)
{
  // Replace with your code
  Eigen::SparseMatrix<double> lapacian;
  igl::cotmatrix(V, F, lapacian);

  Eigen::SparseMatrix<double> mass;
  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, mass);


  Eigen::SparseMatrix<double> mass_inverse;
  igl::invert_diag(mass, mass_inverse);
  Eigen::MatrixXd Hn = mass_inverse * lapacian * V;

  Eigen::MatrixXd normals;
  igl::per_vertex_normals(V, F, normals);

  H.resize(V.rows());
  H.setZero();

  auto get_sign = [](double &signness_dot) {
    if (signness_dot > 0) {
      return 1.0;
    }
    else if (signness_dot < 0) {
      return -1.0;
    }
    else {
      // not sure how often this shows up, but the most correct way of handling should be an epislon check instead of 
      // straight equality
      return 0.0;
    }
  };

  for (int i = 0; i < V.rows(); i++) {
    double signness_dot = Hn.row(i).dot(normals.row(i));
    H(i) = get_sign(signness_dot) * Hn.row(i).norm(); // adjust curving direction accordingly
  }
}
