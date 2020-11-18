#include "../include/angle_defect.h"
#include "../include/internal_angles.h"

#define _USE_MATH_DEFINES

#include <math.h>
#include <igl/squared_edge_lengths.h>
#include <iostream>

void angle_defect(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & D)
{
  D = Eigen::VectorXd::Zero(V.rows());

  Eigen::MatrixXd sqr_edges;
  igl::squared_edge_lengths(V, F, sqr_edges);

  Eigen::MatrixXd angles;
  internal_angles(sqr_edges, angles);

  D.setConstant(2.0 * M_PI); // initial no defect value

  for (int i = 0; i < F.rows(); i++) {
    for (int j = 0; j < 3; j++) {
      D(F(i, j)) -= angles(i, j);
    }
  }
}
