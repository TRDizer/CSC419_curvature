#include "../include/internal_angles.h"
#include <cmath>

#define GET_ADJ1(index)    ((index + 1) % 3)
#define GET_ADJ2(index)    ((index + 2) % 3)

void internal_angles(
  const Eigen::MatrixXd & l_sqr,
  Eigen::MatrixXd & A)
{
  // Add with your code
  A.resizeLike(l_sqr);
  A.setZero();

  double ops, adj1, adj2;
  for (int i = 0; i < l_sqr.rows(); i++) {
    for (int j = 0; j < 3; j++) {
      ops = l_sqr(i, j);
      adj1 = (i, GET_ADJ1(j));
      adj2 = (i, GET_ADJ2(j));

      A(i, j) = std::acos((adj1 + adj2 - ops) / (2.0 * std::sqrt(adj1 * adj2)));
    }
  }
}
