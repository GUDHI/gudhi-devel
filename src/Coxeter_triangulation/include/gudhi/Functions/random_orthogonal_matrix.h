/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2019 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef FUNCTIONS_RANDOM_ORTHOGONAL_MATRIX_H_
#define FUNCTIONS_RANDOM_ORTHOGONAL_MATRIX_H_

#include <cstdlib>  // for std::size_t
#include <cmath>    // for std::cos, std::sin
#include <random>   // for std::uniform_real_distribution, std::random_device

#include <gudhi/math.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SVD>

#include <CGAL/Epick_d.h>
#include <CGAL/point_generators_d.h>

namespace Gudhi {

namespace coxeter_triangulation {

/** \brief Generates a uniform random orthogonal matrix using the "subgroup algorithm" by
 * Diaconis & Shashahani.
 * \details Taken from https://en.wikipedia.org/wiki/Rotation_matrix#Uniform_random_rotation_matrices.
 * The idea: take a random rotation matrix of dimension d-1, embed it
 * as a d*d matrix M with the last column (0,...,0,1).
 * Pick a random vector v on a sphere S^d. rotate the matrix M so that its last column is v.
 * The determinant of the matrix can be either 1 or -1
 */
// Note: the householderQR operation at the end seems to take a lot of time at compilation.
// The CGAL headers are another source of long compilation time.
Eigen::MatrixXd random_orthogonal_matrix(std::size_t d) {
  typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag> Kernel;
  typedef typename Kernel::Point_d Point_d;
  if (d == 1) return Eigen::VectorXd::Constant(1, 1.0);
  if (d == 2) {
    // 0. < alpha < 2 Pi
    std::uniform_real_distribution<double> unif(0., 2 * Gudhi::PI);
    std::random_device rand_dev;
    std::mt19937 rand_engine(rand_dev());
    double alpha = unif(rand_engine);

    Eigen::Matrix2d rot;
    rot << std::cos(alpha), -std::sin(alpha), std::sin(alpha), cos(alpha);
    return rot;
  }
  Eigen::MatrixXd low_dim_rot = random_orthogonal_matrix(d - 1);
  Eigen::MatrixXd rot(d, d);
  Point_d v = *CGAL::Random_points_on_sphere_d<Point_d>(d, 1);
  for (std::size_t i = 0; i < d; ++i) rot(i, 0) = v[i];
  for (std::size_t i = 0; i < d - 1; ++i)
    for (std::size_t j = 1; j < d - 1; ++j) rot(i, j) = low_dim_rot(i, j - 1);
  for (std::size_t j = 1; j < d; ++j) rot(d - 1, j) = 0;
  rot = rot.householderQr()
            .householderQ();  // a way to do Gram-Schmidt, see https://forum.kde.org/viewtopic.php?f=74&t=118568#p297246
  return rot;
}

}  // namespace coxeter_triangulation

}  // namespace Gudhi

#endif
