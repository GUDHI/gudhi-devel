/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2019 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "random_orthogonal_matrix_function"
#include <boost/test/unit_test.hpp>
#include <gudhi/Unitary_tests_utils.h>

#include <gudhi/Functions/random_orthogonal_matrix.h>

#include <string>

#include <random>
#include <cstdlib>

using namespace Gudhi::coxeter_triangulation;

// this test is separated as it requires CGAL
BOOST_AUTO_TEST_CASE(random_orthogonal_matrix_function) {
  // random orthogonal matrix
  Eigen::MatrixXd matrix = random_orthogonal_matrix(5);
  Eigen::MatrixXd id_matrix = matrix.transpose() * matrix;
  for (std::size_t i = 0; i < 5; ++i)
    for (std::size_t j = 0; j < 5; ++j)
      if (i == j)
        GUDHI_TEST_FLOAT_EQUALITY_CHECK(id_matrix(i, j), 1.0, 1e-10);
      else
        GUDHI_TEST_FLOAT_EQUALITY_CHECK(id_matrix(i, j), 0.0, 1e-10);
}
