/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2020 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <CGAL/Epeck_d.h>

#include <vector>
#include <random>
#include <cmath> // for std::fabs

#include <gudhi/Alpha_complex.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Unitary_tests_utils.h>

template<class Kernel>
void do_test() {
  // Check that in exact mode for static dimension 4 the code for dD unweighted and for dD weighted with all weights
  // 0 give exactly the same simplex tree (simplices and filtration values).

  // Random points construction
  using Point_d = typename Kernel::Point_d;
  std::vector<Point_d> points;
  std::uniform_real_distribution<double> rd_pts(-10., 10.);
  std::random_device rand_dev;
  std::mt19937 rand_engine(rand_dev());
  for (int idx = 0; idx < 20; idx++) {
    std::vector<double> point {rd_pts(rand_engine), rd_pts(rand_engine), rd_pts(rand_engine), rd_pts(rand_engine)};
    points.emplace_back(point.begin(), point.end());
  }
  
  // Alpha complex from points
  Gudhi::alpha_complex::Alpha_complex<Kernel, false> alpha_complex_from_points(points);
  Gudhi::Simplex_tree<> simplex;
  Gudhi::Simplex_tree<>::Filtration_value infty = std::numeric_limits<Gudhi::Simplex_tree<>::Filtration_value>::infinity();
  BOOST_CHECK(alpha_complex_from_points.create_complex(simplex, infty, true));
  std::clog << "Iterator on alpha complex simplices in the filtration order, with [filtration value]:"
            << std::endl;
  for (auto f_simplex : simplex.filtration_simplex_range()) {
    std::clog << "   ( ";
    for (auto vertex : simplex.simplex_vertex_range(f_simplex)) {
      std::clog << vertex << " ";
    }
    std::clog << ") -> " << "[" << simplex.filtration(f_simplex) << "] " << std::endl;
  }

  // Alpha complex from zero weighted points
  std::vector<typename Kernel::FT> weights(20, 0.);
  Gudhi::alpha_complex::Alpha_complex<Kernel, true> alpha_complex_from_zero_weighted_points(points, weights);
  Gudhi::Simplex_tree<> zw_simplex;
  BOOST_CHECK(alpha_complex_from_zero_weighted_points.create_complex(zw_simplex, infty, true));

  std::clog << "Iterator on zero weighted alpha complex simplices in the filtration order, with [filtration value]:"
            << std::endl;
  for (auto f_simplex : zw_simplex.filtration_simplex_range()) {
    std::clog << "   ( ";
    for (auto vertex : zw_simplex.simplex_vertex_range(f_simplex)) {
      std::clog << vertex << " ";
    }
    std::clog << ") -> " << "[" << zw_simplex.filtration(f_simplex) << "] " << std::endl;
  }

  BOOST_CHECK(zw_simplex == simplex);
}
