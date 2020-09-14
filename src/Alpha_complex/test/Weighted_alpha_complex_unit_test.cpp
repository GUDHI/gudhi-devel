/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2020 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "weighted_alpha_complex"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <CGAL/Epick_d.h>
#include <CGAL/Epeck_d.h>

#include <cmath>  // float comparison
#include <vector>
#include <random>

#include <gudhi/Alpha_complex.h>
#include <gudhi/Alpha_complex_3d.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Unitary_tests_utils.h>

// Use dynamic_dimension_tag for the user to be able to set dimension
typedef CGAL::Epeck_d< CGAL::Dynamic_dimension_tag > Exact_kernel_d;
// Use static dimension_tag for the user not to be able to set dimension
typedef CGAL::Epeck_d< CGAL::Dimension_tag<4> > Exact_kernel_s;
// Use dynamic_dimension_tag for the user to be able to set dimension
typedef CGAL::Epick_d< CGAL::Dynamic_dimension_tag > Inexact_kernel_d;
// Use static dimension_tag for the user not to be able to set dimension
typedef CGAL::Epick_d< CGAL::Dimension_tag<4> > Inexact_kernel_s;

typedef boost::mpl::list<Exact_kernel_d, Exact_kernel_s, Inexact_kernel_d, Inexact_kernel_s> list_of_kernel_variants;

BOOST_AUTO_TEST_CASE_TEMPLATE(Zero_weighted_alpha_complex, Kernel, list_of_kernel_variants) {
  // Random points construction
  using Point_d = typename Kernel::Point_d;
  std::vector<Point_d> points;
  std::uniform_real_distribution<double> rd_pts(-10., 10.);
  std::random_device rand_dev;
  std::mt19937 rand_engine(rand_dev());
  for (int idx = 0; idx < 20; idx++) {
    std::vector<double> point {rd_pts(rand_engine), rd_pts(rand_engine), rd_pts(rand_engine), rd_pts(rand_engine)};
    points.emplace_back(Point_d(point.begin(), point.end()));
  }
  
  // Alpha complex from points
  Gudhi::alpha_complex::Alpha_complex<Kernel, false> alpha_complex_from_points(points);
  Gudhi::Simplex_tree<> simplex;
  BOOST_CHECK(alpha_complex_from_points.create_complex(simplex));
  std::clog << "Iterator on alpha complex simplices in the filtration order, with [filtration value]:" << std::endl;
  for (auto f_simplex : simplex.filtration_simplex_range()) {
    std::clog << "   ( ";
    for (auto vertex : simplex.simplex_vertex_range(f_simplex)) {
      std::clog << vertex << " ";
    }
    std::clog << ") -> " << "[" << simplex.filtration(f_simplex) << "] ";
    std::clog << std::endl;
  }

  // Alpha complex from zero weighted points
  std::vector<typename Kernel::FT> weights(20, 0.);
  Gudhi::alpha_complex::Alpha_complex<Kernel, true> alpha_complex_from_zero_weighted_points(points, weights);
  Gudhi::Simplex_tree<> zw_simplex;
  BOOST_CHECK(alpha_complex_from_zero_weighted_points.create_complex(zw_simplex));

  std::clog << "Iterator on zero weighted alpha complex simplices in the filtration order, with [filtration value]:" << std::endl;
  for (auto f_simplex : zw_simplex.filtration_simplex_range()) {
    std::clog << "   ( ";
    for (auto vertex : zw_simplex.simplex_vertex_range(f_simplex)) {
      std::clog << vertex << " ";
    }
    std::clog << ") -> " << "[" << zw_simplex.filtration(f_simplex) << "] ";
    std::clog << std::endl;
  }

  BOOST_CHECK(zw_simplex == simplex);
}