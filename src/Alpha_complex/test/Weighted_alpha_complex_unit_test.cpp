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

template <typename Point_d>
bool cgal_3d_point_sort (Point_d a,Point_d b) {
  if (a[0] != b[0])
    return a[0] < b[0];
  if (a[1] != b[1])
    return a[1] < b[1];
  return a[2] < b[2];
}

/*BOOST_AUTO_TEST_CASE(Weighted_alpha_complex_3d_comparison) {
  // Random points construction
  using Kernel_dD = CGAL::Epeck_d< CGAL::Dimension_tag<3> >;
  using Bare_point_d = typename Kernel_dD::Point_d;
  using Weighted_point_d = typename Kernel_dD::Weighted_point_d;
  std::vector<Weighted_point_d> w_points_d;

  using Exact_weighted_alpha_complex_3d =
    Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::EXACT, true, false>;
  using Bare_point_3 = typename Exact_weighted_alpha_complex_3d::Bare_point_3;
  using Weighted_point_3 = typename Exact_weighted_alpha_complex_3d::Weighted_point_3;
  std::vector<Weighted_point_3> w_points_3;

  std::uniform_real_distribution<double> rd_pts(-10., 10.);
  std::uniform_real_distribution<double> rd_wghts(-0.5, 0.5);
  std::random_device rand_dev;
  std::mt19937 rand_engine(rand_dev());
  for (int idx = 0; idx < 20; idx++) {
    std::vector<double> point {rd_pts(rand_engine), rd_pts(rand_engine), rd_pts(rand_engine)};
    double weight = rd_wghts(rand_engine);
    w_points_d.emplace_back(Weighted_point_d(Bare_point_d(point.begin(), point.end()), weight));
    w_points_3.emplace_back(Weighted_point_3(Bare_point_3(point[0], point[1], point[2]), weight));
  }

  // Weighted alpha complex for dD version
  Gudhi::alpha_complex::Alpha_complex<Kernel_dD, true> alpha_complex_dD_from_weighted_points(w_points_d);
  Gudhi::Simplex_tree<> w_simplex_d;
  BOOST_CHECK(alpha_complex_dD_from_weighted_points.create_complex(w_simplex_d));

  std::clog << "Iterator on weighted alpha complex dD simplices in the filtration order, with [filtration value]:" << std::endl;
  for (auto f_simplex : w_simplex_d.filtration_simplex_range()) {
    std::clog << "   ( ";
    std::vector<Bare_point_d> points;
    for (auto vertex : w_simplex_d.simplex_vertex_range(f_simplex)) {
      points.emplace_back(alpha_complex_dD_from_weighted_points.get_point(vertex).point());
    }
    std::sort (points.begin(), points.end(), cgal_3d_point_sort<Bare_point_d>);
    for (auto point : points) {
      std::clog << point[0] << " " << point[1] << " " << point[2] << " | ";
    }
    std::clog << ") -> " << "[" << w_simplex_d.filtration(f_simplex) << "] ";
    std::clog << std::endl;
  }

  // Weighted alpha complex for 3D version
  Exact_weighted_alpha_complex_3d alpha_complex_3D_from_weighted_points(w_points_3);
  Gudhi::Simplex_tree<> w_simplex_3;
  BOOST_CHECK(alpha_complex_3D_from_weighted_points.create_complex(w_simplex_3));

  std::clog << "Iterator on weighted alpha complex 3D simplices in the filtration order, with [filtration value]:" << std::endl;
  for (auto f_simplex : w_simplex_3.filtration_simplex_range()) {
    std::clog << "   ( ";
    std::vector<Bare_point_3> points;
    for (auto vertex : w_simplex_3.simplex_vertex_range(f_simplex)) {
      points.emplace_back(alpha_complex_3D_from_weighted_points.get_point(vertex).point());
    }
    std::sort (points.begin(), points.end());
    for (auto point : points) {
      std::clog << point[0] << " " << point[1] << " " << point[2] << " | ";
    }
    std::clog << ") -> " << "[" << w_simplex_3.filtration(f_simplex) << "] ";
    std::clog << std::endl;
  }

  BOOST_CHECK(w_simplex_d == w_simplex_3);

}*/