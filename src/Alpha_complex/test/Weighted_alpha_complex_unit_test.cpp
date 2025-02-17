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

#include <CGAL/Epeck_d.h>

#include <vector>
#include <random>
#include <array>
#include <cmath> // for std::fabs
#include <stdexcept>

#include <gudhi/Alpha_complex.h>
#include <gudhi/Alpha_complex_3d.h>
#include <gudhi/Simplex_tree.h>

BOOST_AUTO_TEST_CASE(Weighted_alpha_complex_3d_comparison) {
  // check that for random weighted 3d points in safe mode the 3D and dD codes give the same result with some tolerance

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
    w_points_d.emplace_back(Bare_point_d(point.begin(), point.end()), weight);
    w_points_3.emplace_back(Bare_point_3(point[0], point[1], point[2]), weight);
  }

  // Structures necessary for comparison
  using Points = std::vector<std::array<double,3>>;
  using Points_and_filtrations = std::map<Points, double>;
  Points_and_filtrations pts_fltr_dD;
  Points_and_filtrations pts_fltr_3d;

  // Weighted alpha complex for dD version
  Gudhi::alpha_complex::Alpha_complex<Kernel_dD, true> alpha_complex_dD_from_weighted_points(w_points_d);
  Gudhi::Simplex_tree<> w_simplex_d;
  BOOST_CHECK(alpha_complex_dD_from_weighted_points.create_complex(w_simplex_d));

  std::clog << "Iterator on weighted alpha complex dD simplices in the filtration order, with [filtration value]:"
            << std::endl;
  for (auto f_simplex : w_simplex_d.filtration_simplex_range()) {
    Points points;
    for (auto vertex : w_simplex_d.simplex_vertex_range(f_simplex)) {
      CGAL::NT_converter<Kernel_dD::RT, double> cgal_converter;
      Bare_point_d pt = alpha_complex_dD_from_weighted_points.get_point(vertex).point();
      points.push_back({cgal_converter(pt[0]), cgal_converter(pt[1]), cgal_converter(pt[2])});
    }
    std::clog << "   ( ";
    std::sort (points.begin(), points.end());
    for (auto point : points) {
      std::clog << point[0] << " " << point[1] << " " << point[2] << " | ";
    }
    std::clog << ") -> " << "[" << w_simplex_d.filtration(f_simplex) << "] ";
    std::clog << std::endl;
    pts_fltr_dD[points] = w_simplex_d.filtration(f_simplex);
  }

  // Weighted alpha complex for 3D version
  Exact_weighted_alpha_complex_3d alpha_complex_3D_from_weighted_points(w_points_3);
  Gudhi::Simplex_tree<> w_simplex_3;
  BOOST_CHECK(alpha_complex_3D_from_weighted_points.create_complex(w_simplex_3));

  std::clog << "Iterator on weighted alpha complex 3D simplices in the filtration order, with [filtration value]:"
            << std::endl;
  for (auto f_simplex : w_simplex_3.filtration_simplex_range()) {
    Points points;
    for (auto vertex : w_simplex_3.simplex_vertex_range(f_simplex)) {
      Bare_point_3 pt = alpha_complex_3D_from_weighted_points.get_point(vertex).point();
      CGAL::NT_converter<Exact_weighted_alpha_complex_3d::Kernel::RT, double> cgal_converter;
      points.push_back({cgal_converter(pt[0]), cgal_converter(pt[1]), cgal_converter(pt[2])});
    }
    std::clog << "   ( ";
    std::sort (points.begin(), points.end());
    for (auto point : points) {
      std::clog << point[0] << " " << point[1] << " " << point[2] << " | ";
    }
    std::clog << ") -> " << "[" << w_simplex_3.filtration(f_simplex) << "] " << std::endl;
    pts_fltr_3d[points] = w_simplex_d.filtration(f_simplex);
  }

  // Compares structures
  auto d3_itr = pts_fltr_3d.begin();
  auto dD_itr = pts_fltr_dD.begin();
  for (; d3_itr != pts_fltr_3d.end() && dD_itr != pts_fltr_dD.end(); ++d3_itr) {
    if (d3_itr->first != dD_itr->first) {
      for(auto point : d3_itr->first)
        std::clog << point[0] << " " << point[1] << " " << point[2] << " | ";
      std::clog << " versus ";
      for(auto point : dD_itr->first)
        std::clog << point[0] << " " << point[1] << " " << point[2] << " | ";
      std::clog << std::endl;
      BOOST_CHECK(false);
    }
    // In safe mode, relative error is less than 1e-5 (can be changed with set_relative_precision_of_to_double)
    if (std::fabs(d3_itr->second - dD_itr->second) > 1e-5 * (std::fabs(d3_itr->second) + std::fabs(dD_itr->second))) {
      std::clog << d3_itr->second << " versus " << dD_itr->second << " diff "
                << std::fabs(d3_itr->second - dD_itr->second) << std::endl;
      BOOST_CHECK(false);
    }
    ++dD_itr;
  }
}

BOOST_AUTO_TEST_CASE(Is_weighted_alpha_complex_nan) {
  using Kernel = CGAL::Epeck_d< CGAL::Dimension_tag<3> >;
  using Bare_point = Kernel::Point_d;
  using Weighted_point = Kernel::Weighted_point_d;
  using Vector_of_points = std::vector<Weighted_point>;

  Vector_of_points points;
  points.emplace_back(Bare_point(1, -1, -1), 4.);
  points.emplace_back(Bare_point(-1, 1, -1), 4.);
  points.emplace_back(Bare_point(-1, -1, 1), 4.);
  points.emplace_back(Bare_point(1, 1, 1),   4.);
  points.emplace_back(Bare_point(2, 2, 2),   1.);

  Gudhi::alpha_complex::Alpha_complex<Kernel, true> alpha_complex_from_weighted_points(points);

  std::clog << "Weighted alpha complex with squared filtration values - Output_squared_values=true\n";
  Gudhi::Simplex_tree<> stree;
  if (alpha_complex_from_weighted_points.create_complex(stree)) {
    for (auto f_simplex : stree.filtration_simplex_range()) {
      std::clog << "   ( ";
      for (auto vertex : stree.simplex_vertex_range(f_simplex)) {
        std::clog << vertex << " ";
      }
      std::clog << ") -> " << "[" << stree.filtration(f_simplex) << "]\n";

      BOOST_CHECK(!std::isnan(stree.filtration(f_simplex)));
    }
  }
  std::clog << "Weighted alpha complex with square root filtration values - output_squared_values=false\n";
  Gudhi::Simplex_tree<> stree_sqrt;
  // set output_squared_values to false on a weighted alpha complex shall throw an exception
  BOOST_CHECK_THROW(alpha_complex_from_weighted_points.create_complex<false>(stree_sqrt), std::invalid_argument);
}
