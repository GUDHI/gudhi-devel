/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2015 Inria
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "alpha_complex_3d"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <cmath>  // float comparison
#include <limits>
#include <string>
#include <vector>
#include <random>
#include <cstddef>  // for std::size_t

#include <gudhi/Alpha_complex_3d.h>
#include <gudhi/graph_simplicial_complex.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Unitary_tests_utils.h>
// to construct Alpha_complex from a OFF file of points
#include <gudhi/Points_3D_off_io.h>

#include <CGAL/Random.h>
#include <CGAL/point_generators_3.h>

using Fast_weighted_alpha_complex_3d =
    Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::FAST, true, false>;
using Safe_weighted_alpha_complex_3d =
    Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::SAFE, true, false>;
using Exact_weighted_alpha_complex_3d =
    Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::EXACT, true, false>;

typedef boost::mpl::list<Fast_weighted_alpha_complex_3d, Safe_weighted_alpha_complex_3d,
                         Exact_weighted_alpha_complex_3d>
    weighted_variants_type_list;

#ifdef GUDHI_DEBUG
BOOST_AUTO_TEST_CASE_TEMPLATE(Alpha_complex_weighted_throw, Weighted_alpha_complex_3d, weighted_variants_type_list) {
  using Point_3 = typename Weighted_alpha_complex_3d::Point_3;
  std::vector<Point_3> w_points;
  w_points.push_back(Point_3(0.0, 0.0, 0.0));
  w_points.push_back(Point_3(0.0, 0.0, 0.2));
  w_points.push_back(Point_3(0.2, 0.0, 0.2));
  // w_points.push_back(Point_3(0.6, 0.6, 0.0));
  // w_points.push_back(Point_3(0.8, 0.8, 0.2));
  // w_points.push_back(Point_3(0.2, 0.8, 0.6));

  // weights size is different from w_points size to make weighted Alpha_complex_3d throw in debug mode
  std::vector<double> weights = {0.01, 0.005, 0.006, 0.01, 0.009, 0.001};

  std::cout << "Check exception throw in debug mode" << std::endl;
  BOOST_CHECK_THROW(Weighted_alpha_complex_3d wac(w_points, weights), std::invalid_argument);
}
#endif

BOOST_AUTO_TEST_CASE_TEMPLATE(Alpha_complex_weighted, Weighted_alpha_complex_3d, weighted_variants_type_list) {
  std::cout << "Weighted alpha complex 3d from points and weights" << std::endl;
  using Point_3 = typename Weighted_alpha_complex_3d::Point_3;
  std::vector<Point_3> w_points;
  w_points.push_back(Point_3(0.0, 0.0, 0.0));
  w_points.push_back(Point_3(0.0, 0.0, 0.2));
  w_points.push_back(Point_3(0.2, 0.0, 0.2));
  w_points.push_back(Point_3(0.6, 0.6, 0.0));
  w_points.push_back(Point_3(0.8, 0.8, 0.2));
  w_points.push_back(Point_3(0.2, 0.8, 0.6));

  // weights size is different from w_points size to make weighted Alpha_complex_3d throw in debug mode
  std::vector<double> weights = {0.01, 0.005, 0.006, 0.01, 0.009, 0.001};

  Weighted_alpha_complex_3d alpha_complex_p_a_w(w_points, weights);
  Gudhi::Simplex_tree<> stree;
  alpha_complex_p_a_w.create_complex(stree);

  std::cout << "Weighted alpha complex 3d from weighted points" << std::endl;
  using Weighted_point_3 = typename Weighted_alpha_complex_3d::Weighted_point_3;

  std::vector<Weighted_point_3> weighted_points;

  for (std::size_t i = 0; i < w_points.size(); i++) {
    weighted_points.push_back(Weighted_point_3(w_points[i], weights[i]));
  }
  Weighted_alpha_complex_3d alpha_complex_w_p(weighted_points);

  Gudhi::Simplex_tree<> stree_bis;
  alpha_complex_w_p.create_complex(stree_bis);

  // ---------------------
  // Compare both versions
  // ---------------------
  std::cout << "Weighted alpha complex 3d is of dimension " << stree_bis.dimension() << " - versus "
            << stree.dimension() << std::endl;
  BOOST_CHECK(stree_bis.dimension() == stree.dimension());
  std::cout << "Weighted alpha complex 3d num_simplices " << stree_bis.num_simplices() << " - versus "
            << stree.num_simplices() << std::endl;
  BOOST_CHECK(stree_bis.num_simplices() == stree.num_simplices());
  std::cout << "Weighted alpha complex 3d num_vertices " << stree_bis.num_vertices() << " - versus "
            << stree.num_vertices() << std::endl;
  BOOST_CHECK(stree_bis.num_vertices() == stree.num_vertices());

  auto sh = stree.filtration_simplex_range().begin();
  while (sh != stree.filtration_simplex_range().end()) {
    std::vector<int> simplex;
    std::vector<int> exact_simplex;
#ifdef DEBUG_TRACES
    std::cout << " ( ";
#endif
    for (auto vertex : stree.simplex_vertex_range(*sh)) {
      simplex.push_back(vertex);
#ifdef DEBUG_TRACES
      std::cout << vertex << " ";
#endif
    }
#ifdef DEBUG_TRACES
    std::cout << ") -> "
              << "[" << stree.filtration(*sh) << "] ";
    std::cout << std::endl;
#endif

    // Find it in the exact structure
    auto sh_exact = stree_bis.find(simplex);
    BOOST_CHECK(sh_exact != stree_bis.null_simplex());

    // Exact and non-exact version is not exactly the same due to float comparison
    GUDHI_TEST_FLOAT_EQUALITY_CHECK(stree_bis.filtration(sh_exact), stree.filtration(*sh));

    ++sh;
  }
}
