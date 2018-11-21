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

using Fast_alpha_complex_3d =
    Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::FAST, false, false>;
using Safe_alpha_complex_3d =
    Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::SAFE, false, false>;
using Exact_alpha_complex_3d =
    Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::EXACT, false, false>;

BOOST_AUTO_TEST_CASE(Alpha_complex_3d_from_points) {
  // -----------------
  // Fast version
  // -----------------
  std::cout << "Fast alpha complex 3d" << std::endl;
  std::vector<Fast_alpha_complex_3d::Point_3> points;
  points.push_back(Fast_alpha_complex_3d::Point_3(0.0, 0.0, 0.0));
  points.push_back(Fast_alpha_complex_3d::Point_3(0.0, 0.0, 0.2));
  points.push_back(Fast_alpha_complex_3d::Point_3(0.2, 0.0, 0.2));
  points.push_back(Fast_alpha_complex_3d::Point_3(0.6, 0.6, 0.0));
  points.push_back(Fast_alpha_complex_3d::Point_3(0.8, 0.8, 0.2));
  points.push_back(Fast_alpha_complex_3d::Point_3(0.2, 0.8, 0.6));

  Fast_alpha_complex_3d alpha_complex(points);

  Gudhi::Simplex_tree<> stree;
  alpha_complex.create_complex(stree);

  // -----------------
  // Exact version
  // -----------------
  std::cout << "Exact alpha complex 3d" << std::endl;

  Exact_alpha_complex_3d exact_alpha_complex(points);

  Gudhi::Simplex_tree<> exact_stree;
  exact_alpha_complex.create_complex(exact_stree);

  // ---------------------
  // Compare both versions
  // ---------------------
  std::cout << "Exact Alpha complex 3d is of dimension " << exact_stree.dimension() << " - Non exact is "
            << stree.dimension() << std::endl;
  BOOST_CHECK(exact_stree.dimension() == stree.dimension());
  std::cout << "Exact Alpha complex 3d num_simplices " << exact_stree.num_simplices() << " - Non exact is "
            << stree.num_simplices() << std::endl;
  BOOST_CHECK(exact_stree.num_simplices() == stree.num_simplices());
  std::cout << "Exact Alpha complex 3d num_vertices " << exact_stree.num_vertices() << " - Non exact is "
            << stree.num_vertices() << std::endl;
  BOOST_CHECK(exact_stree.num_vertices() == stree.num_vertices());

  auto sh = stree.filtration_simplex_range().begin();
  while (sh != stree.filtration_simplex_range().end()) {
    std::vector<int> simplex;
    std::vector<int> exact_simplex;
    std::cout << "Non-exact ( ";
    for (auto vertex : stree.simplex_vertex_range(*sh)) {
      simplex.push_back(vertex);
      std::cout << vertex << " ";
    }
    std::cout << ") -> "
              << "[" << stree.filtration(*sh) << "] ";
    std::cout << std::endl;

    // Find it in the exact structure
    auto sh_exact = exact_stree.find(simplex);
    BOOST_CHECK(sh_exact != exact_stree.null_simplex());

    // Exact and non-exact version is not exactly the same due to float comparison
    GUDHI_TEST_FLOAT_EQUALITY_CHECK(exact_stree.filtration(sh_exact), stree.filtration(*sh));

    ++sh;
  }
  // -----------------
  // Safe version
  // -----------------
  std::cout << "Safe alpha complex 3d" << std::endl;

  Safe_alpha_complex_3d safe_alpha_complex(points);

  Gudhi::Simplex_tree<> safe_stree;
  safe_alpha_complex.create_complex(safe_stree);

  // ---------------------
  // Compare both versions
  // ---------------------
  std::cout << "Exact Alpha complex 3d is of dimension " << safe_stree.dimension() << " - Non exact is "
            << stree.dimension() << std::endl;
  BOOST_CHECK(safe_stree.dimension() == stree.dimension());
  std::cout << "Exact Alpha complex 3d num_simplices " << safe_stree.num_simplices() << " - Non exact is "
            << stree.num_simplices() << std::endl;
  BOOST_CHECK(safe_stree.num_simplices() == stree.num_simplices());
  std::cout << "Exact Alpha complex 3d num_vertices " << safe_stree.num_vertices() << " - Non exact is "
            << stree.num_vertices() << std::endl;
  BOOST_CHECK(safe_stree.num_vertices() == stree.num_vertices());

  auto safe_sh = stree.filtration_simplex_range().begin();
  while (safe_sh != stree.filtration_simplex_range().end()) {
    std::vector<int> simplex;
    std::vector<int> exact_simplex;
#ifdef DEBUG_TRACES
    std::cout << "Non-exact ( ";
#endif
    for (auto vertex : stree.simplex_vertex_range(*safe_sh)) {
      simplex.push_back(vertex);
#ifdef DEBUG_TRACES
      std::cout << vertex << " ";
#endif
    }
#ifdef DEBUG_TRACES
    std::cout << ") -> "
              << "[" << stree.filtration(*safe_sh) << "] ";
    std::cout << std::endl;
#endif

    // Find it in the exact structure
    auto sh_exact = safe_stree.find(simplex);
    BOOST_CHECK(sh_exact != safe_stree.null_simplex());

    // Exact and non-exact version is not exactly the same due to float comparison
    GUDHI_TEST_FLOAT_EQUALITY_CHECK(safe_stree.filtration(sh_exact), stree.filtration(*safe_sh));

    ++safe_sh;
  }
}
