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

#include <gudhi/Alpha_complex_3d.h>
#include <gudhi/Alpha_complex_3d_options.h>
// to construct a simplex_tree from Delaunay_triangulation
#include <gudhi/graph_simplicial_complex.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Unitary_tests_utils.h>
#include <gudhi/Points_3D_off_io.h>

using Alpha_shapes_3d = Gudhi::alpha_complex::Alpha_shapes_3d;
using Exact_alpha_shapes_3d = Gudhi::alpha_complex::Exact_alpha_shapes_3d;
using Weighted_alpha_shapes_3d = Gudhi::alpha_complex::Weighted_alpha_shapes_3d;
using Periodic_alpha_shapes_3d = Gudhi::alpha_complex::Periodic_alpha_shapes_3d;
using Weighted_periodic_alpha_shapes_3d = Gudhi::alpha_complex::Weighted_periodic_alpha_shapes_3d;

BOOST_AUTO_TEST_CASE(Alpha_complex_3d_from_points) {
  // -----------------
  // Non exact version
  // -----------------
  std::cout << "Alpha complex 3d" << std::endl;
  std::vector<Alpha_shapes_3d::Point_3> points;
  points.push_back(Alpha_shapes_3d::Point_3(0.0, 0.0, 0.0));
  points.push_back(Alpha_shapes_3d::Point_3(0.0, 0.0, 0.2));
  points.push_back(Alpha_shapes_3d::Point_3(0.2, 0.0, 0.2));
  points.push_back(Alpha_shapes_3d::Point_3(0.6, 0.6, 0.0));
  points.push_back(Alpha_shapes_3d::Point_3(0.8, 0.8, 0.2));
  points.push_back(Alpha_shapes_3d::Point_3(0.2, 0.8, 0.6));

  Gudhi::alpha_complex::Alpha_complex_3d<Alpha_shapes_3d> alpha_complex(points);

  Gudhi::Simplex_tree<> stree;
  alpha_complex.create_complex(stree);

  // -----------------
  // Exact version
  // -----------------
  std::cout << "Exact alpha complex 3d" << std::endl;
  using Exact_alpha_shapes_3d = Gudhi::alpha_complex::Exact_alpha_shapes_3d;

  Gudhi::alpha_complex::Alpha_complex_3d<Exact_alpha_shapes_3d> exact_alpha_complex(points);

  Gudhi::Simplex_tree<> exact_stree;
  exact_alpha_complex.create_complex(exact_stree);

  // ---------------------
  // Compare both versions
  // ---------------------
  std::cout << "Exact Alpha complex 3d is of dimension " << exact_stree.dimension()
            << " - Non exact is " << stree.dimension() << std::endl;
  BOOST_CHECK(exact_stree.dimension() == stree.dimension());
  std::cout << "Exact Alpha complex 3d num_simplices " << exact_stree.num_simplices()
            << " - Non exact is " << stree.num_simplices() << std::endl;
  BOOST_CHECK(exact_stree.num_simplices() == stree.num_simplices());
  std::cout << "Exact Alpha complex 3d num_vertices " << exact_stree.num_vertices()
            << " - Non exact is " << stree.num_vertices() << std::endl;
  BOOST_CHECK(exact_stree.num_vertices() == stree.num_vertices());

  auto sh = stree.filtration_simplex_range().begin();
  auto sh_exact = exact_stree.filtration_simplex_range().begin();
  while(sh != stree.filtration_simplex_range().end() && sh_exact != exact_stree.filtration_simplex_range().end()) {
    std::vector<int> simplex;
    std::vector<int> exact_simplex;
    std::cout << "Non-exact ( ";
    for (auto vertex : stree.simplex_vertex_range(*sh)) {
      simplex.push_back(vertex);
      std::cout << vertex << " ";
    }
    std::cout << ") -> " << "[" << stree.filtration(*sh) << "] ";
    std::cout << std::endl;
    std::cout << "Exact     ( ";
    for (auto vertex : exact_stree.simplex_vertex_range(*sh_exact)) {
      exact_simplex.push_back(vertex);
      std::cout << vertex << " ";
    }
    std::cout << ") -> " << "[" << exact_stree.filtration(*sh_exact) << "] ";
    std::cout << std::endl;
    BOOST_CHECK(exact_simplex == simplex);

    // Exact and non-exact version is not exactly the same due to float comparison
    GUDHI_TEST_FLOAT_EQUALITY_CHECK(exact_stree.filtration(*sh_exact), stree.filtration(*sh));
    ++sh;
    ++sh_exact;
  }
}

#ifdef GUDHI_DEBUG
BOOST_AUTO_TEST_CASE(Alpha_complex_weighted_throw) {
  std::vector<Weighted_alpha_shapes_3d::Point_3> w_points;
  w_points.push_back(Weighted_alpha_shapes_3d::Point_3(0.0, 0.0, 0.0));
  w_points.push_back(Weighted_alpha_shapes_3d::Point_3(0.0, 0.0, 0.2));
  w_points.push_back(Weighted_alpha_shapes_3d::Point_3(0.2, 0.0, 0.2));
  // w_points.push_back(Weighted_alpha_shapes_3d::Point_3(0.6, 0.6, 0.0));
  // w_points.push_back(Weighted_alpha_shapes_3d::Point_3(0.8, 0.8, 0.2));
  // w_points.push_back(Weighted_alpha_shapes_3d::Point_3(0.2, 0.8, 0.6));

  // weights size is different from w_points size to make weighted Alpha_complex_3d throw in debug mode
  std::vector<double> weights = {0.01, 0.005, 0.006, 0.01, 0.009, 0.001};

  std::cout << "Check exception throw in debug mode" << std::endl;
  BOOST_CHECK_THROW (Gudhi::alpha_complex::Alpha_complex_3d<Weighted_alpha_shapes_3d> wac(w_points, weights),
                     std::invalid_argument);
}
#endif

BOOST_AUTO_TEST_CASE(Alpha_complex_weighted) {
  std::cout << "Weighted alpha complex 3d" << std::endl;
  using Weighted_alpha_shapes_3d = Gudhi::alpha_complex::Weighted_alpha_shapes_3d;
  std::vector<Weighted_alpha_shapes_3d::Point_3> w_points;
  w_points.push_back(Weighted_alpha_shapes_3d::Point_3(0.0, 0.0, 0.0));
  w_points.push_back(Weighted_alpha_shapes_3d::Point_3(0.0, 0.0, 0.2));
  w_points.push_back(Weighted_alpha_shapes_3d::Point_3(0.2, 0.0, 0.2));
  w_points.push_back(Weighted_alpha_shapes_3d::Point_3(0.6, 0.6, 0.0));
  w_points.push_back(Weighted_alpha_shapes_3d::Point_3(0.8, 0.8, 0.2));
  w_points.push_back(Weighted_alpha_shapes_3d::Point_3(0.2, 0.8, 0.6));

  // weights size is different from w_points size to make weighted Alpha_complex_3d throw in debug mode
  std::vector<double> weights = {0.01, 0.005, 0.006, 0.01, 0.009, 0.001};

  Gudhi::alpha_complex::Alpha_complex_3d<Weighted_alpha_shapes_3d> weighted_alpha_complex(w_points, weights);
  Gudhi::Simplex_tree<> w_stree;
  weighted_alpha_complex.create_complex(w_stree);

  std::cout << "Weighted Alpha complex 3d is of dimension " << w_stree.dimension() << std::endl;
  BOOST_CHECK(w_stree.dimension() == 3);
  std::cout << "    num_simplices " << w_stree.num_simplices() << std::endl;
  BOOST_CHECK(w_stree.num_simplices() == 35);
  std::cout << "    num_vertices " << w_stree.num_vertices() << std::endl;
  BOOST_CHECK(w_stree.num_vertices() == 6);
}

#ifdef GUDHI_DEBUG
BOOST_AUTO_TEST_CASE(Alpha_complex_periodic_throw) {
  std::cout << "Periodic alpha complex 3d exception throw" << std::endl;
  std::vector<Periodic_alpha_shapes_3d::Point_3> p_points;

  // Not important, this is not what we want to check
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.0, 0.0, 0.0));

  std::cout << "Check exception throw in debug mode" << std::endl;
  using Periodic_alpha_complex_3d = Gudhi::alpha_complex::Alpha_complex_3d<Periodic_alpha_shapes_3d>;
  // Check it throws an exception when the cuboid is not iso
  BOOST_CHECK_THROW (Periodic_alpha_complex_3d periodic_alpha_complex(p_points, 0., 0., 0., 0.9, 1., 1.),
                     std::invalid_argument);
  BOOST_CHECK_THROW (Periodic_alpha_complex_3d periodic_alpha_complex(p_points, 0., 0., 0., 1., 0.9, 1.),
                     std::invalid_argument);
  BOOST_CHECK_THROW (Periodic_alpha_complex_3d periodic_alpha_complex(p_points, 0., 0., 0., 1., 1., 0.9),
                     std::invalid_argument);

}
#endif

BOOST_AUTO_TEST_CASE(Alpha_complex_periodic) {
  std::cout << "Periodic alpha complex 3d" << std::endl;
  std::vector<Periodic_alpha_shapes_3d::Point_3> p_points;

  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.0, 0.0, 0.0));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.0, 0.0, 0.2));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.0, 0.0, 0.4));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.0, 0.0, 0.6));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.0, 0.0, 0.8));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.0, 0.2, 0.0));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.0, 0.2, 0.2));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.0, 0.2, 0.4));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.0, 0.2, 0.6));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.0, 0.2, 0.8));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.0, 0.4, 0.0));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.0, 0.4, 0.2));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.0, 0.4, 0.4));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.0, 0.4, 0.6));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.0, 0.4, 0.8));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.0, 0.6, 0.0));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.0, 0.6, 0.2));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.0, 0.6, 0.4));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.0, 0.6, 0.6));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.0, 0.6, 0.8));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.0, 0.8, 0.0));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.0, 0.8, 0.2));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.0, 0.8, 0.4));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.0, 0.8, 0.6));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.0, 0.8, 0.8));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.2, 0.0, 0.0));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.2, 0.0, 0.2));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.2, 0.0, 0.4));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.2, 0.0, 0.6));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.2, 0.0, 0.8));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.2, 0.2, 0.0));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.2, 0.2, 0.2));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.2, 0.2, 0.4));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.2, 0.2, 0.6));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.2, 0.2, 0.8));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.2, 0.4, 0.0));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.2, 0.4, 0.2));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.2, 0.4, 0.4));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.2, 0.4, 0.6));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.2, 0.4, 0.8));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.2, 0.6, 0.0));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.2, 0.6, 0.2));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.2, 0.6, 0.4));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.2, 0.6, 0.6));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.2, 0.6, 0.8));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.2, 0.8, 0.0));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.2, 0.8, 0.2));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.2, 0.8, 0.4));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.2, 0.8, 0.6));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.2, 0.8, 0.8));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.4, 0.0, 0.0));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.4, 0.0, 0.2));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.4, 0.0, 0.4));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.4, 0.0, 0.6));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.4, 0.0, 0.8));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.4, 0.2, 0.0));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.4, 0.2, 0.2));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.4, 0.2, 0.4));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.4, 0.2, 0.6));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.4, 0.2, 0.8));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.4, 0.4, 0.0));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.4, 0.4, 0.2));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.4, 0.4, 0.4));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.4, 0.4, 0.6));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.4, 0.4, 0.8));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.4, 0.6, 0.0));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.4, 0.6, 0.2));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.4, 0.6, 0.4));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.4, 0.6, 0.6));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.4, 0.6, 0.8));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.4, 0.8, 0.0));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.4, 0.8, 0.2));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.4, 0.8, 0.4));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.4, 0.8, 0.6));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.4, 0.8, 0.8));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.6, 0.0, 0.0));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.6, 0.0, 0.2));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.6, 0.0, 0.4));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.6, 0.0, 0.6));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.6, 0.0, 0.8));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.6, 0.1, 0.0));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.6, 0.2, 0.0));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.6, 0.2, 0.2));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.6, 0.2, 0.4));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.6, 0.2, 0.6));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.6, 0.2, 0.8));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.6, 0.4, 0.0));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.6, 0.4, 0.2));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.6, 0.4, 0.4));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.6, 0.4, 0.6));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.6, 0.4, 0.8));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.6, 0.6, 0.0));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.6, 0.6, 0.2));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.6, 0.6, 0.4));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.6, 0.6, 0.6));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.6, 0.6, 0.8));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.6, 0.8, 0.0));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.6, 0.8, 0.2));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.6, 0.8, 0.4));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.6, 0.8, 0.6));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.6, 0.8, 0.8));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.8, 0.0, 0.0));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.8, 0.0, 0.2));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.8, 0.0, 0.4));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.8, 0.0, 0.6));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.8, 0.0, 0.8));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.8, 0.2, 0.0));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.8, 0.2, 0.2));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.8, 0.2, 0.4));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.8, 0.2, 0.6));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.8, 0.2, 0.8));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.8, 0.4, 0.0));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.8, 0.4, 0.2));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.8, 0.4, 0.4));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.8, 0.4, 0.6));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.8, 0.4, 0.8));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.8, 0.6, 0.0));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.8, 0.6, 0.2));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.8, 0.6, 0.4));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.8, 0.6, 0.6));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.8, 0.6, 0.8));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.8, 0.8, 0.0));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.8, 0.8, 0.2));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.8, 0.8, 0.4));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.8, 0.8, 0.6));

  Gudhi::alpha_complex::Alpha_complex_3d<Periodic_alpha_shapes_3d> periodic_alpha_complex(p_points,
                                                                                          0., 0., 0.,
                                                                                          1., 1., 1.);

  Gudhi::Simplex_tree<> p_stree;
  periodic_alpha_complex.create_complex(p_stree);

  std::cout << "Periodic Alpha complex 3d is of dimension " << p_stree.dimension() << std::endl;
  BOOST_CHECK(p_stree.dimension() == 3);
  std::cout << "    num_simplices " << p_stree.num_simplices() << std::endl;
  BOOST_CHECK(p_stree.num_simplices() == 3266);
  std::cout << "    num_vertices " << p_stree.num_vertices() << std::endl;
  BOOST_CHECK(p_stree.num_vertices() == 125);

}

#ifdef GUDHI_DEBUG
BOOST_AUTO_TEST_CASE(Alpha_complex_weighted_periodic_throw) {
  std::cout << "Weighted periodic alpha complex 3d exception throw" << std::endl;

  std::vector<Weighted_periodic_alpha_shapes_3d::Point_3> wp_points;
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.0, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.0, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.0, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.0, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.0, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.2, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.2, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.2, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.2, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.2, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.4, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.4, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.4, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.4, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.4, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.6, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.6, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.6, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.6, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.6, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.8, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.8, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.8, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.8, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.8, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.0, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.0, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.0, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.0, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.0, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.2, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.2, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.2, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.2, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.2, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.4, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.4, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.4, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.4, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.4, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.6, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.6, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.6, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.6, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.6, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.8, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.8, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.8, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.8, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.8, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.0, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.0, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.0, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.0, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.0, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.2, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.2, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.2, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.2, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.2, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.4, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.4, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.4, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.4, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.4, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.6, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.6, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.6, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.6, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.6, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.8, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.8, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.8, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.8, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.8, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.0, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.0, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.0, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.0, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.0, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.1, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.2, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.2, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.2, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.2, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.2, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.4, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.4, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.4, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.4, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.4, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.6, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.6, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.6, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.6, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.6, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.8, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.8, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.8, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.8, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.8, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.0, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.0, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.0, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.0, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.0, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.2, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.2, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.2, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.2, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.2, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.4, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.4, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.4, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.4, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.4, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.6, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.6, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.6, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.6, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.6, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.8, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.8, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.8, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.8, 0.6));

  std::vector<double> p_weights;

  std::random_device rd;
  std::mt19937 mt(rd());
  // Weights must be in range [0, <1/64]
  std::uniform_real_distribution<double> dist(0.0, 0.0156245);

  for (std::size_t i = 0; i < wp_points.size(); ++i) {
    double value = dist(mt);
    p_weights.push_back(value);
  }

  std::cout << "Cuboid is not iso exception" << std::endl;
  using Weighted_periodic_alpha_complex_3d = Gudhi::alpha_complex::Alpha_complex_3d<Weighted_periodic_alpha_shapes_3d>;
  // Check it throws an exception when the cuboid is not iso
  BOOST_CHECK_THROW (Weighted_periodic_alpha_complex_3d periodic_alpha_complex(wp_points, p_weights, 0., 0., 0., 0.9, 1., 1.),
                     std::invalid_argument);
  BOOST_CHECK_THROW (Weighted_periodic_alpha_complex_3d periodic_alpha_complex(wp_points, p_weights, 0., 0., 0., 1., 0.9, 1.),
                     std::invalid_argument);
  BOOST_CHECK_THROW (Weighted_periodic_alpha_complex_3d periodic_alpha_complex(wp_points, p_weights, 0., 0., 0., 1., 1., 0.9),
                     std::invalid_argument);

  std::cout << "0 <= point.weight() < 1/64 * domain_size * domain_size exception" << std::endl;
  // Weights must be in range [0, <1/64]
  p_weights[25] = 1.0;
  BOOST_CHECK_THROW (Weighted_periodic_alpha_complex_3d periodic_alpha_complex(wp_points, p_weights, 0., 0., 0., 1., 1., 1.),
                     std::invalid_argument);
  // Weights must be in range [0, <1/64]
  p_weights[25] = 0.012;
  p_weights[14] = -0.012;
  BOOST_CHECK_THROW (Weighted_periodic_alpha_complex_3d periodic_alpha_complex(wp_points, p_weights, 0., 0., 0., 1., 1., 1.),
                     std::invalid_argument);
  p_weights[14] = 0.005;

  std::cout << "wp_points and p_weights size exception" << std::endl;
  // Weights and points must have the same size
  // + 1
  p_weights.push_back(0.007);
  BOOST_CHECK_THROW (Weighted_periodic_alpha_complex_3d periodic_alpha_complex(wp_points, p_weights, 0., 0., 0., 1., 1., 1.),
                     std::invalid_argument);
  // - 1
  p_weights.pop_back();
  p_weights.pop_back();
  BOOST_CHECK_THROW (Weighted_periodic_alpha_complex_3d periodic_alpha_complex(wp_points, p_weights, 0., 0., 0., 1., 1., 1.),
                     std::invalid_argument);
}
#endif

BOOST_AUTO_TEST_CASE(Alpha_complex_weighted_periodic) {
  std::cout << "Weighted periodic alpha complex 3d" << std::endl;

  std::vector<Weighted_periodic_alpha_shapes_3d::Point_3> wp_points;
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.0, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.0, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.0, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.0, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.0, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.2, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.2, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.2, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.2, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.2, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.4, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.4, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.4, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.4, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.4, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.6, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.6, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.6, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.6, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.6, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.8, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.8, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.8, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.8, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.8, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.0, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.0, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.0, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.0, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.0, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.2, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.2, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.2, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.2, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.2, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.4, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.4, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.4, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.4, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.4, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.6, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.6, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.6, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.6, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.6, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.8, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.8, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.8, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.8, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.8, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.0, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.0, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.0, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.0, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.0, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.2, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.2, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.2, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.2, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.2, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.4, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.4, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.4, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.4, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.4, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.6, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.6, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.6, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.6, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.6, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.8, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.8, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.8, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.8, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.8, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.0, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.0, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.0, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.0, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.0, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.1, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.2, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.2, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.2, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.2, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.2, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.4, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.4, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.4, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.4, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.4, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.6, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.6, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.6, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.6, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.6, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.8, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.8, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.8, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.8, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.8, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.0, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.0, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.0, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.0, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.0, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.2, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.2, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.2, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.2, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.2, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.4, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.4, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.4, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.4, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.4, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.6, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.6, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.6, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.6, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.6, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.8, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.8, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.8, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.8, 0.6));

  std::vector<double> p_weights;

  std::random_device rd;
  std::mt19937 mt(rd());
  // Weights must be in range [0, <1/64]
  std::uniform_real_distribution<double> dist(0.01, 0.0156245);

  for (std::size_t i = 0; i < wp_points.size(); ++i) {
    double value = dist(mt);
    p_weights.push_back(value);
  }

  using Weighted_periodic_alpha_complex_3d = Gudhi::alpha_complex::Alpha_complex_3d<Weighted_periodic_alpha_shapes_3d>;
  Weighted_periodic_alpha_complex_3d weighted_periodic_alpha_complex(wp_points, p_weights, 0., 0., 0., 1., 1., 1.);

  Gudhi::Simplex_tree<> wp_stree;
  weighted_periodic_alpha_complex.create_complex(wp_stree);

  std::cout << "Weighted periodic Alpha complex 3d is of dimension " << wp_stree.dimension() << std::endl;
  BOOST_CHECK(wp_stree.dimension() == 3);
  std::cout << "    num_simplices " << wp_stree.num_simplices() << std::endl;
  BOOST_CHECK(wp_stree.num_simplices() >= 3100);
  std::cout << "    num_vertices " << wp_stree.num_vertices() << std::endl;
  BOOST_CHECK(wp_stree.num_vertices() == 125);
}
