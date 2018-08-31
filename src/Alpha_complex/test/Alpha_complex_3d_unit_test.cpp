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
#include <gudhi/graph_simplicial_complex.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Unitary_tests_utils.h>
#include <gudhi/Points_3D_off_io.h>

using Fast_alpha_complex_3d = Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::fast, false, false>;
using Exact_alpha_complex_3d = Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::exact, false, false>;
using Fast_weighted_alpha_complex_3d = Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::fast, true, false>;
using Exact_weighted_alpha_complex_3d = Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::exact, true, false>;
using Fast_periodic_alpha_complex_3d = Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::fast, false, true>;
using Exact_periodic_alpha_complex_3d = Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::exact, false, true>;
using Fast_weighted_periodic_alpha_complex_3d = Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::fast, true, true>;
using Exact_weighted_periodic_alpha_complex_3d = Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::exact, true, true>;

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
  while(sh != stree.filtration_simplex_range().end()) {
    std::vector<int> simplex;
    std::vector<int> exact_simplex;
    std::cout << "Non-exact ( ";
    for (auto vertex : stree.simplex_vertex_range(*sh)) {
      simplex.push_back(vertex);
      std::cout << vertex << " ";
    }
    std::cout << ") -> " << "[" << stree.filtration(*sh) << "] ";
    std::cout << std::endl;

    // Find it in the exact structure
    auto sh_exact = exact_stree.find(simplex);
    BOOST_CHECK(sh_exact != exact_stree.null_simplex());

    // Exact and non-exact version is not exactly the same due to float comparison
    GUDHI_TEST_FLOAT_EQUALITY_CHECK(exact_stree.filtration(sh_exact), stree.filtration(*sh));

    ++sh;
  }
}

#ifdef GUDHI_DEBUG
typedef boost::mpl::list<Fast_weighted_alpha_complex_3d, Exact_weighted_alpha_complex_3d> weighted_variants_type_list;

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
  BOOST_CHECK_THROW (Weighted_alpha_complex_3d wac(w_points, weights), std::invalid_argument);
}
#endif

BOOST_AUTO_TEST_CASE(Alpha_complex_weighted) {
  // ---------------------
  // Fast weighted version
  // ---------------------
  std::cout << "Fast weighted alpha complex 3d" << std::endl;
  std::vector<Fast_weighted_alpha_complex_3d::Point_3> w_points;
  w_points.push_back(Fast_weighted_alpha_complex_3d::Point_3(0.0, 0.0, 0.0));
  w_points.push_back(Fast_weighted_alpha_complex_3d::Point_3(0.0, 0.0, 0.2));
  w_points.push_back(Fast_weighted_alpha_complex_3d::Point_3(0.2, 0.0, 0.2));
  w_points.push_back(Fast_weighted_alpha_complex_3d::Point_3(0.6, 0.6, 0.0));
  w_points.push_back(Fast_weighted_alpha_complex_3d::Point_3(0.8, 0.8, 0.2));
  w_points.push_back(Fast_weighted_alpha_complex_3d::Point_3(0.2, 0.8, 0.6));

  // weights size is different from w_points size to make weighted Alpha_complex_3d throw in debug mode
  std::vector<double> weights = {0.01, 0.005, 0.006, 0.01, 0.009, 0.001};

  Fast_weighted_alpha_complex_3d weighted_alpha_complex(w_points, weights);
  Gudhi::Simplex_tree<> stree;
  weighted_alpha_complex.create_complex(stree);

  // ----------------------
  // Exact weighted version
  // ----------------------
  std::cout << "Exact weighted alpha complex 3d" << std::endl;

  std::vector<Exact_weighted_alpha_complex_3d::Point_3> e_w_points;
  e_w_points.push_back(Exact_weighted_alpha_complex_3d::Point_3(0.0, 0.0, 0.0));
  e_w_points.push_back(Exact_weighted_alpha_complex_3d::Point_3(0.0, 0.0, 0.2));
  e_w_points.push_back(Exact_weighted_alpha_complex_3d::Point_3(0.2, 0.0, 0.2));
  e_w_points.push_back(Exact_weighted_alpha_complex_3d::Point_3(0.6, 0.6, 0.0));
  e_w_points.push_back(Exact_weighted_alpha_complex_3d::Point_3(0.8, 0.8, 0.2));
  e_w_points.push_back(Exact_weighted_alpha_complex_3d::Point_3(0.2, 0.8, 0.6));
  Exact_weighted_alpha_complex_3d exact_alpha_complex(e_w_points, weights);

  Gudhi::Simplex_tree<> exact_stree;
  exact_alpha_complex.create_complex(exact_stree);

  // ---------------------
  // Compare both versions
  // ---------------------
  std::cout << "Exact weighted alpha complex 3d is of dimension " << exact_stree.dimension()
            << " - Non exact is " << stree.dimension() << std::endl;
  BOOST_CHECK(exact_stree.dimension() == stree.dimension());
  std::cout << "Exact weighted alpha complex 3d num_simplices " << exact_stree.num_simplices()
            << " - Non exact is " << stree.num_simplices() << std::endl;
  BOOST_CHECK(exact_stree.num_simplices() == stree.num_simplices());
  std::cout << "Exact weighted alpha complex 3d num_vertices " << exact_stree.num_vertices()
            << " - Non exact is " << stree.num_vertices() << std::endl;
  BOOST_CHECK(exact_stree.num_vertices() == stree.num_vertices());

  auto sh = stree.filtration_simplex_range().begin();
  while(sh != stree.filtration_simplex_range().end()) {
    std::vector<int> simplex;
    std::vector<int> exact_simplex;
    std::cout << "Non-exact ( ";
    for (auto vertex : stree.simplex_vertex_range(*sh)) {
      simplex.push_back(vertex);
      std::cout << vertex << " ";
    }
    std::cout << ") -> " << "[" << stree.filtration(*sh) << "] ";
    std::cout << std::endl;

    // Find it in the exact structure
    auto sh_exact = exact_stree.find(simplex);
    BOOST_CHECK(sh_exact != exact_stree.null_simplex());

    // Exact and non-exact version is not exactly the same due to float comparison
    GUDHI_TEST_FLOAT_EQUALITY_CHECK(exact_stree.filtration(sh_exact), stree.filtration(*sh));

    ++sh;
  }

}

#ifdef GUDHI_DEBUG
typedef boost::mpl::list<Fast_periodic_alpha_complex_3d, Fast_periodic_alpha_complex_3d> periodic_variants_type_list;

BOOST_AUTO_TEST_CASE_TEMPLATE(Alpha_complex_periodic_throw, Periodic_alpha_complex_3d, periodic_variants_type_list) {
  std::cout << "Periodic alpha complex 3d exception throw" << std::endl;
  using Point_3 = typename Periodic_alpha_complex_3d::Point_3;
  std::vector<Point_3> p_points;

  // Not important, this is not what we want to check
  p_points.push_back(Point_3(0.0, 0.0, 0.0));

  std::cout << "Check exception throw in debug mode" << std::endl;
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
  // ---------------------
  // Fast periodic version
  // ---------------------
  std::cout << "Fast periodic alpha complex 3d" << std::endl;
  std::vector<Fast_periodic_alpha_complex_3d::Point_3> p_points;

  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.0, 0.0, 0.0));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.0, 0.0, 0.2));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.0, 0.0, 0.4));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.0, 0.0, 0.6));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.0, 0.0, 0.8));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.0, 0.2, 0.0));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.0, 0.2, 0.2));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.0, 0.2, 0.4));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.0, 0.2, 0.6));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.0, 0.2, 0.8));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.0, 0.4, 0.0));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.0, 0.4, 0.2));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.0, 0.4, 0.4));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.0, 0.4, 0.6));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.0, 0.4, 0.8));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.0, 0.6, 0.0));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.0, 0.6, 0.2));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.0, 0.6, 0.4));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.0, 0.6, 0.6));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.0, 0.6, 0.8));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.0, 0.8, 0.0));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.0, 0.8, 0.2));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.0, 0.8, 0.4));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.0, 0.8, 0.6));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.0, 0.8, 0.8));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.2, 0.0, 0.0));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.2, 0.0, 0.2));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.2, 0.0, 0.4));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.2, 0.0, 0.6));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.2, 0.0, 0.8));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.2, 0.2, 0.0));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.2, 0.2, 0.2));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.2, 0.2, 0.4));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.2, 0.2, 0.6));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.2, 0.2, 0.8));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.2, 0.4, 0.0));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.2, 0.4, 0.2));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.2, 0.4, 0.4));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.2, 0.4, 0.6));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.2, 0.4, 0.8));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.2, 0.6, 0.0));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.2, 0.6, 0.2));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.2, 0.6, 0.4));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.2, 0.6, 0.6));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.2, 0.6, 0.8));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.2, 0.8, 0.0));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.2, 0.8, 0.2));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.2, 0.8, 0.4));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.2, 0.8, 0.6));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.2, 0.8, 0.8));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.4, 0.0, 0.0));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.4, 0.0, 0.2));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.4, 0.0, 0.4));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.4, 0.0, 0.6));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.4, 0.0, 0.8));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.4, 0.2, 0.0));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.4, 0.2, 0.2));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.4, 0.2, 0.4));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.4, 0.2, 0.6));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.4, 0.2, 0.8));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.4, 0.4, 0.0));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.4, 0.4, 0.2));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.4, 0.4, 0.4));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.4, 0.4, 0.6));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.4, 0.4, 0.8));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.4, 0.6, 0.0));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.4, 0.6, 0.2));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.4, 0.6, 0.4));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.4, 0.6, 0.6));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.4, 0.6, 0.8));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.4, 0.8, 0.0));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.4, 0.8, 0.2));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.4, 0.8, 0.4));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.4, 0.8, 0.6));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.4, 0.8, 0.8));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.6, 0.0, 0.0));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.6, 0.0, 0.2));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.6, 0.0, 0.4));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.6, 0.0, 0.6));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.6, 0.0, 0.8));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.6, 0.1, 0.0));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.6, 0.2, 0.0));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.6, 0.2, 0.2));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.6, 0.2, 0.4));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.6, 0.2, 0.6));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.6, 0.2, 0.8));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.6, 0.4, 0.0));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.6, 0.4, 0.2));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.6, 0.4, 0.4));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.6, 0.4, 0.6));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.6, 0.4, 0.8));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.6, 0.6, 0.0));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.6, 0.6, 0.2));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.6, 0.6, 0.4));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.6, 0.6, 0.6));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.6, 0.6, 0.8));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.6, 0.8, 0.0));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.6, 0.8, 0.2));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.6, 0.8, 0.4));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.6, 0.8, 0.6));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.6, 0.8, 0.8));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.8, 0.0, 0.0));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.8, 0.0, 0.2));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.8, 0.0, 0.4));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.8, 0.0, 0.6));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.8, 0.0, 0.8));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.8, 0.2, 0.0));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.8, 0.2, 0.2));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.8, 0.2, 0.4));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.8, 0.2, 0.6));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.8, 0.2, 0.8));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.8, 0.4, 0.0));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.8, 0.4, 0.2));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.8, 0.4, 0.4));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.8, 0.4, 0.6));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.8, 0.4, 0.8));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.8, 0.6, 0.0));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.8, 0.6, 0.2));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.8, 0.6, 0.4));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.8, 0.6, 0.6));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.8, 0.6, 0.8));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.8, 0.8, 0.0));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.8, 0.8, 0.2));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.8, 0.8, 0.4));
  p_points.push_back(Fast_periodic_alpha_complex_3d::Point_3(0.8, 0.8, 0.6));

  Fast_periodic_alpha_complex_3d periodic_alpha_complex(p_points, 0., 0., 0., 1., 1., 1.);

  Gudhi::Simplex_tree<> stree;
  periodic_alpha_complex.create_complex(stree);

  // ----------------------
  // Exact periodic version
  // ----------------------
  std::cout << "Exact periodic alpha complex 3d" << std::endl;

  std::vector<Exact_periodic_alpha_complex_3d::Point_3> e_p_points;

  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.0, 0.0, 0.0));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.0, 0.0, 0.2));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.0, 0.0, 0.4));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.0, 0.0, 0.6));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.0, 0.0, 0.8));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.0, 0.2, 0.0));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.0, 0.2, 0.2));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.0, 0.2, 0.4));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.0, 0.2, 0.6));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.0, 0.2, 0.8));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.0, 0.4, 0.0));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.0, 0.4, 0.2));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.0, 0.4, 0.4));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.0, 0.4, 0.6));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.0, 0.4, 0.8));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.0, 0.6, 0.0));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.0, 0.6, 0.2));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.0, 0.6, 0.4));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.0, 0.6, 0.6));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.0, 0.6, 0.8));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.0, 0.8, 0.0));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.0, 0.8, 0.2));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.0, 0.8, 0.4));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.0, 0.8, 0.6));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.0, 0.8, 0.8));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.2, 0.0, 0.0));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.2, 0.0, 0.2));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.2, 0.0, 0.4));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.2, 0.0, 0.6));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.2, 0.0, 0.8));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.2, 0.2, 0.0));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.2, 0.2, 0.2));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.2, 0.2, 0.4));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.2, 0.2, 0.6));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.2, 0.2, 0.8));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.2, 0.4, 0.0));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.2, 0.4, 0.2));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.2, 0.4, 0.4));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.2, 0.4, 0.6));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.2, 0.4, 0.8));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.2, 0.6, 0.0));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.2, 0.6, 0.2));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.2, 0.6, 0.4));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.2, 0.6, 0.6));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.2, 0.6, 0.8));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.2, 0.8, 0.0));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.2, 0.8, 0.2));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.2, 0.8, 0.4));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.2, 0.8, 0.6));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.2, 0.8, 0.8));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.4, 0.0, 0.0));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.4, 0.0, 0.2));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.4, 0.0, 0.4));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.4, 0.0, 0.6));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.4, 0.0, 0.8));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.4, 0.2, 0.0));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.4, 0.2, 0.2));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.4, 0.2, 0.4));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.4, 0.2, 0.6));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.4, 0.2, 0.8));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.4, 0.4, 0.0));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.4, 0.4, 0.2));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.4, 0.4, 0.4));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.4, 0.4, 0.6));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.4, 0.4, 0.8));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.4, 0.6, 0.0));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.4, 0.6, 0.2));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.4, 0.6, 0.4));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.4, 0.6, 0.6));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.4, 0.6, 0.8));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.4, 0.8, 0.0));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.4, 0.8, 0.2));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.4, 0.8, 0.4));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.4, 0.8, 0.6));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.4, 0.8, 0.8));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.6, 0.0, 0.0));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.6, 0.0, 0.2));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.6, 0.0, 0.4));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.6, 0.0, 0.6));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.6, 0.0, 0.8));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.6, 0.1, 0.0));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.6, 0.2, 0.0));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.6, 0.2, 0.2));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.6, 0.2, 0.4));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.6, 0.2, 0.6));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.6, 0.2, 0.8));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.6, 0.4, 0.0));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.6, 0.4, 0.2));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.6, 0.4, 0.4));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.6, 0.4, 0.6));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.6, 0.4, 0.8));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.6, 0.6, 0.0));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.6, 0.6, 0.2));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.6, 0.6, 0.4));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.6, 0.6, 0.6));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.6, 0.6, 0.8));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.6, 0.8, 0.0));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.6, 0.8, 0.2));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.6, 0.8, 0.4));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.6, 0.8, 0.6));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.6, 0.8, 0.8));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.8, 0.0, 0.0));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.8, 0.0, 0.2));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.8, 0.0, 0.4));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.8, 0.0, 0.6));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.8, 0.0, 0.8));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.8, 0.2, 0.0));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.8, 0.2, 0.2));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.8, 0.2, 0.4));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.8, 0.2, 0.6));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.8, 0.2, 0.8));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.8, 0.4, 0.0));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.8, 0.4, 0.2));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.8, 0.4, 0.4));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.8, 0.4, 0.6));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.8, 0.4, 0.8));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.8, 0.6, 0.0));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.8, 0.6, 0.2));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.8, 0.6, 0.4));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.8, 0.6, 0.6));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.8, 0.6, 0.8));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.8, 0.8, 0.0));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.8, 0.8, 0.2));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.8, 0.8, 0.4));
  e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(0.8, 0.8, 0.6));

  Exact_periodic_alpha_complex_3d exact_alpha_complex(e_p_points, 0., 0., 0., 1., 1., 1.);

  Gudhi::Simplex_tree<> exact_stree;
  exact_alpha_complex.create_complex(exact_stree);

  // ---------------------
  // Compare both versions
  // ---------------------
  std::cout << "Exact periodic alpha complex 3d is of dimension " << exact_stree.dimension()
            << " - Non exact is " << stree.dimension() << std::endl;
  BOOST_CHECK(exact_stree.dimension() == stree.dimension());
  std::cout << "Exact periodic alpha complex 3d num_simplices " << exact_stree.num_simplices()
            << " - Non exact is " << stree.num_simplices() << std::endl;
  BOOST_CHECK(exact_stree.num_simplices() == stree.num_simplices());
  std::cout << "Exact periodic alpha complex 3d num_vertices " << exact_stree.num_vertices()
            << " - Non exact is " << stree.num_vertices() << std::endl;
  BOOST_CHECK(exact_stree.num_vertices() == stree.num_vertices());

  auto sh = stree.filtration_simplex_range().begin();
  while(sh != stree.filtration_simplex_range().end()) {
    std::vector<int> simplex;
    std::vector<int> exact_simplex;
    std::cout << "Non-exact ( ";
    for (auto vertex : stree.simplex_vertex_range(*sh)) {
      simplex.push_back(vertex);
      std::cout << vertex << " ";
    }
    std::cout << ") -> " << "[" << stree.filtration(*sh) << "] ";
    std::cout << std::endl;

    // Find it in the exact structure
    auto sh_exact = exact_stree.find(simplex);
    // TODO(VR): BOOST_CHECK(sh_exact != exact_stree.null_simplex());

    // Exact and non-exact version is not exactly the same due to float comparison
    // TODO(VR): GUDHI_TEST_FLOAT_EQUALITY_CHECK(exact_stree.filtration(sh_exact), stree.filtration(*sh));
    ++sh;
  }


}

#ifdef GUDHI_DEBUG
typedef boost::mpl::list<Fast_weighted_periodic_alpha_complex_3d, Fast_weighted_periodic_alpha_complex_3d> wp_variants_type_list;

BOOST_AUTO_TEST_CASE_TEMPLATE(Alpha_complex_weighted_periodic_throw, Weighted_periodic_alpha_complex_3d, wp_variants_type_list) {
  std::cout << "Weighted periodic alpha complex 3d exception throw" << std::endl;

  using Point_3 = typename Weighted_periodic_alpha_complex_3d::Point_3;
  std::vector<Point_3> wp_points;
  wp_points.push_back(Point_3(0.0, 0.0, 0.0));
  wp_points.push_back(Point_3(0.0, 0.0, 0.2));
  wp_points.push_back(Point_3(0.0, 0.0, 0.4));
  wp_points.push_back(Point_3(0.0, 0.0, 0.6));
  wp_points.push_back(Point_3(0.0, 0.0, 0.8));
  wp_points.push_back(Point_3(0.0, 0.2, 0.0));
  wp_points.push_back(Point_3(0.0, 0.2, 0.2));
  wp_points.push_back(Point_3(0.0, 0.2, 0.4));
  wp_points.push_back(Point_3(0.0, 0.2, 0.6));
  wp_points.push_back(Point_3(0.0, 0.2, 0.8));
  wp_points.push_back(Point_3(0.0, 0.4, 0.0));
  wp_points.push_back(Point_3(0.0, 0.4, 0.2));
  wp_points.push_back(Point_3(0.0, 0.4, 0.4));
  wp_points.push_back(Point_3(0.0, 0.4, 0.6));
  wp_points.push_back(Point_3(0.0, 0.4, 0.8));
  wp_points.push_back(Point_3(0.0, 0.6, 0.0));
  wp_points.push_back(Point_3(0.0, 0.6, 0.2));
  wp_points.push_back(Point_3(0.0, 0.6, 0.4));
  wp_points.push_back(Point_3(0.0, 0.6, 0.6));
  wp_points.push_back(Point_3(0.0, 0.6, 0.8));
  wp_points.push_back(Point_3(0.0, 0.8, 0.0));
  wp_points.push_back(Point_3(0.0, 0.8, 0.2));
  wp_points.push_back(Point_3(0.0, 0.8, 0.4));
  wp_points.push_back(Point_3(0.0, 0.8, 0.6));
  wp_points.push_back(Point_3(0.0, 0.8, 0.8));
  wp_points.push_back(Point_3(0.2, 0.0, 0.0));
  wp_points.push_back(Point_3(0.2, 0.0, 0.2));
  wp_points.push_back(Point_3(0.2, 0.0, 0.4));
  wp_points.push_back(Point_3(0.2, 0.0, 0.6));
  wp_points.push_back(Point_3(0.2, 0.0, 0.8));
  wp_points.push_back(Point_3(0.2, 0.2, 0.0));
  wp_points.push_back(Point_3(0.2, 0.2, 0.2));
  wp_points.push_back(Point_3(0.2, 0.2, 0.4));
  wp_points.push_back(Point_3(0.2, 0.2, 0.6));
  wp_points.push_back(Point_3(0.2, 0.2, 0.8));
  wp_points.push_back(Point_3(0.2, 0.4, 0.0));
  wp_points.push_back(Point_3(0.2, 0.4, 0.2));
  wp_points.push_back(Point_3(0.2, 0.4, 0.4));
  wp_points.push_back(Point_3(0.2, 0.4, 0.6));
  wp_points.push_back(Point_3(0.2, 0.4, 0.8));
  wp_points.push_back(Point_3(0.2, 0.6, 0.0));
  wp_points.push_back(Point_3(0.2, 0.6, 0.2));
  wp_points.push_back(Point_3(0.2, 0.6, 0.4));
  wp_points.push_back(Point_3(0.2, 0.6, 0.6));
  wp_points.push_back(Point_3(0.2, 0.6, 0.8));
  wp_points.push_back(Point_3(0.2, 0.8, 0.0));
  wp_points.push_back(Point_3(0.2, 0.8, 0.2));
  wp_points.push_back(Point_3(0.2, 0.8, 0.4));
  wp_points.push_back(Point_3(0.2, 0.8, 0.6));
  wp_points.push_back(Point_3(0.2, 0.8, 0.8));
  wp_points.push_back(Point_3(0.4, 0.0, 0.0));
  wp_points.push_back(Point_3(0.4, 0.0, 0.2));
  wp_points.push_back(Point_3(0.4, 0.0, 0.4));
  wp_points.push_back(Point_3(0.4, 0.0, 0.6));
  wp_points.push_back(Point_3(0.4, 0.0, 0.8));
  wp_points.push_back(Point_3(0.4, 0.2, 0.0));
  wp_points.push_back(Point_3(0.4, 0.2, 0.2));
  wp_points.push_back(Point_3(0.4, 0.2, 0.4));
  wp_points.push_back(Point_3(0.4, 0.2, 0.6));
  wp_points.push_back(Point_3(0.4, 0.2, 0.8));
  wp_points.push_back(Point_3(0.4, 0.4, 0.0));
  wp_points.push_back(Point_3(0.4, 0.4, 0.2));
  wp_points.push_back(Point_3(0.4, 0.4, 0.4));
  wp_points.push_back(Point_3(0.4, 0.4, 0.6));
  wp_points.push_back(Point_3(0.4, 0.4, 0.8));
  wp_points.push_back(Point_3(0.4, 0.6, 0.0));
  wp_points.push_back(Point_3(0.4, 0.6, 0.2));
  wp_points.push_back(Point_3(0.4, 0.6, 0.4));
  wp_points.push_back(Point_3(0.4, 0.6, 0.6));
  wp_points.push_back(Point_3(0.4, 0.6, 0.8));
  wp_points.push_back(Point_3(0.4, 0.8, 0.0));
  wp_points.push_back(Point_3(0.4, 0.8, 0.2));
  wp_points.push_back(Point_3(0.4, 0.8, 0.4));
  wp_points.push_back(Point_3(0.4, 0.8, 0.6));
  wp_points.push_back(Point_3(0.4, 0.8, 0.8));
  wp_points.push_back(Point_3(0.6, 0.0, 0.0));
  wp_points.push_back(Point_3(0.6, 0.0, 0.2));
  wp_points.push_back(Point_3(0.6, 0.0, 0.4));
  wp_points.push_back(Point_3(0.6, 0.0, 0.6));
  wp_points.push_back(Point_3(0.6, 0.0, 0.8));
  wp_points.push_back(Point_3(0.6, 0.1, 0.0));
  wp_points.push_back(Point_3(0.6, 0.2, 0.0));
  wp_points.push_back(Point_3(0.6, 0.2, 0.2));
  wp_points.push_back(Point_3(0.6, 0.2, 0.4));
  wp_points.push_back(Point_3(0.6, 0.2, 0.6));
  wp_points.push_back(Point_3(0.6, 0.2, 0.8));
  wp_points.push_back(Point_3(0.6, 0.4, 0.0));
  wp_points.push_back(Point_3(0.6, 0.4, 0.2));
  wp_points.push_back(Point_3(0.6, 0.4, 0.4));
  wp_points.push_back(Point_3(0.6, 0.4, 0.6));
  wp_points.push_back(Point_3(0.6, 0.4, 0.8));
  wp_points.push_back(Point_3(0.6, 0.6, 0.0));
  wp_points.push_back(Point_3(0.6, 0.6, 0.2));
  wp_points.push_back(Point_3(0.6, 0.6, 0.4));
  wp_points.push_back(Point_3(0.6, 0.6, 0.6));
  wp_points.push_back(Point_3(0.6, 0.6, 0.8));
  wp_points.push_back(Point_3(0.6, 0.8, 0.0));
  wp_points.push_back(Point_3(0.6, 0.8, 0.2));
  wp_points.push_back(Point_3(0.6, 0.8, 0.4));
  wp_points.push_back(Point_3(0.6, 0.8, 0.6));
  wp_points.push_back(Point_3(0.6, 0.8, 0.8));
  wp_points.push_back(Point_3(0.8, 0.0, 0.0));
  wp_points.push_back(Point_3(0.8, 0.0, 0.2));
  wp_points.push_back(Point_3(0.8, 0.0, 0.4));
  wp_points.push_back(Point_3(0.8, 0.0, 0.6));
  wp_points.push_back(Point_3(0.8, 0.0, 0.8));
  wp_points.push_back(Point_3(0.8, 0.2, 0.0));
  wp_points.push_back(Point_3(0.8, 0.2, 0.2));
  wp_points.push_back(Point_3(0.8, 0.2, 0.4));
  wp_points.push_back(Point_3(0.8, 0.2, 0.6));
  wp_points.push_back(Point_3(0.8, 0.2, 0.8));
  wp_points.push_back(Point_3(0.8, 0.4, 0.0));
  wp_points.push_back(Point_3(0.8, 0.4, 0.2));
  wp_points.push_back(Point_3(0.8, 0.4, 0.4));
  wp_points.push_back(Point_3(0.8, 0.4, 0.6));
  wp_points.push_back(Point_3(0.8, 0.4, 0.8));
  wp_points.push_back(Point_3(0.8, 0.6, 0.0));
  wp_points.push_back(Point_3(0.8, 0.6, 0.2));
  wp_points.push_back(Point_3(0.8, 0.6, 0.4));
  wp_points.push_back(Point_3(0.8, 0.6, 0.6));
  wp_points.push_back(Point_3(0.8, 0.6, 0.8));
  wp_points.push_back(Point_3(0.8, 0.8, 0.0));
  wp_points.push_back(Point_3(0.8, 0.8, 0.2));
  wp_points.push_back(Point_3(0.8, 0.8, 0.4));
  wp_points.push_back(Point_3(0.8, 0.8, 0.6));

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
  // Check it throws an exception when the cuboid is not iso
  BOOST_CHECK_THROW (Weighted_periodic_alpha_complex_3d wp_alpha_complex(wp_points, p_weights, 0., 0., 0., 0.9, 1., 1.),
                     std::invalid_argument);
  BOOST_CHECK_THROW (Weighted_periodic_alpha_complex_3d wp_alpha_complex(wp_points, p_weights, 0., 0., 0., 1., 0.9, 1.),
                     std::invalid_argument);
  BOOST_CHECK_THROW (Weighted_periodic_alpha_complex_3d wp_alpha_complex(wp_points, p_weights, 0., 0., 0., 1., 1., 0.9),
                     std::invalid_argument);

  std::cout << "0 <= point.weight() < 1/64 * domain_size * domain_size exception" << std::endl;
  // Weights must be in range [0, <1/64]
  p_weights[25] = 1.0;
  BOOST_CHECK_THROW (Weighted_periodic_alpha_complex_3d wp_alpha_complex(wp_points, p_weights, 0., 0., 0., 1., 1., 1.),
                     std::invalid_argument);
  // Weights must be in range [0, <1/64]
  p_weights[25] = 0.012;
  p_weights[14] = -0.012;
  BOOST_CHECK_THROW (Weighted_periodic_alpha_complex_3d wp_alpha_complex(wp_points, p_weights, 0., 0., 0., 1., 1., 1.),
                     std::invalid_argument);
  p_weights[14] = 0.005;

  std::cout << "wp_points and p_weights size exception" << std::endl;
  // Weights and points must have the same size
  // + 1
  p_weights.push_back(0.007);
  BOOST_CHECK_THROW (Weighted_periodic_alpha_complex_3d wp_alpha_complex(wp_points, p_weights, 0., 0., 0., 1., 1., 1.),
                     std::invalid_argument);
  // - 1
  p_weights.pop_back();
  p_weights.pop_back();
  BOOST_CHECK_THROW (Weighted_periodic_alpha_complex_3d wp_alpha_complex(wp_points, p_weights, 0., 0., 0., 1., 1., 1.),
                     std::invalid_argument);
}
#endif

BOOST_AUTO_TEST_CASE(Alpha_complex_weighted_periodic) {
  // ------------------------------
  // Fast weighted periodic version
  // ------------------------------
  std::cout << "Fast weighted periodic alpha complex 3d" << std::endl;

  std::vector<Fast_weighted_periodic_alpha_complex_3d::Point_3> wp_points;
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.0, 0.0, 0.2));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.0, 0.0, 0.4));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.0, 0.0, 0.6));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.0, 0.0, 0.8));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.0, 0.2, 0.0));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.0, 0.2, 0.2));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.0, 0.2, 0.4));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.0, 0.2, 0.6));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.0, 0.2, 0.8));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.0, 0.4, 0.0));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.0, 0.4, 0.2));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.0, 0.4, 0.4));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.0, 0.4, 0.6));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.0, 0.4, 0.8));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.0, 0.6, 0.0));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.0, 0.6, 0.2));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.0, 0.6, 0.4));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.0, 0.6, 0.6));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.0, 0.6, 0.8));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.0, 0.8, 0.0));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.0, 0.8, 0.2));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.0, 0.8, 0.4));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.0, 0.8, 0.6));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.0, 0.8, 0.8));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.2, 0.0, 0.0));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.2, 0.0, 0.2));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.2, 0.0, 0.4));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.2, 0.0, 0.6));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.2, 0.0, 0.8));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.2, 0.2, 0.0));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.2, 0.2, 0.2));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.2, 0.2, 0.4));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.2, 0.2, 0.6));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.2, 0.2, 0.8));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.2, 0.4, 0.0));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.2, 0.4, 0.2));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.2, 0.4, 0.4));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.2, 0.4, 0.6));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.2, 0.4, 0.8));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.2, 0.6, 0.0));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.2, 0.6, 0.2));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.2, 0.6, 0.4));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.2, 0.6, 0.6));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.2, 0.6, 0.8));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.0, 0.0, 0.0));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.2, 0.8, 0.0));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.2, 0.8, 0.2));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.2, 0.8, 0.4));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.2, 0.8, 0.6));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.2, 0.8, 0.8));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.4, 0.0, 0.0));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.4, 0.0, 0.2));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.4, 0.0, 0.4));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.4, 0.0, 0.6));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.4, 0.0, 0.8));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.4, 0.2, 0.0));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.4, 0.2, 0.2));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.4, 0.2, 0.4));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.4, 0.2, 0.6));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.4, 0.2, 0.8));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.4, 0.4, 0.0));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.4, 0.4, 0.2));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.4, 0.4, 0.4));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.4, 0.4, 0.6));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.4, 0.4, 0.8));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.4, 0.6, 0.0));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.4, 0.6, 0.2));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.4, 0.6, 0.4));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.4, 0.6, 0.6));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.4, 0.6, 0.8));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.4, 0.8, 0.0));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.4, 0.8, 0.2));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.4, 0.8, 0.4));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.4, 0.8, 0.6));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.4, 0.8, 0.8));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.6, 0.0, 0.0));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.6, 0.0, 0.2));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.6, 0.0, 0.4));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.6, 0.0, 0.6));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.6, 0.0, 0.8));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.6, 0.1, 0.0));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.6, 0.2, 0.0));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.6, 0.2, 0.2));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.6, 0.2, 0.4));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.6, 0.2, 0.6));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.6, 0.2, 0.8));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.6, 0.4, 0.0));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.6, 0.4, 0.2));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.6, 0.4, 0.4));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.6, 0.4, 0.6));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.6, 0.4, 0.8));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.6, 0.6, 0.0));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.6, 0.6, 0.2));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.6, 0.6, 0.4));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.6, 0.6, 0.6));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.6, 0.6, 0.8));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.6, 0.8, 0.0));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.6, 0.8, 0.2));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.6, 0.8, 0.4));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.6, 0.8, 0.6));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.6, 0.8, 0.8));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.8, 0.0, 0.0));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.8, 0.0, 0.2));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.8, 0.0, 0.4));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.8, 0.0, 0.6));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.8, 0.0, 0.8));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.8, 0.2, 0.0));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.8, 0.2, 0.2));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.8, 0.2, 0.4));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.8, 0.2, 0.6));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.8, 0.2, 0.8));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.8, 0.4, 0.0));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.8, 0.4, 0.2));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.8, 0.4, 0.4));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.8, 0.4, 0.6));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.8, 0.4, 0.8));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.8, 0.6, 0.0));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.8, 0.6, 0.2));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.8, 0.6, 0.4));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.8, 0.6, 0.6));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.8, 0.6, 0.8));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.8, 0.8, 0.0));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.8, 0.8, 0.2));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.8, 0.8, 0.4));
  wp_points.push_back(Fast_weighted_periodic_alpha_complex_3d::Point_3(0.8, 0.8, 0.6));

  std::vector<double> p_weights;

  std::random_device rd;
  std::mt19937 mt(rd());
  // Weights must be in range [0, <1/64]
  std::uniform_real_distribution<double> dist(0.01, 0.0156245);

  for (std::size_t i = 0; i < wp_points.size(); ++i) {
    double value = dist(mt);
    p_weights.push_back(value);
  }

  Fast_weighted_periodic_alpha_complex_3d weighted_periodic_alpha_complex(wp_points, p_weights, 0., 0., 0., 1., 1., 1.);

  Gudhi::Simplex_tree<> stree;
  weighted_periodic_alpha_complex.create_complex(stree);

  // -------------------------------
  // Exact weighted periodic version
  // -------------------------------
  std::cout << "Exact weighted periodic alpha complex 3d" << std::endl;

  std::vector<Exact_weighted_periodic_alpha_complex_3d::Point_3> e_wp_points;
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.0, 0.0, 0.4));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.0, 0.0, 0.6));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.0, 0.0, 0.8));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.0, 0.2, 0.0));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.0, 0.2, 0.2));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.0, 0.2, 0.4));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.0, 0.2, 0.6));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.0, 0.2, 0.8));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.0, 0.4, 0.0));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.0, 0.4, 0.2));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.0, 0.4, 0.4));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.0, 0.4, 0.6));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.0, 0.4, 0.8));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.0, 0.6, 0.0));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.0, 0.6, 0.2));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.0, 0.6, 0.4));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.0, 0.0, 0.2));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.0, 0.6, 0.6));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.0, 0.6, 0.8));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.0, 0.8, 0.0));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.0, 0.8, 0.2));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.0, 0.8, 0.4));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.0, 0.8, 0.6));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.0, 0.8, 0.8));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.2, 0.0, 0.0));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.2, 0.0, 0.2));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.2, 0.0, 0.4));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.2, 0.0, 0.6));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.2, 0.0, 0.8));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.2, 0.2, 0.0));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.2, 0.2, 0.2));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.2, 0.2, 0.4));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.2, 0.2, 0.6));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.2, 0.2, 0.8));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.2, 0.4, 0.0));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.2, 0.4, 0.2));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.2, 0.4, 0.4));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.2, 0.4, 0.6));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.2, 0.4, 0.8));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.2, 0.6, 0.0));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.2, 0.6, 0.2));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.2, 0.6, 0.4));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.2, 0.6, 0.6));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.2, 0.6, 0.8));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.0, 0.0, 0.0));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.2, 0.8, 0.0));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.2, 0.8, 0.2));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.2, 0.8, 0.4));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.2, 0.8, 0.6));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.2, 0.8, 0.8));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.4, 0.0, 0.0));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.4, 0.0, 0.2));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.4, 0.0, 0.4));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.4, 0.0, 0.6));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.4, 0.0, 0.8));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.4, 0.2, 0.0));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.4, 0.2, 0.2));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.4, 0.2, 0.4));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.4, 0.2, 0.6));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.4, 0.2, 0.8));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.4, 0.4, 0.0));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.4, 0.4, 0.2));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.4, 0.4, 0.4));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.4, 0.4, 0.6));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.4, 0.4, 0.8));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.4, 0.6, 0.0));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.4, 0.6, 0.2));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.4, 0.6, 0.4));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.4, 0.6, 0.6));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.4, 0.6, 0.8));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.4, 0.8, 0.0));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.4, 0.8, 0.2));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.4, 0.8, 0.4));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.4, 0.8, 0.6));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.4, 0.8, 0.8));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.6, 0.0, 0.0));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.6, 0.0, 0.2));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.6, 0.0, 0.4));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.6, 0.0, 0.6));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.6, 0.0, 0.8));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.6, 0.1, 0.0));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.6, 0.2, 0.0));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.6, 0.2, 0.2));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.6, 0.2, 0.4));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.6, 0.2, 0.6));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.6, 0.2, 0.8));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.6, 0.4, 0.0));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.6, 0.4, 0.2));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.6, 0.4, 0.4));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.6, 0.4, 0.6));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.6, 0.4, 0.8));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.6, 0.6, 0.0));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.6, 0.6, 0.2));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.6, 0.6, 0.4));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.6, 0.6, 0.6));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.6, 0.6, 0.8));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.6, 0.8, 0.0));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.6, 0.8, 0.2));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.6, 0.8, 0.4));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.6, 0.8, 0.6));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.6, 0.8, 0.8));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.8, 0.0, 0.0));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.8, 0.0, 0.2));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.8, 0.0, 0.4));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.8, 0.0, 0.6));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.8, 0.0, 0.8));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.8, 0.2, 0.0));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.8, 0.2, 0.2));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.8, 0.2, 0.4));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.8, 0.2, 0.6));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.8, 0.2, 0.8));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.8, 0.4, 0.0));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.8, 0.4, 0.2));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.8, 0.4, 0.4));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.8, 0.4, 0.6));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.8, 0.4, 0.8));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.8, 0.6, 0.0));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.8, 0.6, 0.2));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.8, 0.6, 0.4));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.8, 0.6, 0.6));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.8, 0.6, 0.8));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.8, 0.8, 0.0));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.8, 0.8, 0.2));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.8, 0.8, 0.4));
  e_wp_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(0.8, 0.8, 0.6));

  Exact_weighted_periodic_alpha_complex_3d e_weighted_periodic_alpha_complex(e_wp_points, p_weights, 0., 0., 0., 1., 1., 1.);

  Gudhi::Simplex_tree<> exact_stree;
  e_weighted_periodic_alpha_complex.create_complex(exact_stree);

  // ---------------------
  // Compare both versions
  // ---------------------
  std::cout << "Exact periodic alpha complex 3d is of dimension " << exact_stree.dimension()
            << " - Non exact is " << stree.dimension() << std::endl;
  BOOST_CHECK(exact_stree.dimension() == stree.dimension());
  std::cout << "Exact periodic alpha complex 3d num_simplices " << exact_stree.num_simplices()
            << " - Non exact is " << stree.num_simplices() << std::endl;
  // TODO(VR): BOOST_CHECK(exact_stree.num_simplices() == stree.num_simplices());
  std::cout << "Exact periodic alpha complex 3d num_vertices " << exact_stree.num_vertices()
            << " - Non exact is " << stree.num_vertices() << std::endl;
  BOOST_CHECK(exact_stree.num_vertices() == stree.num_vertices());

  auto sh = stree.filtration_simplex_range().begin();
  while(sh != stree.filtration_simplex_range().end()) {
    std::vector<int> simplex;
    std::vector<int> exact_simplex;
    std::cout << "Non-exact ( ";
    for (auto vertex : stree.simplex_vertex_range(*sh)) {
      simplex.push_back(vertex);
      std::cout << vertex << " ";
    }
    std::cout << ") -> " << "[" << stree.filtration(*sh) << "] ";
    std::cout << std::endl;

    // Find it in the exact structure
    auto sh_exact = exact_stree.find(simplex);
    // TODO(VR): BOOST_CHECK(sh_exact != exact_stree.null_simplex());

    // Exact and non-exact version is not exactly the same due to float comparison
    // TODO(VR): GUDHI_TEST_FLOAT_EQUALITY_CHECK(exact_stree.filtration(sh_exact), stree.filtration(*sh));
    ++sh;
  }

}
