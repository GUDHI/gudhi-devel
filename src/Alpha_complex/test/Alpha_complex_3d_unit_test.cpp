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

using Fast_alpha_complex_3d =
    Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::FAST, false, false>;
using Safe_alpha_complex_3d =
    Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::SAFE, false, false>;
using Exact_alpha_complex_3d =
    Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::EXACT, false, false>;

using Fast_weighted_alpha_complex_3d =
    Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::FAST, true, false>;
using Safe_weighted_alpha_complex_3d =
    Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::SAFE, true, false>;
using Exact_weighted_alpha_complex_3d =
    Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::EXACT, true, false>;

using Fast_periodic_alpha_complex_3d =
    Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::FAST, false, true>;
using Safe_periodic_alpha_complex_3d =
    Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::SAFE, false, true>;
using Exact_periodic_alpha_complex_3d =
    Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::EXACT, false, true>;

using Fast_weighted_periodic_alpha_complex_3d =
    Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::FAST, true, true>;
using Safe_weighted_periodic_alpha_complex_3d =
    Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::SAFE, true, true>;
using Exact_weighted_periodic_alpha_complex_3d =
    Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::EXACT, true, true>;

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
  using Weighted_point_3 = typename Weighted_alpha_complex_3d::Triangulation_3::Weighted_point;

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

#ifdef GUDHI_DEBUG
typedef boost::mpl::list<Fast_periodic_alpha_complex_3d, Safe_periodic_alpha_complex_3d,
                         Exact_periodic_alpha_complex_3d>
    periodic_variants_type_list;

BOOST_AUTO_TEST_CASE_TEMPLATE(Alpha_complex_periodic_throw, Periodic_alpha_complex_3d, periodic_variants_type_list) {
  std::cout << "Periodic alpha complex 3d exception throw" << std::endl;
  using Point_3 = typename Periodic_alpha_complex_3d::Point_3;
  std::vector<Point_3> p_points;

  // Not important, this is not what we want to check
  p_points.push_back(Point_3(0.0, 0.0, 0.0));

  std::cout << "Check exception throw in debug mode" << std::endl;
  // Check it throws an exception when the cuboid is not iso
  BOOST_CHECK_THROW(Periodic_alpha_complex_3d periodic_alpha_complex(p_points, 0., 0., 0., 0.9, 1., 1.),
                    std::invalid_argument);
  BOOST_CHECK_THROW(Periodic_alpha_complex_3d periodic_alpha_complex(p_points, 0., 0., 0., 1., 0.9, 1.),
                    std::invalid_argument);
  BOOST_CHECK_THROW(Periodic_alpha_complex_3d periodic_alpha_complex(p_points, 0., 0., 0., 1., 1., 0.9),
                    std::invalid_argument);
  BOOST_CHECK_THROW(Periodic_alpha_complex_3d periodic_alpha_complex(p_points, 0., 0., 0., 1.1, 1., 1.),
                    std::invalid_argument);
  BOOST_CHECK_THROW(Periodic_alpha_complex_3d periodic_alpha_complex(p_points, 0., 0., 0., 1., 1.1, 1.),
                    std::invalid_argument);
  BOOST_CHECK_THROW(Periodic_alpha_complex_3d periodic_alpha_complex(p_points, 0., 0., 0., 1., 1., 1.1),
                    std::invalid_argument);
}
#endif

BOOST_AUTO_TEST_CASE(Alpha_complex_periodic) {
  // ---------------------
  // Fast periodic version
  // ---------------------
  std::cout << "Fast periodic alpha complex 3d" << std::endl;

  using Creator = CGAL::Creator_uniform_3<double, Fast_periodic_alpha_complex_3d::Point_3>;
  CGAL::Random random(7);
  CGAL::Random_points_in_cube_3<Fast_periodic_alpha_complex_3d::Point_3, Creator> in_cube(1, random);
  std::vector<Fast_periodic_alpha_complex_3d::Point_3> p_points;

  for (int i = 0; i < 50; i++) {
    Fast_periodic_alpha_complex_3d::Point_3 p = *in_cube++;
    p_points.push_back(p);
  }

  Fast_periodic_alpha_complex_3d periodic_alpha_complex(p_points, -1., -1., -1., 1., 1., 1.);

  Gudhi::Simplex_tree<> stree;
  periodic_alpha_complex.create_complex(stree);

  // ----------------------
  // Exact periodic version
  // ----------------------
  std::cout << "Exact periodic alpha complex 3d" << std::endl;

  std::vector<Exact_periodic_alpha_complex_3d::Point_3> e_p_points;

  for (auto p : p_points) {
    e_p_points.push_back(Exact_periodic_alpha_complex_3d::Point_3(p[0], p[1], p[2]));
  }

  Exact_periodic_alpha_complex_3d exact_alpha_complex(e_p_points, -1., -1., -1., 1., 1., 1.);

  Gudhi::Simplex_tree<> exact_stree;
  exact_alpha_complex.create_complex(exact_stree);

  // ---------------------
  // Compare both versions
  // ---------------------
  std::cout << "Exact periodic alpha complex 3d is of dimension " << exact_stree.dimension() << " - Non exact is "
            << stree.dimension() << std::endl;
  BOOST_CHECK(exact_stree.dimension() == stree.dimension());
  std::cout << "Exact periodic alpha complex 3d num_simplices " << exact_stree.num_simplices() << " - Non exact is "
            << stree.num_simplices() << std::endl;
  BOOST_CHECK(exact_stree.num_simplices() == stree.num_simplices());
  std::cout << "Exact periodic alpha complex 3d num_vertices " << exact_stree.num_vertices() << " - Non exact is "
            << stree.num_vertices() << std::endl;
  BOOST_CHECK(exact_stree.num_vertices() == stree.num_vertices());

  // We cannot compare as objects from dispatcher on the alpha shape is not deterministic.
  // cf. https://github.com/CGAL/cgal/issues/3346
  auto sh = stree.filtration_simplex_range().begin();
  auto sh_exact = exact_stree.filtration_simplex_range().begin();

  while (sh != stree.filtration_simplex_range().end() || sh_exact != exact_stree.filtration_simplex_range().end()) {
    GUDHI_TEST_FLOAT_EQUALITY_CHECK(stree.filtration(*sh), exact_stree.filtration(*sh_exact), 1e-14);

    std::vector<int> vh(stree.simplex_vertex_range(*sh).begin(), stree.simplex_vertex_range(*sh).end());
    std::vector<int> exact_vh(exact_stree.simplex_vertex_range(*sh_exact).begin(),
                              exact_stree.simplex_vertex_range(*sh_exact).end());

    BOOST_CHECK(vh.size() == exact_vh.size());
    ++sh;
    ++sh_exact;
  }

  BOOST_CHECK(sh == stree.filtration_simplex_range().end());
  BOOST_CHECK(sh_exact == exact_stree.filtration_simplex_range().end());

  // ----------------------
  // Safe periodic version
  // ----------------------
  std::cout << "Safe periodic alpha complex 3d" << std::endl;

  std::vector<Safe_periodic_alpha_complex_3d::Point_3> s_p_points;

  for (auto p : p_points) {
    s_p_points.push_back(Safe_periodic_alpha_complex_3d::Point_3(p[0], p[1], p[2]));
  }

  Safe_periodic_alpha_complex_3d safe_alpha_complex(s_p_points, -1., -1., -1., 1., 1., 1.);

  Gudhi::Simplex_tree<> safe_stree;
  safe_alpha_complex.create_complex(safe_stree);

  // ---------------------
  // Compare both versions
  // ---------------------
  // We cannot compare as objects from dispatcher on the alpha shape is not deterministic.
  // cf. https://github.com/CGAL/cgal/issues/3346
  sh = stree.filtration_simplex_range().begin();
  auto sh_safe = safe_stree.filtration_simplex_range().begin();

  while (sh != stree.filtration_simplex_range().end() || sh_safe != safe_stree.filtration_simplex_range().end()) {
    GUDHI_TEST_FLOAT_EQUALITY_CHECK(stree.filtration(*sh), safe_stree.filtration(*sh_safe), 1e-14);

    std::vector<int> vh(stree.simplex_vertex_range(*sh).begin(), stree.simplex_vertex_range(*sh).end());
    std::vector<int> safe_vh(safe_stree.simplex_vertex_range(*sh_safe).begin(),
                             safe_stree.simplex_vertex_range(*sh_safe).end());

    BOOST_CHECK(vh.size() == safe_vh.size());
    ++sh;
    ++sh_safe;
  }

  BOOST_CHECK(sh == stree.filtration_simplex_range().end());
  BOOST_CHECK(sh_safe == safe_stree.filtration_simplex_range().end());
}

typedef boost::mpl::list<Fast_weighted_periodic_alpha_complex_3d, Exact_weighted_periodic_alpha_complex_3d,
                         Safe_weighted_periodic_alpha_complex_3d>
    wp_variants_type_list;

#ifdef GUDHI_DEBUG
BOOST_AUTO_TEST_CASE_TEMPLATE(Alpha_complex_weighted_periodic_throw, Weighted_periodic_alpha_complex_3d,
                              wp_variants_type_list) {
  std::cout << "Weighted periodic alpha complex 3d exception throw" << std::endl;

  using Creator = CGAL::Creator_uniform_3<double, Weighted_periodic_alpha_complex_3d::Point_3>;
  CGAL::Random random(7);
  CGAL::Random_points_in_cube_3<Weighted_periodic_alpha_complex_3d::Point_3, Creator> in_cube(1, random);
  std::vector<Weighted_periodic_alpha_complex_3d::Point_3> wp_points;

  for (int i = 0; i < 50; i++) {
    Weighted_periodic_alpha_complex_3d::Point_3 p = *in_cube++;
    wp_points.push_back(p);
  }
  std::vector<double> p_weights;
  // Weights must be in range ]0, 1/64 = 0.015625[
  for (std::size_t i = 0; i < wp_points.size(); ++i) {
    p_weights.push_back(random.get_double(0., 0.01));
  }

  std::cout << "Cuboid is not iso exception" << std::endl;
  // Check it throws an exception when the cuboid is not iso
  BOOST_CHECK_THROW(
      Weighted_periodic_alpha_complex_3d wp_alpha_complex(wp_points, p_weights, -1., -1., -1., 0.9, 1., 1.),
      std::invalid_argument);
  BOOST_CHECK_THROW(
      Weighted_periodic_alpha_complex_3d wp_alpha_complex(wp_points, p_weights, -1., -1., -1., 1., 0.9, 1.),
      std::invalid_argument);
  BOOST_CHECK_THROW(
      Weighted_periodic_alpha_complex_3d wp_alpha_complex(wp_points, p_weights, -1., -1., -1., 1., 1., 0.9),
      std::invalid_argument);
  BOOST_CHECK_THROW(
      Weighted_periodic_alpha_complex_3d wp_alpha_complex(wp_points, p_weights, -1., -1., -1., 1.1, 1., 1.),
      std::invalid_argument);
  BOOST_CHECK_THROW(
      Weighted_periodic_alpha_complex_3d wp_alpha_complex(wp_points, p_weights, -1., -1., -1., 1., 1.1, 1.),
      std::invalid_argument);
  BOOST_CHECK_THROW(
      Weighted_periodic_alpha_complex_3d wp_alpha_complex(wp_points, p_weights, -1., -1., -1., 1., 1., 1.1),
      std::invalid_argument);

  std::cout << "0 <= point.weight() < 1/64 * domain_size * domain_size exception" << std::endl;
  // Weights must be in range ]0, 1/64 = 0.015625[
  double temp = p_weights[25];
  p_weights[25] = 1.0;
  BOOST_CHECK_THROW(Weighted_periodic_alpha_complex_3d wp_alpha_complex(wp_points, p_weights, 0., 0., 0., 1., 1., 1.),
                    std::invalid_argument);
  // Weights must be in range ]0, 1/64 = 0.015625[
  p_weights[25] = temp;
  temp = p_weights[14];
  p_weights[14] = -1e-10;
  BOOST_CHECK_THROW(Weighted_periodic_alpha_complex_3d wp_alpha_complex(wp_points, p_weights, 0., 0., 0., 1., 1., 1.),
                    std::invalid_argument);
  p_weights[14] = temp;

  std::cout << "wp_points and p_weights size exception" << std::endl;
  // Weights and points must have the same size
  // + 1
  p_weights.push_back(1e-10);
  BOOST_CHECK_THROW(Weighted_periodic_alpha_complex_3d wp_alpha_complex(wp_points, p_weights, 0., 0., 0., 1., 1., 1.),
                    std::invalid_argument);
  // - 1
  p_weights.pop_back();
  p_weights.pop_back();
  BOOST_CHECK_THROW(Weighted_periodic_alpha_complex_3d wp_alpha_complex(wp_points, p_weights, 0., 0., 0., 1., 1., 1.),
                    std::invalid_argument);
}
#endif

BOOST_AUTO_TEST_CASE(Alpha_complex_weighted_periodic) {
  // ---------------------
  // Fast weighted periodic version
  // ---------------------
  std::cout << "Fast weighted periodic alpha complex 3d" << std::endl;

  using Creator = CGAL::Creator_uniform_3<double, Fast_weighted_periodic_alpha_complex_3d::Point_3>;
  CGAL::Random random(7);
  CGAL::Random_points_in_cube_3<Fast_weighted_periodic_alpha_complex_3d::Point_3, Creator> in_cube(1, random);
  std::vector<Fast_weighted_periodic_alpha_complex_3d::Point_3> p_points;

  for (int i = 0; i < 50; i++) {
    Fast_weighted_periodic_alpha_complex_3d::Point_3 p = *in_cube++;
    p_points.push_back(p);
  }
  std::vector<double> p_weights;
  // Weights must be in range ]0, 1/64 = 0.015625[
  for (std::size_t i = 0; i < p_points.size(); ++i) {
    p_weights.push_back(random.get_double(0., 0.01));
  }

  Fast_weighted_periodic_alpha_complex_3d periodic_alpha_complex(p_points, p_weights, -1., -1., -1., 1., 1., 1.);

  Gudhi::Simplex_tree<> stree;
  periodic_alpha_complex.create_complex(stree);

  // ----------------------
  // Exact weighted periodic version
  // ----------------------
  std::cout << "Exact weighted periodic alpha complex 3d" << std::endl;

  std::vector<Exact_weighted_periodic_alpha_complex_3d::Point_3> e_p_points;

  for (auto p : p_points) {
    e_p_points.push_back(Exact_weighted_periodic_alpha_complex_3d::Point_3(p[0], p[1], p[2]));
  }

  Exact_weighted_periodic_alpha_complex_3d exact_alpha_complex(e_p_points, p_weights, -1., -1., -1., 1., 1., 1.);

  Gudhi::Simplex_tree<> exact_stree;
  exact_alpha_complex.create_complex(exact_stree);

  // ---------------------
  // Compare both versions
  // ---------------------
  std::cout << "Exact weighted periodic alpha complex 3d is of dimension " << exact_stree.dimension()
            << " - Non exact is " << stree.dimension() << std::endl;
  BOOST_CHECK(exact_stree.dimension() == stree.dimension());
  std::cout << "Exact weighted periodic alpha complex 3d num_simplices " << exact_stree.num_simplices()
            << " - Non exact is " << stree.num_simplices() << std::endl;
  BOOST_CHECK(exact_stree.num_simplices() == stree.num_simplices());
  std::cout << "Exact weighted periodic alpha complex 3d num_vertices " << exact_stree.num_vertices()
            << " - Non exact is " << stree.num_vertices() << std::endl;
  BOOST_CHECK(exact_stree.num_vertices() == stree.num_vertices());

  // We cannot compare as objects from dispatcher on the alpha shape is not deterministic.
  // cf. https://github.com/CGAL/cgal/issues/3346
  auto sh = stree.filtration_simplex_range().begin();
  auto sh_exact = exact_stree.filtration_simplex_range().begin();

  while (sh != stree.filtration_simplex_range().end() || sh_exact != exact_stree.filtration_simplex_range().end()) {
    GUDHI_TEST_FLOAT_EQUALITY_CHECK(stree.filtration(*sh), exact_stree.filtration(*sh_exact), 1e-14);

    std::vector<int> vh(stree.simplex_vertex_range(*sh).begin(), stree.simplex_vertex_range(*sh).end());
    std::vector<int> exact_vh(exact_stree.simplex_vertex_range(*sh_exact).begin(),
                              exact_stree.simplex_vertex_range(*sh_exact).end());

    BOOST_CHECK(vh.size() == exact_vh.size());
    ++sh;
    ++sh_exact;
  }

  BOOST_CHECK(sh == stree.filtration_simplex_range().end());
  BOOST_CHECK(sh_exact == exact_stree.filtration_simplex_range().end());

  // ----------------------
  // Safe weighted periodic version
  // ----------------------
  std::cout << "Safe weighted periodic alpha complex 3d" << std::endl;

  std::vector<Safe_weighted_periodic_alpha_complex_3d::Point_3> s_p_points;

  for (auto p : p_points) {
    s_p_points.push_back(Safe_weighted_periodic_alpha_complex_3d::Point_3(p[0], p[1], p[2]));
  }

  Safe_weighted_periodic_alpha_complex_3d safe_alpha_complex(s_p_points, p_weights, -1., -1., -1., 1., 1., 1.);

  Gudhi::Simplex_tree<> safe_stree;
  safe_alpha_complex.create_complex(safe_stree);

  // ---------------------
  // Compare both versions
  // ---------------------
  // We cannot compare as objects from dispatcher on the alpha shape is not deterministic.
  // cf. https://github.com/CGAL/cgal/issues/3346
  sh = stree.filtration_simplex_range().begin();
  auto sh_safe = safe_stree.filtration_simplex_range().begin();

  while (sh != stree.filtration_simplex_range().end() || sh_safe != safe_stree.filtration_simplex_range().end()) {
    GUDHI_TEST_FLOAT_EQUALITY_CHECK(stree.filtration(*sh), safe_stree.filtration(*sh_safe), 1e-14);

    std::vector<int> vh(stree.simplex_vertex_range(*sh).begin(), stree.simplex_vertex_range(*sh).end());
    std::vector<int> safe_vh(safe_stree.simplex_vertex_range(*sh_safe).begin(),
                             safe_stree.simplex_vertex_range(*sh_safe).end());

    BOOST_CHECK(vh.size() == safe_vh.size());
    ++sh;
    ++sh_safe;
  }

  BOOST_CHECK(sh == stree.filtration_simplex_range().end());
  BOOST_CHECK(sh_safe == safe_stree.filtration_simplex_range().end());
}
