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

using Fast_weighted_periodic_alpha_complex_3d =
    Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::FAST, true, true>;
using Safe_weighted_periodic_alpha_complex_3d =
    Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::SAFE, true, true>;
using Exact_weighted_periodic_alpha_complex_3d =
    Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::EXACT, true, true>;


typedef boost::mpl::list<Fast_weighted_periodic_alpha_complex_3d, Exact_weighted_periodic_alpha_complex_3d,
                         Safe_weighted_periodic_alpha_complex_3d>
    wp_variants_type_list;

#ifdef GUDHI_DEBUG
BOOST_AUTO_TEST_CASE_TEMPLATE(Alpha_complex_weighted_periodic_throw, Weighted_periodic_alpha_complex_3d,
                              wp_variants_type_list) {
  std::cout << "Weighted periodic alpha complex 3d exception throw" << std::endl;

  using Creator = CGAL::Creator_uniform_3<double, typename Weighted_periodic_alpha_complex_3d::Point_3>;
  CGAL::Random random(7);
  CGAL::Random_points_in_cube_3<typename Weighted_periodic_alpha_complex_3d::Point_3, Creator> in_cube(1, random);
  std::vector<typename Weighted_periodic_alpha_complex_3d::Point_3> wp_points;

  for (int i = 0; i < 50; i++) {
    typename Weighted_periodic_alpha_complex_3d::Point_3 p = *in_cube++;
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
