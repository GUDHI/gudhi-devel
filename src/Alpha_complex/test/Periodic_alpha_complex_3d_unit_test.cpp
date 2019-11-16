/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2015 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
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

using Fast_periodic_alpha_complex_3d =
    Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::FAST, false, true>;
using Safe_periodic_alpha_complex_3d =
    Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::SAFE, false, true>;
using Exact_periodic_alpha_complex_3d =
    Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::EXACT, false, true>;

#ifdef GUDHI_DEBUG
typedef boost::mpl::list<Fast_periodic_alpha_complex_3d, Safe_periodic_alpha_complex_3d,
                         Exact_periodic_alpha_complex_3d>
    periodic_variants_type_list;

BOOST_AUTO_TEST_CASE_TEMPLATE(Alpha_complex_periodic_throw, Periodic_alpha_complex_3d, periodic_variants_type_list) {
  std::cout << "Periodic alpha complex 3d exception throw" << std::endl;
  using Bare_point_3 = typename Periodic_alpha_complex_3d::Bare_point_3;
  std::vector<Bare_point_3> p_points;

  // Not important, this is not what we want to check
  p_points.push_back(Bare_point_3(0.0, 0.0, 0.0));

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

  using Creator = CGAL::Creator_uniform_3<double, Fast_periodic_alpha_complex_3d::Bare_point_3>;
  CGAL::Random random(7);
  CGAL::Random_points_in_cube_3<Fast_periodic_alpha_complex_3d::Bare_point_3, Creator> in_cube(1, random);
  std::vector<Fast_periodic_alpha_complex_3d::Bare_point_3> p_points;

  for (int i = 0; i < 50; i++) {
    Fast_periodic_alpha_complex_3d::Bare_point_3 p = *in_cube++;
    p_points.push_back(p);
  }

  Fast_periodic_alpha_complex_3d periodic_alpha_complex(p_points, -1., -1., -1., 1., 1., 1.);

  Gudhi::Simplex_tree<> stree;
  periodic_alpha_complex.create_complex(stree);

  for (std::size_t index = 0; index < p_points.size(); index++) {
    bool found = false;
    Fast_periodic_alpha_complex_3d::Bare_point_3 ap = periodic_alpha_complex.get_point(index);
    for (auto point : p_points) {
      if ((point.x() == ap.x()) && (point.y() == ap.y()) && (point.z() == ap.z())) {
        found = true;
        break;
      }
    }
    // Check all points from alpha complex are found in the input point cloud
    BOOST_CHECK(found);
  }
  // Exception if we go out of range
  BOOST_CHECK_THROW(periodic_alpha_complex.get_point(p_points.size()), std::out_of_range);

  // ----------------------
  // Exact periodic version
  // ----------------------
  std::cout << "Exact periodic alpha complex 3d" << std::endl;

  std::vector<Exact_periodic_alpha_complex_3d::Bare_point_3> e_p_points;

  for (auto p : p_points) {
    e_p_points.push_back(Exact_periodic_alpha_complex_3d::Bare_point_3(p[0], p[1], p[2]));
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

  std::vector<Safe_periodic_alpha_complex_3d::Bare_point_3> s_p_points;

  for (auto p : p_points) {
    s_p_points.push_back(Safe_periodic_alpha_complex_3d::Bare_point_3(p[0], p[1], p[2]));
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
