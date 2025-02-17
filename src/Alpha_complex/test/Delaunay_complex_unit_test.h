/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2020 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <vector>
#include <limits>  // NaN
#include <cmath>

#include <gudhi/Alpha_complex.h>
// to construct a simplex_tree from Delaunay_triangulation
#include <gudhi/graph_simplicial_complex.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Unitary_tests_utils.h>
#include <gudhi/random_point_generators.h>

using Simplex_tree = Gudhi::Simplex_tree<>;
using Simplex_handle = Simplex_tree::Simplex_handle;

template<class CGAL_kernel>
void compare_delaunay_complex_simplices() {
  std::clog << "***************************************************************************************************\n";
  using Point = typename CGAL_kernel::Point_d;
  std::vector<Point> points;
  // 50 points on a 4-sphere
  points = Gudhi::generate_points_on_sphere_d<CGAL_kernel>(10, 5, 1.);

  Gudhi::alpha_complex::Alpha_complex<CGAL_kernel> alpha_complex(points);

  // Alpha complex
  Simplex_tree stree_from_alpha_complex;
  BOOST_CHECK(alpha_complex.create_complex(stree_from_alpha_complex));

  // Delaunay complex
  Simplex_tree stree_from_delaunay_complex;
  BOOST_CHECK(alpha_complex.create_complex(stree_from_delaunay_complex, 0., false, true));

  std::clog << "Check all the simplices from alpha complex are in the Delaunay complex\n";
  for (auto f_simplex : stree_from_alpha_complex.complex_simplex_range()) {
    Simplex_handle sh = stree_from_delaunay_complex.find(stree_from_alpha_complex.simplex_vertex_range(f_simplex));
    BOOST_CHECK(std::isnan(stree_from_delaunay_complex.filtration(sh)));
    BOOST_CHECK(sh != stree_from_delaunay_complex.null_simplex());
  }

  std::clog << "Alpha complex with square root filtration values\n";
  Simplex_tree stree_from_alpha_sqrt;
  // set Output_squared_values to false
  BOOST_CHECK(alpha_complex.template create_complex<false>(stree_from_alpha_sqrt));

  std::clog << "Check simplices from alpha complex filtration values when Output_squared_values is true\n";
  // Check all the simplices from alpha complex are in the Delaunay complex
  for (auto f_simplex : stree_from_alpha_complex.complex_simplex_range()) {
    std::clog << "[";
    for (auto vertex : stree_from_alpha_complex.simplex_vertex_range(f_simplex)) std::clog << vertex << " ";
    std::clog << "] - ";

    Simplex_handle sh = stree_from_alpha_sqrt.find(stree_from_alpha_complex.simplex_vertex_range(f_simplex));
    BOOST_CHECK(sh != stree_from_alpha_sqrt.null_simplex());
    std::clog << "alpha_sqrt = " << stree_from_alpha_sqrt.filtration(sh);
    std::clog << " vs. sqrt(alpha) = " << std::sqrt(stree_from_alpha_complex.filtration(f_simplex)) << "\n";
    GUDHI_TEST_FLOAT_EQUALITY_CHECK(stree_from_alpha_sqrt.filtration(sh),
                                    std::sqrt(stree_from_alpha_complex.filtration(f_simplex)));
  }

}
