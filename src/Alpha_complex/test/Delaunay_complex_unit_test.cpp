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
#define BOOST_TEST_MODULE "delaunay_complex"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <CGAL/Epick_d.h>
#include <CGAL/Epeck_d.h>

#include <vector>
#include <limits>  // NaN
#include <cmath>

#include <gudhi/Alpha_complex.h>
// to construct a simplex_tree from Delaunay_triangulation
#include <gudhi/graph_simplicial_complex.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Unitary_tests_utils.h>
#include <gudhi/random_point_generators.h>

// Use dynamic_dimension_tag for the user to be able to set dimension
typedef CGAL::Epeck_d< CGAL::Dynamic_dimension_tag > Exact_kernel_d;
// Use static dimension_tag for the user not to be able to set dimension
typedef CGAL::Epeck_d< CGAL::Dimension_tag<5> > Exact_kernel_s;
// Use dynamic_dimension_tag for the user to be able to set dimension
typedef CGAL::Epick_d< CGAL::Dynamic_dimension_tag > Inexact_kernel_d;
// Use static dimension_tag for the user not to be able to set dimension
typedef CGAL::Epick_d< CGAL::Dimension_tag<5> > Inexact_kernel_s;
// The triangulation uses the default instantiation of the TriangulationDataStructure template parameter

typedef boost::mpl::list<Exact_kernel_d, Exact_kernel_s, Inexact_kernel_d, Inexact_kernel_s> list_of_kernel_variants;

using Simplex_tree = Gudhi::Simplex_tree<>;
using Simplex_handle = Simplex_tree::Simplex_handle;

BOOST_AUTO_TEST_CASE_TEMPLATE(Alpha_complex_from_OFF_file, TestedKernel, list_of_kernel_variants) {
  std::cout << "*****************************************************************************************************";
  using Point = typename TestedKernel::Point_d;
  std::vector<Point> points;
  // 50 points on a 4-sphere
  points = Gudhi::generate_points_on_sphere_d<TestedKernel>(10, 5, 1.);

  Gudhi::alpha_complex::Alpha_complex<TestedKernel> alpha_complex(points);

  // Alpha complex
  Simplex_tree stree_from_alpha_complex;
  BOOST_CHECK(alpha_complex.create_complex(stree_from_alpha_complex));

  // Delaunay complex
  Simplex_tree stree_from_delaunay_complex;
  BOOST_CHECK(alpha_complex.create_complex(stree_from_delaunay_complex, 0., false, true));

  // Check all the simplices from alpha complex are in the Delaunay complex
  for (auto f_simplex : stree_from_alpha_complex.complex_simplex_range()) {
    Simplex_handle sh = stree_from_delaunay_complex.find(stree_from_alpha_complex.simplex_vertex_range(f_simplex));
    BOOST_CHECK(std::isnan(stree_from_delaunay_complex.filtration(sh)));
    BOOST_CHECK(sh != stree_from_delaunay_complex.null_simplex());
  }
}
