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
#define BOOST_TEST_MODULE "alpha_complex_dim3"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <CGAL/Epick_d.h>
#include <CGAL/Epeck_d.h>

#include <stdexcept> // std::out_of_range
#include <string>
#include <vector>

#include <gudhi/Alpha_complex.h>
#include <gudhi/Simplex_tree.h>

// Use dynamic_dimension_tag for the user to be able to set dimension
typedef CGAL::Epeck_d< CGAL::Dynamic_dimension_tag > Exact_kernel_d;
// Use static dimension_tag for the user not to be able to set dimension
typedef CGAL::Epeck_d< CGAL::Dimension_tag<3> > Exact_kernel_s;
// Use dynamic_dimension_tag for the user to be able to set dimension
typedef CGAL::Epick_d< CGAL::Dynamic_dimension_tag > Inexact_kernel_d;
// Use static dimension_tag for the user not to be able to set dimension
typedef CGAL::Epick_d< CGAL::Dimension_tag<3> > Inexact_kernel_s;
// The triangulation uses the default instantiation of the TriangulationDataStructure template parameter

typedef boost::mpl::list<Exact_kernel_d, Exact_kernel_s, Inexact_kernel_d, Inexact_kernel_s> list_of_kernel_variants;

BOOST_AUTO_TEST_CASE_TEMPLATE(Alpha_complex_from_OFF_file, TestedKernel, list_of_kernel_variants) {
  // ----------------------------------------------------------------------------
  //
  // Init of an alpha-complex from a OFF file
  //
  // ----------------------------------------------------------------------------
  std::string off_file_name("alphacomplexdoc.off");
  double max_alpha_square_value = 60.0;
  std::clog << "========== OFF FILE NAME = " << off_file_name << " - alpha²=" <<
      max_alpha_square_value << "==========" << std::endl;

  Gudhi::alpha_complex::Alpha_complex<TestedKernel> alpha_complex_from_file(off_file_name);

  Gudhi::Simplex_tree<> simplex_tree_60;
  BOOST_CHECK(alpha_complex_from_file.create_complex(simplex_tree_60, max_alpha_square_value));

  std::clog << "alpha_complex_from_file.num_vertices()=" << alpha_complex_from_file.num_vertices() << std::endl;
  BOOST_CHECK(alpha_complex_from_file.num_vertices() == 7);

  std::clog << "simplex_tree_60.dimension()=" << simplex_tree_60.dimension() << std::endl;
  BOOST_CHECK(simplex_tree_60.dimension() == 2);

  std::clog << "simplex_tree_60.num_vertices()=" << simplex_tree_60.num_vertices() << std::endl;
  BOOST_CHECK(simplex_tree_60.num_vertices() == 7);

  std::clog << "simplex_tree_60.num_simplices()=" << simplex_tree_60.num_simplices() << std::endl;
  BOOST_CHECK(simplex_tree_60.num_simplices() == 25);

  max_alpha_square_value = 59.0;
  std::clog << "========== OFF FILE NAME = " << off_file_name << " - alpha²=" <<
      max_alpha_square_value << "==========" << std::endl;

  Gudhi::Simplex_tree<> simplex_tree_59;
  BOOST_CHECK(alpha_complex_from_file.create_complex(simplex_tree_59, max_alpha_square_value));
  
  std::clog << "alpha_complex_from_file.num_vertices()=" << alpha_complex_from_file.num_vertices() << std::endl;
  BOOST_CHECK(alpha_complex_from_file.num_vertices() == 7);

  std::clog << "simplex_tree_59.dimension()=" << simplex_tree_59.dimension() << std::endl;
  BOOST_CHECK(simplex_tree_59.dimension() == 2);

  std::clog << "simplex_tree_59.num_vertices()=" << simplex_tree_59.num_vertices() << std::endl;
  BOOST_CHECK(simplex_tree_59.num_vertices() == 7);

  std::clog << "simplex_tree_59.num_simplices()=" << simplex_tree_59.num_simplices() << std::endl;
  BOOST_CHECK(simplex_tree_59.num_simplices() == 23);
}


BOOST_AUTO_TEST_CASE_TEMPLATE(Alpha_complex_from_empty_points, TestedKernel, list_of_kernel_variants) {
  std::clog << "========== Alpha_complex_from_empty_points ==========" << std::endl;

  // ----------------------------------------------------------------------------
  // Init of an empty list of points
  // ----------------------------------------------------------------------------
  std::vector<typename TestedKernel::Point_d> points;

  // ----------------------------------------------------------------------------
  // Init of an alpha complex from the list of points
  // ----------------------------------------------------------------------------
  Gudhi::alpha_complex::Alpha_complex<TestedKernel> alpha_complex_from_points(points);

  std::clog << "alpha_complex_from_points.num_vertices()=" << alpha_complex_from_points.num_vertices() << std::endl;
  BOOST_CHECK(alpha_complex_from_points.num_vertices() == points.size());

  // Test to the limit
  BOOST_CHECK_THROW (alpha_complex_from_points.get_point(0), std::out_of_range);

  Gudhi::Simplex_tree<> simplex_tree;
  BOOST_CHECK(!alpha_complex_from_points.create_complex(simplex_tree));
  
  std::clog << "simplex_tree.num_simplices()=" << simplex_tree.num_simplices() << std::endl;
  BOOST_CHECK(simplex_tree.num_simplices() == 0);

  std::clog << "simplex_tree.dimension()=" << simplex_tree.dimension() << std::endl;
  BOOST_CHECK(simplex_tree.dimension() == -1);
  
  std::clog << "simplex_tree.num_vertices()=" << simplex_tree.num_vertices() << std::endl;
  BOOST_CHECK(simplex_tree.num_vertices() == points.size());
}
