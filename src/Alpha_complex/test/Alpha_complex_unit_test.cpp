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
#define BOOST_TEST_MODULE "alpha_complex"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <CGAL/Epick_d.h>
#include <CGAL/Epeck_d.h>

#include <stdexcept> // std::out_of_range
#include <string>
#include <vector>

#include <gudhi/Alpha_complex.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Unitary_tests_utils.h>

// Use static dimension_tag for the user not to be able to set dimension
typedef CGAL::Epeck_d< CGAL::Dimension_tag<4> > Kernel_4;
typedef Kernel_4::Point_d Point_4;
typedef std::vector<Point_4> Vector_4_Points;

bool is_point_in_list(Vector_4_Points points_list, Point_4 point) {
  for (auto& point_in_list : points_list) {
    if (point_in_list == point) {
      return true;  // point found
    }
  }
  return false;  // point not found
}

BOOST_AUTO_TEST_CASE(Alpha_complex_from_points) {
  // ----------------------------------------------------------------------------
  // Init of a list of points
  // ----------------------------------------------------------------------------
  Vector_4_Points points;
  std::vector<double> coords = { 0.0, 0.0, 0.0, 1.0 };
  points.push_back(Point_4(coords.begin(), coords.end()));
  coords = { 0.0, 0.0, 1.0, 0.0 };
  points.push_back(Point_4(coords.begin(), coords.end()));
  coords = { 0.0, 1.0, 0.0, 0.0 };
  points.push_back(Point_4(coords.begin(), coords.end()));
  coords = { 1.0, 0.0, 0.0, 0.0 };
  points.push_back(Point_4(coords.begin(), coords.end()));

  // ----------------------------------------------------------------------------
  // Init of an alpha complex from the list of points
  // ----------------------------------------------------------------------------
  Gudhi::alpha_complex::Alpha_complex<Kernel_4> alpha_complex_from_points(points);

  std::clog << "========== Alpha_complex_from_points ==========" << std::endl;

  Gudhi::Simplex_tree<> simplex_tree;
  BOOST_CHECK(alpha_complex_from_points.create_complex(simplex_tree));
  
  std::clog << "alpha_complex_from_points.num_vertices()=" << alpha_complex_from_points.num_vertices() << std::endl;
  BOOST_CHECK(alpha_complex_from_points.num_vertices() == points.size());

  // Another way to check num_simplices
  std::clog << "Iterator on alpha complex simplices in the filtration order, with [filtration value]:" << std::endl;
  int num_simplices = 0;
  for (auto f_simplex : simplex_tree.filtration_simplex_range()) {
    num_simplices++;
    std::clog << "   ( ";
    for (auto vertex : simplex_tree.simplex_vertex_range(f_simplex)) {
      std::clog << vertex << " ";
    }
    std::clog << ") -> " << "[" << simplex_tree.filtration(f_simplex) << "] ";
    std::clog << std::endl;
  }
  BOOST_CHECK(num_simplices == 15);
  std::clog << "simplex_tree.num_simplices()=" << simplex_tree.num_simplices() << std::endl;
  BOOST_CHECK(simplex_tree.num_simplices() == 15);

  std::clog << "simplex_tree.dimension()=" << simplex_tree.dimension() << std::endl;
  BOOST_CHECK(simplex_tree.dimension() == 3);
  std::clog << "simplex_tree.num_vertices()=" << simplex_tree.num_vertices() << std::endl;
  BOOST_CHECK(simplex_tree.num_vertices() == points.size());

  for (auto f_simplex : simplex_tree.filtration_simplex_range()) {
    switch (simplex_tree.dimension(f_simplex)) {
      case 0:
        GUDHI_TEST_FLOAT_EQUALITY_CHECK(simplex_tree.filtration(f_simplex), 0.0);
        break;
      case 1:
        GUDHI_TEST_FLOAT_EQUALITY_CHECK(simplex_tree.filtration(f_simplex), 1.0/2.0);
        break;
      case 2:
        GUDHI_TEST_FLOAT_EQUALITY_CHECK(simplex_tree.filtration(f_simplex), 2.0/3.0);
        break;
      case 3:
        GUDHI_TEST_FLOAT_EQUALITY_CHECK(simplex_tree.filtration(f_simplex), 3.0/4.0);
        break;
      default:
        BOOST_CHECK(false);  // Shall not happen
        break;
    }
  }

  Point_4 p0 = alpha_complex_from_points.get_point(0);
  std::clog << "alpha_complex_from_points.get_point(0)=" << p0 << std::endl;
  BOOST_CHECK(4 == p0.dimension());
  BOOST_CHECK(is_point_in_list(points, p0));

  Point_4 p1 = alpha_complex_from_points.get_point(1);
  std::clog << "alpha_complex_from_points.get_point(1)=" << p1 << std::endl;
  BOOST_CHECK(4 == p1.dimension());
  BOOST_CHECK(is_point_in_list(points, p1));

  Point_4 p2 = alpha_complex_from_points.get_point(2);
  std::clog << "alpha_complex_from_points.get_point(2)=" << p2 << std::endl;
  BOOST_CHECK(4 == p2.dimension());
  BOOST_CHECK(is_point_in_list(points, p2));

  Point_4 p3 = alpha_complex_from_points.get_point(3);
  std::clog << "alpha_complex_from_points.get_point(3)=" << p3 << std::endl;
  BOOST_CHECK(4 == p3.dimension());
  BOOST_CHECK(is_point_in_list(points, p3));

  // Test to the limit
  BOOST_CHECK_THROW (alpha_complex_from_points.get_point(4), std::out_of_range);
  BOOST_CHECK_THROW (alpha_complex_from_points.get_point(-1), std::out_of_range);
  BOOST_CHECK_THROW (alpha_complex_from_points.get_point(1234), std::out_of_range);
  
  // Test after prune_above_filtration
  bool modified = simplex_tree.prune_above_filtration(0.6);
  BOOST_CHECK(modified);
  
  // Another way to check num_simplices
  std::clog << "Iterator on alpha complex simplices in the filtration order, with [filtration value]:" << std::endl;
  num_simplices = 0;
  for (auto f_simplex : simplex_tree.filtration_simplex_range()) {
    num_simplices++;
    std::clog << "   ( ";
    for (auto vertex : simplex_tree.simplex_vertex_range(f_simplex)) {
      std::clog << vertex << " ";
    }
    std::clog << ") -> " << "[" << simplex_tree.filtration(f_simplex) << "] ";
    std::clog << std::endl;
  }
  BOOST_CHECK(num_simplices == 10);
  std::clog << "simplex_tree.num_simplices()=" << simplex_tree.num_simplices() << std::endl;
  BOOST_CHECK(simplex_tree.num_simplices() == 10);

  std::clog << "simplex_tree.dimension()=" << simplex_tree.dimension() << std::endl;
  BOOST_CHECK(simplex_tree.dimension() == 1);
  std::clog << "simplex_tree.num_vertices()=" << simplex_tree.num_vertices() << std::endl;
  BOOST_CHECK(simplex_tree.num_vertices() == 4);

  for (auto f_simplex : simplex_tree.filtration_simplex_range()) {
    switch (simplex_tree.dimension(f_simplex)) {
      case 0:
        GUDHI_TEST_FLOAT_EQUALITY_CHECK(simplex_tree.filtration(f_simplex), 0.0);
        break;
      case 1:
        GUDHI_TEST_FLOAT_EQUALITY_CHECK(simplex_tree.filtration(f_simplex), 1.0/2.0);
        break;
      default:
        BOOST_CHECK(false);  // Shall not happen
        break;
    }
  }

}


using Inexact_kernel_2 = CGAL::Epick_d< CGAL::Dimension_tag<2> >;
using Exact_kernel_2 = CGAL::Epeck_d< CGAL::Dimension_tag<2> >;
using list_of_kernel_2_variants = boost::mpl::list<Inexact_kernel_2, Exact_kernel_2>;

BOOST_AUTO_TEST_CASE_TEMPLATE(Alpha_complex_with_duplicated_points, TestedKernel, list_of_kernel_2_variants) {
  std::clog << "========== Alpha_complex_with_duplicated_points ==========" << std::endl;

  using Point = typename TestedKernel::Point_d;
  using Vector_of_points = std::vector<Point>;

  // ----------------------------------------------------------------------------
  // Init of a list of points
  // ----------------------------------------------------------------------------
  Vector_of_points points;
  points.push_back(Point(1.0, 1.0));
  points.push_back(Point(7.0, 0.0));
  points.push_back(Point(4.0, 6.0));
  points.push_back(Point(9.0, 6.0));
  points.push_back(Point(0.0, 14.0));
  points.push_back(Point(2.0, 19.0));
  points.push_back(Point(9.0, 17.0));
  // duplicated points
  points.push_back(Point(1.0, 1.0));
  points.push_back(Point(7.0, 0.0));

  // ----------------------------------------------------------------------------
  // Init of an alpha complex from the list of points
  // ----------------------------------------------------------------------------
  std::clog << "Init" << std::endl;
  Gudhi::alpha_complex::Alpha_complex<TestedKernel> alpha_complex_from_points(points);

  Gudhi::Simplex_tree<> simplex_tree;
  std::clog << "create_complex" << std::endl;
  BOOST_CHECK(alpha_complex_from_points.create_complex(simplex_tree));
  
  std::clog << "alpha_complex_from_points.num_vertices()=" << alpha_complex_from_points.num_vertices() << std::endl;
  BOOST_CHECK(alpha_complex_from_points.num_vertices() < points.size());

  std::clog << "simplex_tree.num_vertices()=" << simplex_tree.num_vertices()
      << std::endl;
  BOOST_CHECK(simplex_tree.num_vertices() < points.size());

  int found = 0;
  for (int i = -1; i < (int)points.size() + 2; ++i)
    try {
      alpha_complex_from_points.get_point(i);
      ++found;
    } catch (...) {}
  BOOST_CHECK(found == simplex_tree.num_vertices());
}
