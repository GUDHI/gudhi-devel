/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2015  INRIA Saclay (France)
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
#define BOOST_TEST_MODULE "alpha_complex"
#include <boost/test/unit_test.hpp>

#include <CGAL/Delaunay_triangulation.h>
#include <CGAL/Epick_d.h>

#include <cmath>  // float comparison
#include <limits>
#include <string>
#include <vector>

#include <gudhi/Alpha_complex.h>

// Use dynamic_dimension_tag for the user to be able to set dimension
typedef CGAL::Epick_d< CGAL::Dynamic_dimension_tag > Kernel_d;
// The triangulation uses the default instantiation of the TriangulationDataStructure template parameter

BOOST_AUTO_TEST_CASE(ALPHA_DOC_OFF_file) {
  // ----------------------------------------------------------------------------
  //
  // Init of an alpha-complex from a OFF file
  //
  // ----------------------------------------------------------------------------
  std::string off_file_name("alphacomplexdoc.off");
  double max_alpha_square_value = 60.0;
  std::cout << "========== OFF FILE NAME = " << off_file_name << " - alpha²=" <<
      max_alpha_square_value << "==========" << std::endl;

  Gudhi::alphacomplex::Alpha_complex<Kernel_d> alpha_complex_from_file(off_file_name, max_alpha_square_value);

  const int DIMENSION = 2;
  std::cout << "alpha_complex_from_file.dimension()=" << alpha_complex_from_file.dimension() << std::endl;
  BOOST_CHECK(alpha_complex_from_file.dimension() == DIMENSION);

  const int NUMBER_OF_VERTICES = 7;
  std::cout << "alpha_complex_from_file.num_vertices()=" << alpha_complex_from_file.num_vertices() << std::endl;
  BOOST_CHECK(alpha_complex_from_file.num_vertices() == NUMBER_OF_VERTICES);

  const int NUMBER_OF_SIMPLICES = 25;
  std::cout << "alpha_complex_from_file.num_simplices()=" << alpha_complex_from_file.num_simplices() << std::endl;
  BOOST_CHECK(alpha_complex_from_file.num_simplices() == NUMBER_OF_SIMPLICES);

}

BOOST_AUTO_TEST_CASE(ALPHA_DOC_OFF_file_filtered) {
  // ----------------------------------------------------------------------------
  //
  // Init of an alpha-complex from a OFF file
  //
  // ----------------------------------------------------------------------------
  std::string off_file_name("alphacomplexdoc.off");
  double max_alpha_square_value = 59.0;
  std::cout << "========== OFF FILE NAME = " << off_file_name << " - alpha²=" <<
      max_alpha_square_value << "==========" << std::endl;

  // Use of the default dynamic kernel
  Gudhi::alphacomplex::Alpha_complex<> alpha_complex_from_file(off_file_name, max_alpha_square_value);

  const int DIMENSION = 2;
  std::cout << "alpha_complex_from_file.dimension()=" << alpha_complex_from_file.dimension() << std::endl;
  BOOST_CHECK(alpha_complex_from_file.dimension() == DIMENSION);

  const int NUMBER_OF_VERTICES = 7;
  std::cout << "alpha_complex_from_file.num_vertices()=" << alpha_complex_from_file.num_vertices() << std::endl;
  BOOST_CHECK(alpha_complex_from_file.num_vertices() == NUMBER_OF_VERTICES);

  const int NUMBER_OF_SIMPLICES = 23;
  std::cout << "alpha_complex_from_file.num_simplices()=" << alpha_complex_from_file.num_simplices() << std::endl;
  BOOST_CHECK(alpha_complex_from_file.num_simplices() == NUMBER_OF_SIMPLICES);
}

bool are_almost_the_same(float a, float b) {
  return std::fabs(a - b) < std::numeric_limits<float>::epsilon();
}

// Use dynamic_dimension_tag for the user to be able to set dimension
typedef CGAL::Epick_d< CGAL::Dimension_tag<4> > Kernel_s;
typedef Kernel_s::Point_d Point;
typedef std::vector<Point> Vector_of_points;


bool is_point_in_list(Vector_of_points points_list, Point point) {
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
  Vector_of_points points;
  std::vector<double> coords = { 0.0, 0.0, 0.0, 1.0 };
  points.push_back(Point(coords.begin(), coords.end()));
  coords = { 0.0, 0.0, 1.0, 0.0 };
  points.push_back(Point(coords.begin(), coords.end()));
  coords = { 0.0, 1.0, 0.0, 0.0 };
  points.push_back(Point(coords.begin(), coords.end()));
  coords = { 1.0, 0.0, 0.0, 0.0 };
  points.push_back(Point(coords.begin(), coords.end()));

  // ----------------------------------------------------------------------------
  // Init of an alpha complex from the list of points
  // ----------------------------------------------------------------------------
  Gudhi::alphacomplex::Alpha_complex<Kernel_s> alpha_complex_from_points(points);

  std::cout << "========== Alpha_complex_from_points ==========" << std::endl;

  // Another way to check num_simplices
  std::cout << "Iterator on alpha complex simplices in the filtration order, with [filtration value]:" << std::endl;
  int num_simplices = 0;
  for (auto f_simplex : alpha_complex_from_points.filtration_simplex_range()) {
    num_simplices++;
    std::cout << "   ( ";
    for (auto vertex : alpha_complex_from_points.simplex_vertex_range(f_simplex)) {
      std::cout << vertex << " ";
    }
    std::cout << ") -> " << "[" << alpha_complex_from_points.filtration(f_simplex) << "] ";
    std::cout << std::endl;
  }
  BOOST_CHECK(num_simplices == 15);
  std::cout << "alpha_complex_from_points.num_simplices()=" << alpha_complex_from_points.num_simplices() << std::endl;
  BOOST_CHECK(alpha_complex_from_points.num_simplices() == 15);

  std::cout << "alpha_complex_from_points.dimension()=" << alpha_complex_from_points.dimension() << std::endl;
  BOOST_CHECK(alpha_complex_from_points.dimension() == 4);
  std::cout << "alpha_complex_from_points.num_vertices()=" << alpha_complex_from_points.num_vertices() << std::endl;
  BOOST_CHECK(alpha_complex_from_points.num_vertices() == 4);

  for (auto f_simplex : alpha_complex_from_points.filtration_simplex_range()) {
    switch (alpha_complex_from_points.dimension(f_simplex)) {
      case 0:
        BOOST_CHECK(are_almost_the_same(alpha_complex_from_points.filtration(f_simplex), 0.0));
        break;
      case 1:
        BOOST_CHECK(are_almost_the_same(alpha_complex_from_points.filtration(f_simplex), 1.0/2.0));
        break;
      case 2:
        BOOST_CHECK(are_almost_the_same(alpha_complex_from_points.filtration(f_simplex), 2.0/3.0));
        break;
      case 3:
        BOOST_CHECK(are_almost_the_same(alpha_complex_from_points.filtration(f_simplex), 3.0/4.0));
        break;
      default:
        BOOST_CHECK(false);  // Shall not happen
        break;
    }
  }

  Point p0 = alpha_complex_from_points.get_point(0);
  std::cout << "alpha_complex_from_points.get_point(0)=" << p0 << std::endl;
  BOOST_CHECK(4 == p0.dimension());
  BOOST_CHECK(is_point_in_list(points, p0));

  Point p1 = alpha_complex_from_points.get_point(1);
  std::cout << "alpha_complex_from_points.get_point(1)=" << p1 << std::endl;
  BOOST_CHECK(4 == p1.dimension());
  BOOST_CHECK(is_point_in_list(points, p1));

  Point p2 = alpha_complex_from_points.get_point(2);
  std::cout << "alpha_complex_from_points.get_point(2)=" << p2 << std::endl;
  BOOST_CHECK(4 == p2.dimension());
  BOOST_CHECK(is_point_in_list(points, p2));

  Point p3 = alpha_complex_from_points.get_point(3);
  std::cout << "alpha_complex_from_points.get_point(3)=" << p3 << std::endl;
  BOOST_CHECK(4 == p3.dimension());
  BOOST_CHECK(is_point_in_list(points, p3));

  // Test to the limit
  BOOST_CHECK_THROW (alpha_complex_from_points.get_point(4), std::out_of_range);
  BOOST_CHECK_THROW (alpha_complex_from_points.get_point(-1), std::out_of_range);
  BOOST_CHECK_THROW (alpha_complex_from_points.get_point(1234), std::out_of_range);
  
  // Test after prune_above_filtration
  bool modified = alpha_complex_from_points.prune_above_filtration(0.6);
  if (modified) {
    alpha_complex_from_points.initialize_filtration();
  }
  BOOST_CHECK(modified);
  
  // Another way to check num_simplices
  std::cout << "Iterator on alpha complex simplices in the filtration order, with [filtration value]:" << std::endl;
  num_simplices = 0;
  for (auto f_simplex : alpha_complex_from_points.filtration_simplex_range()) {
    num_simplices++;
    std::cout << "   ( ";
    for (auto vertex : alpha_complex_from_points.simplex_vertex_range(f_simplex)) {
      std::cout << vertex << " ";
    }
    std::cout << ") -> " << "[" << alpha_complex_from_points.filtration(f_simplex) << "] ";
    std::cout << std::endl;
  }
  BOOST_CHECK(num_simplices == 10);
  std::cout << "alpha_complex_from_points.num_simplices()=" << alpha_complex_from_points.num_simplices() << std::endl;
  BOOST_CHECK(alpha_complex_from_points.num_simplices() == 10);

  std::cout << "alpha_complex_from_points.dimension()=" << alpha_complex_from_points.dimension() << std::endl;
  BOOST_CHECK(alpha_complex_from_points.dimension() == 4);
  std::cout << "alpha_complex_from_points.num_vertices()=" << alpha_complex_from_points.num_vertices() << std::endl;
  BOOST_CHECK(alpha_complex_from_points.num_vertices() == 4);

  for (auto f_simplex : alpha_complex_from_points.filtration_simplex_range()) {
    switch (alpha_complex_from_points.dimension(f_simplex)) {
      case 0:
        BOOST_CHECK(are_almost_the_same(alpha_complex_from_points.filtration(f_simplex), 0.0));
        break;
      case 1:
        BOOST_CHECK(are_almost_the_same(alpha_complex_from_points.filtration(f_simplex), 1.0/2.0));
        break;
      default:
        BOOST_CHECK(false);  // Shall not happen
        break;
    }
  }

}
