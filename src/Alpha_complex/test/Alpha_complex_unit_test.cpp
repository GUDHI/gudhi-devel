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

#define BOOST_TEST_MODULE alpha_complex
#include <boost/test/included/unit_test.hpp>
#include <boost/system/error_code.hpp>
#include <boost/chrono/thread_clock.hpp>
// to construct a Delaunay_triangulation from a OFF file
#include "gudhi/Delaunay_triangulation_off_io.h"
#include "gudhi/Alpha_complex.h"

#include <CGAL/Delaunay_triangulation.h>
#include <CGAL/Epick_d.h>

#include <cmath> // float comparison
#include <limits>

// Use dynamic_dimension_tag for the user to be able to set dimension
typedef CGAL::Epick_d< CGAL::Dynamic_dimension_tag > Kernel;
typedef Kernel::Point_d Point;
typedef std::vector<Point> Vector_of_points;
// The triangulation uses the default instantiation of the TriangulationDataStructure template parameter

BOOST_AUTO_TEST_CASE(S4_100_OFF_file) {
  // ----------------------------------------------------------------------------
  //
  // Init of an alpha-complex from a OFF file
  //
  // ----------------------------------------------------------------------------
  std::string off_file_name("S4_100.off");
  std::cout << "========== OFF FILE NAME = " << off_file_name << " ==========" << std::endl;

  Gudhi::alphacomplex::Alpha_complex alpha_complex_from_file(off_file_name);

  const int DIMENSION = 4;
  std::cout << "alpha_complex_from_file.dimension()=" << alpha_complex_from_file.dimension() << std::endl;
  BOOST_CHECK(alpha_complex_from_file.dimension() == DIMENSION);

  const int NUMBER_OF_VERTICES = 100;
  std::cout << "alpha_complex_from_file.num_vertices()=" << alpha_complex_from_file.num_vertices() << std::endl;
  BOOST_CHECK(alpha_complex_from_file.num_vertices() == NUMBER_OF_VERTICES);

  const int NUMBER_OF_SIMPLICES = 6879;
  std::cout << "alpha_complex_from_file.num_simplices()=" << alpha_complex_from_file.num_simplices() << std::endl;
  BOOST_CHECK(alpha_complex_from_file.num_simplices() == NUMBER_OF_SIMPLICES);

}

BOOST_AUTO_TEST_CASE(S8_10_OFF_file) {
  // ----------------------------------------------------------------------------
  //
  // Init of an alpha-complex from a OFF file
  //
  // ----------------------------------------------------------------------------
  std::string off_file_name("S8_10.off");
  std::cout << "========== OFF FILE NAME = " << off_file_name << " ==========" << std::endl;

  Gudhi::alphacomplex::Alpha_complex alpha_complex_from_file(off_file_name);

  const int DIMENSION = 8;
  std::cout << "alpha_complex_from_file.dimension()=" << alpha_complex_from_file.dimension() << std::endl;
  BOOST_CHECK(alpha_complex_from_file.dimension() == DIMENSION);

  const int NUMBER_OF_VERTICES = 10;
  std::cout << "alpha_complex_from_file.num_vertices()=" << alpha_complex_from_file.num_vertices() << std::endl;
  BOOST_CHECK(alpha_complex_from_file.num_vertices() == NUMBER_OF_VERTICES);

  const int NUMBER_OF_SIMPLICES = 1007;
  std::cout << "alpha_complex_from_file.num_simplices()=" << alpha_complex_from_file.num_simplices() << std::endl;
  BOOST_CHECK(alpha_complex_from_file.num_simplices() == NUMBER_OF_SIMPLICES);
}

bool are_almost_the_same(float a, float b) {
  return std::fabs(a - b) < std::numeric_limits<float>::epsilon();
}

BOOST_AUTO_TEST_CASE(Alpha_complex_from_points) {

  // ----------------------------------------------------------------------------
  // Init of a list of points
  // ----------------------------------------------------------------------------
  Vector_of_points points;
  std::vector<double> coords;

  coords.clear();
  coords.push_back(0.0);
  coords.push_back(0.0);
  coords.push_back(0.0);
  coords.push_back(1.0);
  points.push_back(Point(coords.begin(), coords.end()));
  coords.clear();
  coords.push_back(0.0);
  coords.push_back(0.0);
  coords.push_back(1.0);
  coords.push_back(0.0);
  points.push_back(Point(coords.begin(), coords.end()));
  coords.clear();
  coords.push_back(0.0);
  coords.push_back(1.0);
  coords.push_back(0.0);
  coords.push_back(0.0);
  points.push_back(Point(coords.begin(), coords.end()));
  coords.clear();
  coords.push_back(1.0);
  coords.push_back(0.0);
  coords.push_back(0.0);
  coords.push_back(0.0);
  points.push_back(Point(coords.begin(), coords.end()));

  // ----------------------------------------------------------------------------
  // Init of an alpha complex from the list of points
  // ----------------------------------------------------------------------------
  Gudhi::alphacomplex::Alpha_complex alpha_complex_from_points(3, points.size(), points.begin(), points.end());

  std::cout << "========== Alpha_complex_from_points ==========" << std::endl;

  std::cout << "alpha_complex_from_points.dimension()=" << alpha_complex_from_points.dimension() << std::endl;
  BOOST_CHECK(alpha_complex_from_points.dimension() == 3);
  std::cout << "alpha_complex_from_points.num_simplices()=" << alpha_complex_from_points.num_simplices() << std::endl;
  BOOST_CHECK(alpha_complex_from_points.num_simplices() == 15);
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
        BOOST_CHECK(false); // Shall not happen
        break;
    }
  }

}
