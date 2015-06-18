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

#include <cmath> // float comparison

// Use dynamic_dimension_tag for the user to be able to set dimension
typedef CGAL::Epick_d< CGAL::Dynamic_dimension_tag > K;
typedef CGAL::Delaunay_triangulation<K> T;
// The triangulation uses the default instantiation of the 
// TriangulationDataStructure template parameter

BOOST_AUTO_TEST_CASE( S4_100_OFF_file ) {
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

BOOST_AUTO_TEST_CASE( S8_10_OFF_file ) {
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
