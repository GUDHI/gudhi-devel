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

#define BOOST_TEST_MODULE alpha_shapes
#include <boost/test/included/unit_test.hpp>
#include <boost/system/error_code.hpp>
#include <boost/chrono/thread_clock.hpp>
// to construct a Delaunay_triangulation from a OFF file
#include "gudhi/Alpha_shapes/Delaunay_triangulation_off_io.h"
#include "gudhi/Alpha_shapes.h"

// to construct a simplex_tree from Delaunay_triangulation
#include "gudhi/graph_simplicial_complex.h"
#include "gudhi/Simplex_tree.h"

#include <CGAL/Delaunay_triangulation.h>
#include <CGAL/Epick_d.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/algorithm.h>
#include <CGAL/assertions.h>

#include <iostream>
#include <iterator>

#include <stdio.h>
#include <stdlib.h>
#include <string>

// Use dynamic_dimension_tag for the user to be able to set dimension
typedef CGAL::Epick_d< CGAL::Dynamic_dimension_tag > K;
typedef CGAL::Delaunay_triangulation<K> T;
// The triangulation uses the default instanciation of the 
// TriangulationDataStructure template parameter

BOOST_AUTO_TEST_CASE( OFF_file ) {
  // ----------------------------------------------------------------------------
  //
  // Init of an alpha-shape from a OFF file
  //
  // ----------------------------------------------------------------------------
  std::string off_file_name("S4_100.off");
  std::cout << "========== OFF FILE NAME = " << off_file_name << " ==========" << std::endl;
  
  Gudhi::alphashapes::Alpha_shapes alpha_shapes_from_file(off_file_name, 4);

  const int DIMENSION = 4;
  std::cout << "alpha_shapes_from_file.dimension()=" << alpha_shapes_from_file.dimension() << std::endl;
  BOOST_CHECK(alpha_shapes_from_file.dimension() == DIMENSION);
  
  const double FILTRATION = 0.0;
  std::cout << "alpha_shapes_from_file.filtration()=" << alpha_shapes_from_file.filtration() << std::endl;
  BOOST_CHECK(alpha_shapes_from_file.filtration() == FILTRATION);
  
  const int NUMBER_OF_VERTICES = 100;
  std::cout << "alpha_shapes_from_file.num_vertices()=" << alpha_shapes_from_file.num_vertices() << std::endl;
  BOOST_CHECK(alpha_shapes_from_file.num_vertices() == NUMBER_OF_VERTICES);
  
  const int NUMBER_OF_SIMPLICES = 6779;
  std::cout << "alpha_shapes_from_file.num_simplices()=" << alpha_shapes_from_file.num_simplices() << std::endl;
  BOOST_CHECK(alpha_shapes_from_file.num_simplices() == NUMBER_OF_SIMPLICES);

}

BOOST_AUTO_TEST_CASE( Delaunay_triangulation ) {
  // ----------------------------------------------------------------------------
  //
  // Init of an alpha-shape from a Delauny triangulation
  //
  // ----------------------------------------------------------------------------
  T dt(8);
  std::string off_file_name("S8_10.off");
  std::cout << "========== OFF FILE NAME = " << off_file_name << " ==========" << std::endl;
  
  Gudhi::alphashapes::Delaunay_triangulation_off_reader<T> off_reader(off_file_name, dt, true, true);
  std::cout << "off_reader.is_valid()=" << off_reader.is_valid() << std::endl;
  BOOST_CHECK(off_reader.is_valid());

  const int NUMBER_OF_VERTICES = 10;
  std::cout << "dt.number_of_vertices()=" << dt.number_of_vertices() << std::endl;
  BOOST_CHECK(dt.number_of_vertices() == NUMBER_OF_VERTICES);
  
  const int NUMBER_OF_FULL_CELLS = 30;
  std::cout << "dt.number_of_full_cells()=" << dt.number_of_full_cells() << std::endl;
  BOOST_CHECK(dt.number_of_full_cells() == NUMBER_OF_FULL_CELLS);
  
  const int NUMBER_OF_FINITE_FULL_CELLS = 6;
  std::cout << "dt.number_of_finite_full_cells()=" << dt.number_of_finite_full_cells() << std::endl;
  BOOST_CHECK(dt.number_of_finite_full_cells() == NUMBER_OF_FINITE_FULL_CELLS);

  Gudhi::alphashapes::Alpha_shapes alpha_shapes_from_dt(dt);

  const int DIMENSION = 8;
  std::cout << "alpha_shapes_from_dt.dimension()=" << alpha_shapes_from_dt.dimension() << std::endl;
  BOOST_CHECK(alpha_shapes_from_dt.dimension() == DIMENSION);
  
  const double FILTRATION = 0.0;
  std::cout << "alpha_shapes_from_dt.filtration()=" << alpha_shapes_from_dt.filtration() << std::endl;
  BOOST_CHECK(alpha_shapes_from_dt.filtration() == FILTRATION);
  
  std::cout << "alpha_shapes_from_dt.num_vertices()=" << alpha_shapes_from_dt.num_vertices() << std::endl;
  BOOST_CHECK(alpha_shapes_from_dt.num_vertices() == NUMBER_OF_VERTICES);
  
  const int NUMBER_OF_SIMPLICES = 997;
  std::cout << "alpha_shapes_from_dt.num_simplices()=" << alpha_shapes_from_dt.num_simplices() << std::endl;
  BOOST_CHECK(alpha_shapes_from_dt.num_simplices() == NUMBER_OF_SIMPLICES);
}

