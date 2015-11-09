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

// to construct a Delaunay_triangulation from a OFF file
#include "gudhi/Delaunay_triangulation_off_io.h"

#include <CGAL/Delaunay_triangulation.h>
#include <CGAL/Epick_d.h>

#include <stdlib.h>

#include <iostream>
#include <string>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "delaunay_triangulation_off_read_write"
#include <boost/test/unit_test.hpp>

// Use dynamic_dimension_tag for the user to be able to set dimension
typedef CGAL::Epick_d< CGAL::Dynamic_dimension_tag > K;
typedef CGAL::Delaunay_triangulation<K> T;

BOOST_AUTO_TEST_CASE( Delaunay_triangulation_doc_test )
{
  // Read the OFF file (input file name given as parameter) and triangulates points
  Gudhi::Delaunay_triangulation_off_reader<T> off_reader("alphacomplexdoc.off");
  // Check the read operation was correct
  BOOST_CHECK(off_reader.is_valid());
  
  // Retrieve the triangulation
  T* triangulation = off_reader.get_complex();
  BOOST_CHECK(triangulation != nullptr);
  // Operations on triangulation
  BOOST_CHECK(triangulation->number_of_vertices() == 7);
  BOOST_CHECK(triangulation->number_of_finite_full_cells() == 6);

  // Write the OFF file (output file name given as parameter) with the points and triangulated cells as faces
  Gudhi::Delaunay_triangulation_off_writer<T> off_writer("UT.off", triangulation);

  // Check the write operation was correct
  BOOST_CHECK(off_writer.is_valid());
  
  delete triangulation;
}

BOOST_AUTO_TEST_CASE( Delaunay_triangulation_unexisting_file_read_test )
{
  Gudhi::Delaunay_triangulation_off_reader<T> off_reader("pouetpouet_tralala.off");
  // Check the read operation was correct
  BOOST_CHECK(!off_reader.is_valid());
  T* triangulation = off_reader.get_complex();
  BOOST_CHECK(triangulation == nullptr);
}

BOOST_AUTO_TEST_CASE( Delaunay_triangulation_unexisting_file_write_test )
{
  // Read the OFF file (input file name given as parameter) and triangulates points
  Gudhi::Delaunay_triangulation_off_reader<T> off_reader("alphacomplexdoc.off");
  
  // Retrieve the triangulation
  T* triangulation = off_reader.get_complex();

  // Write the OFF file (output file name given as parameter) with the points and triangulated cells as faces
  Gudhi::Delaunay_triangulation_off_writer<T> off_writer("/pouetpouet_tralala/pouetpouet_tralala/pouetpouet_tralala.off", triangulation);

  // Check the write operation was correct
  BOOST_CHECK(!off_writer.is_valid());
  
  delete triangulation;
}

