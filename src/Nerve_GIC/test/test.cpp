/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Mathieu Carrière
 *
 *    Copyright (C) 2017  INRIA Saclay (France)
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
#define BOOST_TEST_MODULE "graph_induced_complex"

#include <boost/test/unit_test.hpp>
#include <cmath>  // float comparison
#include <limits>
#include <string>
#include <vector>
#include <algorithm>    // std::max
#include <gudhi/GIC.h>
#include <gudhi/distance_functions.h>
#include <gudhi/reader_utils.h>

bool are_almost_the_same(float a, float b) {
  return std::fabs(a - b) < std::numeric_limits<float>::epsilon();
}

BOOST_AUTO_TEST_CASE(check_nerve) {

  Gudhi::graph_induced_complex::Graph_induced_complex GIC;
  std::string graph_file_name("data/graph"); GIC.set_graph_from_file(graph_file_name);
  std::string cover_file_name("data/cover"); GIC.set_cover_from_file(cover_file_name);
  std::string cloud_file_name("data/cloud"); GIC.set_color_from_file(cloud_file_name);
  GIC.find_Nerve_simplices(); Simplex_tree stree; GIC.create_complex(stree);

  BOOST_CHECK(stree.num_vertices() == 3);
  BOOST_CHECK((stree.num_simplices()-stree.num_vertices()) == 0);
  BOOST_CHECK(stree.dimension() == 0);
}

BOOST_AUTO_TEST_CASE(check_GICMAP) {

  Gudhi::graph_induced_complex::Graph_induced_complex GIC; GIC.set_verbose(1);
  std::string graph_file_name("data/graph"); GIC.set_graph_from_file(graph_file_name);
  std::string cover_file_name("data/cover"); GIC.set_cover_from_file(cover_file_name);
  std::string cloud_file_name("data/cloud"); GIC.set_color_from_file(cloud_file_name);
  GIC.find_GICMAP_simplices_with_functional_minimal_cover(); Simplex_tree stree; GIC.create_complex(stree);

  BOOST_CHECK(stree.num_vertices() == 3);
  BOOST_CHECK((stree.num_simplices()-stree.num_vertices()) == 2);
  BOOST_CHECK(stree.dimension() == 1);
}

BOOST_AUTO_TEST_CASE(check_GIC) {

  Gudhi::graph_induced_complex::Graph_induced_complex GIC;
  std::string graph_file_name("data/graph"); GIC.set_graph_from_file(graph_file_name);
  std::string cover_file_name("data/cover"); GIC.set_cover_from_file(cover_file_name);
  std::string cloud_file_name("data/cloud"); GIC.set_color_from_file(cloud_file_name);
  GIC.find_GIC_simplices(); Simplex_tree stree; GIC.create_complex(stree);

  BOOST_CHECK(stree.num_vertices() == 3);
  BOOST_CHECK((stree.num_simplices()-stree.num_vertices()) == 4);
  BOOST_CHECK(stree.dimension() == 2);
}


