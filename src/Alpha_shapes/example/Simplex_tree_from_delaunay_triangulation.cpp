/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2014  INRIA Saclay (France)
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

void usage(char * const progName) {
  std::cerr << "Usage: " << progName << " filename.off" << std::endl;
  exit(-1); // ----- >>
}

int main(int argc, char **argv) {
  if (argc != 2) {
    std::cerr << "Error: Number of arguments (" << argc << ") is not correct" << std::endl;
    usage(argv[0]);
  }

  std::string off_file_name(argv[1]);

  // ----------------------------------------------------------------------------
  //
  // Init of an alpha-shape from a OFF file
  //
  // ----------------------------------------------------------------------------
  Gudhi::alphashapes::Alpha_shapes alpha_shapes_from_file(off_file_name);
  //std::cout << alpha_shapes_from_file << std::endl;
  
  std::cout << "alpha_shapes_from_file.dimension()=" << alpha_shapes_from_file.dimension() << std::endl;
  std::cout << "alpha_shapes_from_file.filtration()=" << alpha_shapes_from_file.filtration() << std::endl;
  std::cout << "alpha_shapes_from_file.num_simplices()=" << alpha_shapes_from_file.num_simplices() << std::endl;
  std::cout << "alpha_shapes_from_file.num_vertices()=" << alpha_shapes_from_file.num_vertices() << std::endl;

  return 0;
}