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
  std::cerr << "Usage: " << progName << " filename.off dimension" << std::endl;
  exit(-1); // ----- >>
}

int main(int argc, char **argv) {
  if (argc != 3) {
    std::cerr << "Error: Number of arguments (" << argc << ") is not correct" << std::endl;
    usage(argv[0]);
  }

  int dimension = 0;
  int returnedScanValue = sscanf(argv[2], "%d", &dimension);
  if ((returnedScanValue == EOF) || (dimension <= 0)) {
    std::cerr << "Error: " << argv[2] << " is not correct" << std::endl;
    usage(argv[0]);
  }

  // ----------------------------------------------------------------------------
  //
  // Init of an alpha-shape from a Delaunay triangulation
  //
  // ----------------------------------------------------------------------------
  T dt(dimension);
  std::string off_file_name(argv[1]);

  Gudhi::alphashapes::Delaunay_triangulation_off_reader<T> off_reader(off_file_name, dt, false, false);
  if (!off_reader.is_valid()) {
    std::cerr << "Unable to read file " << off_file_name << std::endl;
    exit(-1); // ----- >>
  }

  std::cout << "number of vertices=" << dt.number_of_vertices() << std::endl;
  std::cout << "number of full cells=" << dt.number_of_full_cells() << std::endl;
  std::cout << "number of finite full cells=" << dt.number_of_finite_full_cells() << std::endl;

  Gudhi::alphashapes::Alpha_shapes alpha_shapes_from_dt(dt);
  //std::cout << alpha_shapes_from_dt << std::endl;

  // ----------------------------------------------------------------------------
  //
  // Init of an alpha-shape from a OFF file
  //
  // ----------------------------------------------------------------------------
  Gudhi::alphashapes::Alpha_shapes alpha_shapes_from_file(off_file_name, dimension);
  //std::cout << alpha_shapes_from_file << std::endl;
  
  std::cout << "alpha_shapes_from_file.dimension()=" << alpha_shapes_from_file.dimension() << std::endl;
  std::cout << "alpha_shapes_from_file.filtration()=" << alpha_shapes_from_file.filtration() << std::endl;
  std::cout << "alpha_shapes_from_file.num_simplices()=" << alpha_shapes_from_file.num_simplices() << std::endl;
  std::cout << "alpha_shapes_from_file.num_vertices()=" << alpha_shapes_from_file.num_vertices() << std::endl;

  return 0;
}