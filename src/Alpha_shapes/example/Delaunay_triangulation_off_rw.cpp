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

#include <CGAL/Delaunay_triangulation.h>
#include <CGAL/Epick_d.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/algorithm.h>
#include <CGAL/assertions.h>

#include <iostream>
#include <iterator>
#include <vector>

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

  T dt(dimension);
  std::string offFileName(argv[1]);
  Gudhi::alphashapes::Delaunay_triangulation_off_reader<T> off_reader(offFileName, dt, true, true);
  if (!off_reader.is_valid()) {
    std::cerr << "Unable to read file " << offFileName << std::endl;
    exit(-1); // ----- >>
  }

  std::cout << "number of vertices=" << dt.number_of_vertices() << std::endl;
  std::cout << "number of full cells=" << dt.number_of_full_cells() << std::endl;
  std::cout << "number of finite full cells=" << dt.number_of_finite_full_cells() << std::endl;

  // Points list
  /*for (T::Vertex_iterator vit = dt.vertices_begin(); vit != dt.vertices_end(); ++vit) {
        for (auto Coord = vit->point().cartesian_begin(); Coord != vit->point().cartesian_end(); ++Coord) {
          std::cout << *Coord << " ";
        }
        std::cout << std::endl;
  }
  std::cout << std::endl;*/
 
  int i = 0, j = 0;
  typedef T::Full_cell_iterator Full_cell_iterator;
  typedef T::Facet Facet;

  for (Full_cell_iterator cit = dt.full_cells_begin(); cit != dt.full_cells_end(); ++cit) {
    if (!dt.is_infinite(cit)) {
      j++;
      continue;
    }
    Facet fct(cit, cit->index(dt.infinite_vertex()));
    i++;
  }
  std::cout << "There are " << i << " facets on the convex hull." << std::endl;
  std::cout << "There are " << j << " facets not on the convex hull." << std::endl;


  std::string offOutputFile("out.off");
  Gudhi::alphashapes::Delaunay_triangulation_off_writer<T> off_writer(offOutputFile, dt);

  return 0;
}