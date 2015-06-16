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
  std::cerr << "Usage: " << progName << " inputFile.off dimension outputFile.off" << std::endl;
  exit(-1); // ----- >>
}

int main(int argc, char **argv) {
  if (argc != 4) {
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
  Gudhi::alphacomplex::Delaunay_triangulation_off_reader<T> off_reader(offFileName, dt);
  if (!off_reader.is_valid()) {
    std::cerr << "Unable to read file " << offFileName << std::endl;
    exit(-1); // ----- >>
  }
  
  std::cout << "number_of_finite_full_cells= " << dt.number_of_finite_full_cells() << std::endl;

  std::string outFileName(argv[3]);
  std::string offOutputFile(outFileName);
  Gudhi::alphacomplex::Delaunay_triangulation_off_writer<T> off_writer(offOutputFile, dt);

  return 0;
}