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
#include "gudhi/Delaunay_triangulation_off_io.h"
#include "gudhi/Alpha_complex.h"

// to construct a simplex_tree from Delaunay_triangulation
#include "gudhi/graph_simplicial_complex.h"
#include "gudhi/Simplex_tree.h"

#include <iostream>
#include <iterator>

#include <stdio.h>
#include <stdlib.h>
#include <string>

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
  // Init of an alpha-complex from an OFF file
  //
  // ----------------------------------------------------------------------------
  Gudhi::alphacomplex::Alpha_complex alpha_complex_from_file(off_file_name);
  
  std::cout << "alpha_complex_from_file.dimension()=" << alpha_complex_from_file.dimension() << std::endl;
  std::cout << "alpha_complex_from_file.filtration()=" << alpha_complex_from_file.filtration() << std::endl;
  std::cout << "alpha_complex_from_file.num_simplices()=" << alpha_complex_from_file.num_simplices() << std::endl;
  std::cout << "alpha_complex_from_file.num_vertices()=" << alpha_complex_from_file.num_vertices() << std::endl;
  
  std::cout << "Iterator on Simplices in the filtration order, with [filtration value]:" << std::endl;
  for (auto f_simplex : alpha_complex_from_file.filtration_simplex_range()) {
    std::cout << "   " << "[" << alpha_complex_from_file.filtration(f_simplex) << "] ";
    for (auto vertex : alpha_complex_from_file.simplex_vertex_range(f_simplex)) {
      std::cout << vertex << " ";
    }
    std::cout << std::endl;
  }
  return 0;
}