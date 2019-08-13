/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/Persistence_vectors.h>

#include <iostream>
#include <sstream>
#include <limits>
#include <vector>

using Euclidean_distance = Gudhi::Euclidean_distance;
using Vector_distances_in_diagram = Gudhi::Persistence_representations::Vector_distances_in_diagram<Euclidean_distance>;

int main(int argc, char** argv) {
  std::cout << "This program creates persistence vectors files (*.vect) of persistence diagrams files (*.pers) "
            << "provided as an input.\n"
            << "The first parameter of this program is a dimension of persistence that will be used in creation of "
            << "the persistence heat maps."
            << "If your input files contains persistence pairs of various dimension, as a first parameter of the "
            << "procedure please provide the dimension of persistence you want to use."
            << "If in your files there are only birth-death pairs of the same dimension, set the first parameter to "
            << "-1.\n"
            << "The remaining parameters are the names of files with persistence diagrams. \n";

  if (argc < 3) {
    std::cout << "Wrong parameter list, the program will now terminate \n";
    return 1;
  }

  std::cout << "The remaining parameters are the names of files with persistence diagrams. \n";
  int dim = atoi(argv[1]);
  unsigned dimension = std::numeric_limits<unsigned>::max();
  if (dim >= 0) {
    dimension = (unsigned)dim;
  }

  std::vector<const char*> filenames;
  for (int i = 2; i < argc; ++i) {
    filenames.push_back(argv[i]);
  }

  for (size_t i = 0; i != filenames.size(); ++i) {
    std::cerr << "Creating persistence vectors based on a file : " << filenames[i] << std::endl;
    Vector_distances_in_diagram l(filenames[i], dimension);
    std::stringstream ss;
    ss << filenames[i] << ".vect";
    l.print_to_file(ss.str().c_str());
  }
  return 0;
}
