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

#include <gudhi/Persistence_heat_maps.h>

#include <iostream>
#include <sstream>
#include <limits>
#include <vector>

using distance_from_diagonal_scaling = Gudhi::Persistence_representations::distance_from_diagonal_scaling;
using Persistence_heat_maps = Gudhi::Persistence_representations::Persistence_heat_maps<distance_from_diagonal_scaling>;

int main(int argc, char** argv) {
  std::cout << "This program creates persistence heat map files (*.mps) of persistence diagrams files (*.pers) "
            << "provided as an input.The Gaussian kernels are weighted by the distance of a center from the "
            << "diagonal.\n"
            << "The first parameter of a program is an integer, a size of a grid.\n"
            << "The second and third parameters are min and max of the grid. If you want those numbers to be computed "
            << "based on the data, set them both to -1 \n"
            << "The fourth parameter is an integer, the standard deviation of a Gaussian kernel expressed in a number "
            << "of pixels.\n"
            << "The fifth parameter of this program is a dimension of persistence that will be used in creation of "
            << "the persistence heat maps."
            << "If your input files contains persistence pairs of various dimension, as a fifth parameter of the "
            << "procedure please provide the dimension of persistence you want to use."
            << "If in your files there are only birth-death pairs of the same dimension, set the fifth parameter to "
            << "-1.\n"
            << "The remaining parameters are the names of files with persistence diagrams. \n";

  if (argc < 7) {
    std::cout << "Wrong parameter list, the program will now terminate \n";
    return 1;
  }

  size_t size_of_grid = (size_t)atoi(argv[1]);
  double min_ = atof(argv[2]);
  double max_ = atof(argv[3]);
  size_t stdiv = atof(argv[4]);

  unsigned dimension = std::numeric_limits<unsigned>::max();
  int dim = atoi(argv[5]);
  if (dim >= 0) {
    dimension = (unsigned)dim;
  }

  std::vector<const char*> filenames;
  for (int i = 6; i != argc; ++i) {
    filenames.push_back(argv[i]);
  }

  std::vector<std::vector<double> > filter = Gudhi::Persistence_representations::create_Gaussian_filter(stdiv, 1);
  for (size_t i = 0; i != filenames.size(); ++i) {
    std::cout << "Creating a heat map based on a file : " << filenames[i] << std::endl;
    Persistence_heat_maps l(filenames[i], filter, false, size_of_grid, min_, max_, dimension);

    std::stringstream ss;
    ss << filenames[i] << ".mps";
    l.print_to_file(ss.str().c_str());
  }
  return 0;
}
