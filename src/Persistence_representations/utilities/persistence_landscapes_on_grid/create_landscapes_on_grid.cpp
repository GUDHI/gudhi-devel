/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/Persistence_landscape_on_grid.h>

#include <iostream>
#include <sstream>
#include <limits>
#include <vector>

using Persistence_landscape_on_grid = Gudhi::Persistence_representations::Persistence_landscape_on_grid;

int main(int argc, char** argv) {
  std::clog << "This program creates persistence landscapes on grid files (*.g_land) of persistence diagrams files "
            << "(*.pers) provided as an input.\n"
            << "The first parameter of a program is an integer, a size of a grid.\n"
            << "The second and third parameters are min and max of the grid. If you want those numbers to be computed "
            << "based on the data, set them both to -1 \n"
            << "The fourth parameter of this program is a dimension of persistence that will be used in creation of "
            << "the persistence heat maps."
            << "If your input files contains persistence pairs of various dimension, as a fourth parameter of the "
            << "procedure please provide the dimension of persistence you want to use."
            << "If in your files there are only birth-death pairs of the same dimension, set the fourth parameter to "
            << "-1.\n"
            << "The remaining parameters are the names of files with persistence diagrams. \n";

  if (argc < 6) {
    std::clog << "Wrong parameter list, the program will now terminate \n";
    return 1;
  }

  size_t size_of_grid = (size_t)atoi(argv[1]);
  double min_ = atof(argv[2]);
  double max_ = atof(argv[3]);
  int dim = atoi(argv[4]);
  unsigned dimension = std::numeric_limits<unsigned>::max();
  if (dim >= 0) {
    dimension = (unsigned)dim;
  }

  std::vector<const char*> filenames;
  for (int i = 5; i < argc; ++i) {
    filenames.push_back(argv[i]);
  }

  for (size_t i = 0; i != filenames.size(); ++i) {
    std::clog << "Creating persistence landscape on a grid based on a file : " << filenames[i] << std::endl;
    Persistence_landscape_on_grid l;
    if ((min_ != -1) || (max_ != -1)) {
      l = Persistence_landscape_on_grid(filenames[i], min_, max_, size_of_grid, dimension);
    } else {
      // (min_ == -1) && (max_ == -1), in this case the program will find min_ and max_ based on the data.
      l = Persistence_landscape_on_grid(filenames[i], size_of_grid, dimension);
    }
    std::stringstream ss;
    ss << filenames[i] << ".g_land";
    l.print_to_file(ss.str().c_str());
  }
  return 0;
}
