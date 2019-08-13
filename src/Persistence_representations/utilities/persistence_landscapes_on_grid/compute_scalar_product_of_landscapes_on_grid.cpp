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

#include <gudhi/Persistence_landscape_on_grid.h>

#include <iostream>
#include <sstream>
#include <vector>

using Persistence_landscape_on_grid = Gudhi::Persistence_representations::Persistence_landscape_on_grid;

int main(int argc, char** argv) {
  std::cout
      << "This program computes scalar product of persistence landscapes on grid stored in a file (the file needs to "
      << "be created beforehand). \n"
      << "The parameters of this programs are names of files with persistence landscapes on grid.\n";

  if (argc < 3) {
    std::cout << "Wrong number of parameters, the program will now terminate \n";
    return 1;
  }

  std::vector<const char*> filenames;
  for (int i = 1; i < argc; ++i) {
    filenames.push_back(argv[i]);
  }
  std::vector<Persistence_landscape_on_grid> landscaspes;
  landscaspes.reserve(filenames.size());
  for (size_t file_no = 0; file_no != filenames.size(); ++file_no) {
    Persistence_landscape_on_grid l;
    l.load_landscape_from_file(filenames[file_no]);
    landscaspes.push_back(l);
  }

  // and now we will compute the scalar product of landscapes.

  // first we prepare an array:
  std::vector<std::vector<double> > scalar_product(filenames.size());
  for (size_t i = 0; i != filenames.size(); ++i) {
    std::vector<double> v(filenames.size(), 0);
    scalar_product[i] = v;
  }

  // and now we can compute the scalar product:
  for (size_t i = 0; i != landscaspes.size(); ++i) {
    for (size_t j = i; j != landscaspes.size(); ++j) {
      scalar_product[i][j] = scalar_product[j][i] = compute_inner_product(landscaspes[i], landscaspes[j]);
    }
  }

  // and now output the result to the screen and a file:
  std::ofstream out;
  out.open("scalar_product.g_land");
  for (size_t i = 0; i != scalar_product.size(); ++i) {
    for (size_t j = 0; j != scalar_product.size(); ++j) {
      std::cout << scalar_product[i][j] << " ";
      out << scalar_product[i][j] << " ";
    }
    std::cout << std::endl;
    out << std::endl;
  }
  out.close();

  std::cout << "Distance can be found in 'scalar_product.g_land' file\n";
  return 0;
}
