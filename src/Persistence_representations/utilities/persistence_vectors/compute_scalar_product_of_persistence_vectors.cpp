/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
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
  std::clog << "This program computes scalar product of persistence vectors stored in a file (the file needs to "
            << "be created beforehand). \n"
            << "The parameters of this programs are names of files with persistence vectors.\n";

  if (argc < 3) {
    std::clog << "Wrong number of parameters, the program will now terminate \n";
    return 1;
  }

  std::vector<const char*> filenames;
  for (int i = 1; i < argc; ++i) {
    filenames.push_back(argv[i]);
  }
  std::vector<Vector_distances_in_diagram> vectors;
  vectors.reserve(filenames.size());
  for (size_t file_no = 0; file_no != filenames.size(); ++file_no) {
    Vector_distances_in_diagram l;
    l.load_from_file(filenames[file_no]);
    vectors.push_back(l);
  }

  // and now we will compute the scalar product of landscapes.

  // first we prepare an array:
  std::vector<std::vector<double> > scalar_product(filenames.size());
  for (size_t i = 0; i != filenames.size(); ++i) {
    std::vector<double> v(filenames.size(), 0);
    scalar_product[i] = v;
  }

  // and now we can compute the scalar product:
  for (size_t i = 0; i != vectors.size(); ++i) {
    for (size_t j = i; j != vectors.size(); ++j) {
      scalar_product[i][j] = scalar_product[j][i] = vectors[i].compute_scalar_product(vectors[j]);
    }
  }

  // and now output the result to the screen and a file:
  std::ofstream out;
  out.open("scalar_product.vect");
  for (size_t i = 0; i != scalar_product.size(); ++i) {
    for (size_t j = 0; j != scalar_product.size(); ++j) {
      std::clog << scalar_product[i][j] << " ";
      out << scalar_product[i][j] << " ";
    }
    std::clog << std::endl;
    out << std::endl;
  }
  out.close();

  std::clog << "Distance can be found in 'scalar_product.vect' file\n";
  return 0;
}
