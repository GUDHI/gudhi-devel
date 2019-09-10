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

using Euclidean_distance = Gudhi::Euclidean_distance;
using Vector_distances_in_diagram = Gudhi::Persistence_representations::Vector_distances_in_diagram<Euclidean_distance>;

int main(int argc, char** argv) {
  std::cout << "This program create a Gnuplot script to plot persistence vector. Please call this program with the "
               "name of file with persistence vector. \n";
  if (argc != 2) {
    std::cout << "Wrong number of parameters, the program will now terminate. \n";
    return 1;
  }
  Vector_distances_in_diagram l;
  l.load_from_file(argv[1]);
  l.plot(argv[1]);

  return 0;
}
