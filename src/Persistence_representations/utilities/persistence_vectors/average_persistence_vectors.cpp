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
#include <vector>

using Euclidean_distance = Gudhi::Euclidean_distance;
using Vector_distances_in_diagram = Gudhi::Persistence_representations::Vector_distances_in_diagram<Euclidean_distance>;

int main(int argc, char** argv) {
  std::cout << "This program computes average of persistence vectors stored in files (the files needs to "
            << "be created beforehand).\n"
            << "The parameters of this programs are names of files with persistence vectors.\n";

  if (argc < 3) {
    std::cout << "Wrong number of parameters, the program will now terminate \n";
    return 1;
  }

  std::vector<const char*> filenames;
  for (int i = 1; i < argc; ++i) {
    filenames.push_back(argv[i]);
  }

  std::vector<Vector_distances_in_diagram*> lands;
  for (size_t i = 0; i != filenames.size(); ++i) {
    Vector_distances_in_diagram* l = new Vector_distances_in_diagram;
    l->load_from_file(filenames[i]);
    lands.push_back(l);
  }

  Vector_distances_in_diagram av;
  av.compute_average(lands);

  av.print_to_file("average.vect");

  for (size_t i = 0; i != filenames.size(); ++i) {
    delete lands[i];
  }

  std::cout << "Done \n";

  return 0;
}
