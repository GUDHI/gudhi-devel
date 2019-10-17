/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/Persistence_heat_maps.h>

#include <iostream>
#include <vector>

using constant_scaling_function = Gudhi::Persistence_representations::constant_scaling_function;
using Persistence_heat_maps = Gudhi::Persistence_representations::Persistence_heat_maps<constant_scaling_function>;

int main(int argc, char** argv) {
  std::cout << "This program computes average of persistence heat maps stored in files (the files needs to be "
            << "created beforehand).\n"
            << "The parameters of this programs are names of files with persistence heat maps.\n";

  if (argc < 3) {
    std::cout << "Wrong number of parameters, the program will now terminate \n";
    return 1;
  }

  std::vector<const char*> filenames;
  for (int i = 1; i < argc; ++i) {
    filenames.push_back(argv[i]);
  }

  std::vector<Persistence_heat_maps*> maps;
  for (size_t i = 0; i != filenames.size(); ++i) {
    Persistence_heat_maps* l = new Persistence_heat_maps;
    l->load_from_file(filenames[i]);
    maps.push_back(l);
  }

  Persistence_heat_maps av;
  av.compute_average(maps);
  av.print_to_file("average.mps");

  for (size_t i = 0; i != filenames.size(); ++i) {
    delete maps[i];
  }

  std::cout << "Average can be found in 'average.mps' file\n";
  return 0;
}
