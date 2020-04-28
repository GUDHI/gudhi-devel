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
#include <vector>

using Persistence_landscape_on_grid = Gudhi::Persistence_representations::Persistence_landscape_on_grid;

int main(int argc, char** argv) {
  std::clog << "This program computes average of persistence landscapes on grid stored in files (the files needs to "
            << "be created beforehand).\n"
            << "The parameters of this programs are names of files with persistence landscapes on grid.\n";

  if (argc < 3) {
    std::clog << "Wrong number of parameters, the program will now terminate \n";
    return 1;
  }

  std::vector<const char*> filenames;
  for (int i = 1; i < argc; ++i) {
    filenames.push_back(argv[i]);
  }

  std::vector<Persistence_landscape_on_grid*> lands;
  for (size_t i = 0; i != filenames.size(); ++i) {
    Persistence_landscape_on_grid* l = new Persistence_landscape_on_grid;
    l->load_landscape_from_file(filenames[i]);
    lands.push_back(l);
  }

  Persistence_landscape_on_grid av;
  av.compute_average(lands);

  av.print_to_file("average.g_land");

  for (size_t i = 0; i != filenames.size(); ++i) {
    delete lands[i];
  }

  std::clog << "Average can be found in 'average.g_land' file\n";
  return 0;
}
