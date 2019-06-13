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

#include <gudhi/Persistence_landscape.h>

#include <iostream>
#include <vector>

using Persistence_landscape = Gudhi::Persistence_representations::Persistence_landscape;

int main(int argc, char** argv) {
  std::cout << "This program computes average of persistence landscapes stored in files (the files needs to be "
            << "created beforehand).\n"
            << "The parameters of this programs are names of files with persistence landscapes.\n";
  std::vector<const char*> filenames;

  if (argc < 3) {
    std::cout << "Wrong number of parameters, the program will now terminate \n";
    return 1;
  }

  for (int i = 1; i < argc; ++i) {
    filenames.push_back(argv[i]);
  }

  std::vector<Persistence_landscape*> lands;
  for (size_t i = 0; i != filenames.size(); ++i) {
    Persistence_landscape* l = new Persistence_landscape;
    l->load_landscape_from_file(filenames[i]);
    lands.push_back(l);
  }

  Persistence_landscape av;
  av.compute_average(lands);

  av.print_to_file("average.land");

  for (size_t i = 0; i != filenames.size(); ++i) {
    delete lands[i];
  }

  std::cout << "Average can be found in 'average.land' file\n";
  return 0;
}
