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
#include <sstream>
#include <limits>
#include <vector>

using constant_scaling_function = Gudhi::Persistence_representations::constant_scaling_function;
using Persistence_heat_maps = Gudhi::Persistence_representations::Persistence_heat_maps<constant_scaling_function>;

int main(int argc, char** argv) {
  std::clog << "This program computes distance of persistence heat maps stored in files (the files needs to be "
            << "created beforehand).\n"
            << "The first parameter of a program is an integer p. The program compute L^p distance of the two heat "
            << "maps. For L^infty distance choose p = -1. \n"
            << "The remaining parameters of this program are names of files with persistence heat maps.\n";

  if (argc < 3) {
    std::clog << "Wrong number of parameters, the program will now terminate \n";
    return 1;
  }

  int pp = atoi(argv[1]);
  double p = std::numeric_limits<double>::max();
  if (pp != -1) {
    p = pp;
  }

  std::vector<const char*> filenames;
  for (int i = 2; i < argc; ++i) {
    filenames.push_back(argv[i]);
  }
  std::vector<Persistence_heat_maps> maps;
  maps.reserve(filenames.size());
  for (size_t file_no = 0; file_no != filenames.size(); ++file_no) {
    Persistence_heat_maps l;
    l.load_from_file(filenames[file_no]);
    maps.push_back(l);
  }

  // and now we will compute the scalar product of landscapes.

  // first we prepare an array:
  std::vector<std::vector<double> > distance(filenames.size());
  for (size_t i = 0; i != filenames.size(); ++i) {
    std::vector<double> v(filenames.size(), 0);
    distance[i] = v;
  }

  // and now we can compute the distances:
  for (size_t i = 0; i != filenames.size(); ++i) {
    for (size_t j = i; j != filenames.size(); ++j) {
      distance[i][j] = distance[j][i] = maps[i].distance(maps[j], p);
    }
  }

  // and now output the result to the screen and a file:
  std::ofstream out;
  out.open("distance.mps");
  for (size_t i = 0; i != distance.size(); ++i) {
    for (size_t j = 0; j != distance.size(); ++j) {
      std::clog << distance[i][j] << " ";
      out << distance[i][j] << " ";
    }
    std::clog << std::endl;
    out << std::endl;
  }
  out.close();

  std::clog << "Distance can be found in 'distance.mps' file\n";
  return 0;
}
