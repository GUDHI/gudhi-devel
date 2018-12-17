/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2016 Inria
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
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
