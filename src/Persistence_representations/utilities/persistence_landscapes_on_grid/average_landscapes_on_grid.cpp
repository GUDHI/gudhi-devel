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

#include <gudhi/Persistence_landscape_on_grid.h>

#include <iostream>
#include <vector>

using Persistence_landscape_on_grid = Gudhi::Persistence_representations::Persistence_landscape_on_grid;

int main(int argc, char** argv) {
  std::cout << "This program computes average of persistence landscapes on grid stored in files (the files needs to "
            << "be created beforehand).\n"
            << "The parameters of this programs are names of files with persistence landscapes on grid.\n";

  if (argc < 3) {
    std::cout << "Wrong number of parameters, the program will now terminate \n";
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

  std::cout << "Average can be found in 'average.g_land' file\n";
  return 0;
}
