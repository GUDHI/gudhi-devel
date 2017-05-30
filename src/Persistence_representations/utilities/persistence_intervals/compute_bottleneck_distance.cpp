/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2016  INRIA (France)
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

#include <gudhi/Persistence_intervals_with_distances.h>
#include <gudhi/read_persistence_from_file.h>

using namespace Gudhi;
using namespace Gudhi::Persistence_representations;

#include <iostream>
#include <sstream>

int main(int argc, char** argv) {
  std::cout << "This program compute the bottleneck distance of persistence diagrams stored in a files. \n";
  std::cout << "The first parameter of the program is the dimension of persistence to be used to construct persistence "
               "landscapes. If your file contains ";
  std::cout << "the information about dimension of persistence pairs, please provide here the dimension of persistence "
               "pairs you want to use. If your input files consist only ";
  std::cout << "of birth-death pairs, please set this first parameter to -1 \n";
  std::cout << "The remaining parameters of this programs are names of files with persistence diagrams.\n";

  if (argc < 3) {
    std::cout << "Wrong number of parameters, the program will now terminate \n";
    return 1;
  }

  int dim = atoi(argv[1]);
  unsigned dimension = std::numeric_limits<unsigned>::max();
  if (dim >= 0) {
    dimension = (unsigned)dim;
  }

  std::vector<const char*> filenames;
  for (int i = 2; i < argc; ++i) {
    filenames.push_back(argv[i]);
  }

  // reading the persistence intervals:
  std::vector<Persistence_intervals_with_distances> persistence_intervals;
  for (size_t i = 0; i != filenames.size(); ++i) {
    Persistence_intervals_with_distances pers(filenames[i], dimension);
    persistence_intervals.push_back(pers);
  }

  // and now we will compute the scalar product of landscapes.

  // first we prepare an array:
  std::vector<std::vector<double> > distance(filenames.size());
  for (size_t i = 0; i != filenames.size(); ++i) {
    std::vector<double> v(filenames.size(), 0);
    distance[i] = v;
  }

  // and now we can compute the distances:
  for (size_t i = 0; i != persistence_intervals.size(); ++i) {
    for (size_t j = i + 1; j != persistence_intervals.size(); ++j) {
      distance[i][j] = distance[j][i] = persistence_intervals[i].distance(persistence_intervals[j]);
    }
  }

  // and now output the result to the screen and a file:
  std::ofstream out;
  out.open("distance");
  for (size_t i = 0; i != distance.size(); ++i) {
    for (size_t j = 0; j != distance.size(); ++j) {
      std::cout << distance[i][j] << " ";
      out << distance[i][j] << " ";
    }
    std::cout << std::endl;
    out << std::endl;
  }
  out.close();

  return 0;
}
