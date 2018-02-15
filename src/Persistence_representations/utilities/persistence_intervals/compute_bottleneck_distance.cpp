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

#include <iostream>
#include <sstream>
#include <limits>
#include <vector>

using Persistence_intervals_with_distances = Gudhi::Persistence_representations::Persistence_intervals_with_distances;

int main(int argc, char** argv) {
  std::cout << "This program computes the bottleneck distance of persistence pairs in diagrams provided as "
            << "an input.\n"
            << "The first parameter is the dimension of persistence to be used to create persistence intervals. "
            << "If your file contains the information about dimension of persistence pairs, please provide here the "
            << "dimension of persistence pairs you want to use. "
            << "If your input files consist only of birth-death pairs, please set this first parameter to -1.\n"
            << "The remaining parameters of the program are the names of files with persistence diagrams.\n";

  if (argc < 3) {
    std::cout << "Wrong number of parameters, the program will now terminate \n";
    return 1;
  }

  unsigned dimension = std::numeric_limits<unsigned>::max();
  int dim = atoi(argv[1]);
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
  out.open("distance.itv");
  for (size_t i = 0; i != distance.size(); ++i) {
    for (size_t j = 0; j != distance.size(); ++j) {
      std::cout << distance[i][j] << " ";
      out << distance[i][j] << " ";
    }
    std::cout << std::endl;
    out << std::endl;
  }
  out.close();

  std::cout << "Distance can be found in 'distance.itv' file\n";
  return 0;
}
