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

#include <gudhi/Persistence_intervals.h>

#include <iostream>
#include <limits>
#include <vector>
#include <utility>

using Persistence_intervals = Gudhi::Persistence_representations::Persistence_intervals;

int main(int argc, char** argv) {
  std::cout << "This program compute the dominant intervals. A number of intervals to be displayed is a parameter of "
               "this program. \n";
  if (argc != 4) {
    std::cout << "To run this program, please provide the name of a file with persistence diagram, dimension of "
                 "intervals that should be taken into account (if your file contains only persistence pairs in a "
                 "single dimension, set it up to -1) and number of dominant intervals you would like to get \n";
    return 1;
  }
  int dim = atoi(argv[2]);
  unsigned dimension = std::numeric_limits<unsigned>::max();
  if (dim >= 0) {
    dimension = (unsigned)dim;
  }
  Persistence_intervals p(argv[1], dimension);
  std::vector<std::pair<double, double> > dominant_intervals = p.dominant_intervals(atoi(argv[3]));
  std::cout << "Here are the dominant intervals : " << std::endl;
  for (size_t i = 0; i != dominant_intervals.size(); ++i) {
    std::cout << " " << dominant_intervals[i].first << "," << dominant_intervals[i].second << " " << std::endl;
  }

  return 0;
}
