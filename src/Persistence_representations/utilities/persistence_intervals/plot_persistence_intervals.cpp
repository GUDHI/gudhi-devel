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

#include <gudhi/reader_utils.h>
#include <gudhi/Persistence_intervals.h>
#include <gudhi/read_persistence_from_file.h>

#include <iostream>
#include <limits>
#include <vector>
#include <utility>

using namespace Gudhi;
using namespace Gudhi::Persistence_representations;

double epsilon = 0.0000005;

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cout << "To run this program, please provide the name of a file with persistence diagram \n";
    std::cout << "The second optional parameter of a program is the dimension of the persistence that is to be used. "
                 "If your file contains only birth-death pairs, you can skip this parameter\n";
    return 1;
  }
  unsigned dimension = std::numeric_limits<unsigned>::max();
  int dim = -1;
  if (argc > 2) {
    dim = atoi(argv[2]);
  }
  if (dim >= 0) {
    dimension = (unsigned)dim;
  }
  std::vector<std::pair<double, double> > intervals =
      read_persistence_intervals_in_one_dimension_from_file(argv[1], dimension);
  Persistence_intervals b(intervals);
  b.plot(argv[1]);
  return 0;
}
