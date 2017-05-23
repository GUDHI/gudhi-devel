/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Authors:       Francois Godi, small modifications by Pawel Dlotko
 *
 *    Copyright (C) 2015  INRIA
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

#define CGAL_HAS_THREADS

#include <gudhi/Bottleneck.h>
#include <gudhi/reader_utils.h>
#include <iostream>
#include <vector>
#include <utility>  // for pair
#include <fstream>
#include <sstream>
#include <string>

struct Persistence_interval
  : std::pair<double, double>
{
  Persistence_interval(std::tuple<int, double, double> data)
    : std::pair<double, double>(std::make_pair(std::get<1>(data), std::get<2>(data)))
  {}
};

int main(int argc, char** argv) {
  if (argc < 3) {
    std::cout << "To run this program please provide as an input two files with persistence diagrams. Each file " <<
        "should contain a birth-death pair per line. Third, optional parameter is an error bound on a bottleneck" <<
        " distance (set by default to zero). The program will now terminate \n";
  }
  std::vector<Persistence_interval> diag1;
  std::vector<Persistence_interval> diag2;
  read_persistence_diagram_from_file(argv[1], std::back_inserter(diag1));
  read_persistence_diagram_from_file(argv[2], std::back_inserter(diag2));

  double tolerance = 0.;
  if (argc == 4) {
    tolerance = atof(argv[3]);
  }
  double b = Gudhi::persistence_diagram::bottleneck_distance(diag1, diag2, tolerance);
  std::cout << "The distance between the diagrams is : " << b << ". The tolerance is : " << tolerance << std::endl;
}
