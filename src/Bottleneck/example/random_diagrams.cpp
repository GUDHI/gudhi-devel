/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Francois Godi
 *
 *    Copyright (C) 2015  INRIA Saclay (France)
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

#include "gudhi/Graph_matching.h"
#include <iostream>

using namespace Gudhi::bottleneck;

int main() {
  int n = 100;
  std::vector< std::pair<double, double> > v1, v2;
  for (int i = 0; i < n; i++) {
    int a = rand() % n;
    v1.emplace_back(a, a + rand() % (n - a));
    int b = rand() % n;
    v2.emplace_back(b, b + rand() % (n - b));
  }
  // v1 and v2 are persistence diagrams containing each 100 randoms points.
  double b = bottleneck_distance(v1, v2, 0);
  std::cout << b << std::endl;
}
