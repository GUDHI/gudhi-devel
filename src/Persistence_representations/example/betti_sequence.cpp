/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Mathieu Carriere
 *
 *    Copyright (C) 2018  INRIA (France)
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

#include <gudhi/Betti_sequence.h>

#include <iostream>
#include <vector>
#include <utility>

using Persistence_diagram = Gudhi::Persistence_representations::Persistence_diagram;
using BS = Gudhi::Persistence_representations::Betti_sequence;

int main(int argc, char** argv) {

  Persistence_diagram persistence;

  persistence.push_back(std::make_pair(1, 2));
  persistence.push_back(std::make_pair(6, 8));
  persistence.push_back(std::make_pair(0, 4));
  persistence.push_back(std::make_pair(3, 8));

  double min_x = 0; double max_x = 8; int res_x = 1000;

  BS bs(persistence, min_x, max_x, res_x);
  std::vector<int> B = bs.vectorize();

  for(int i = 0; i < res_x; i++)  std::cout << B[i] << ", ";

  return 0;
}
