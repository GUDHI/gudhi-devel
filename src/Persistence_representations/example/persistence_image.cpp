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

#include <gudhi/Persistence_image.h>
#include <gudhi/Persistence_weighted_gaussian.h>

#include <iostream>
#include <vector>
#include <utility>
#include <string>

using PI = Gudhi::Persistence_representations::Persistence_image;
using Weight = std::function<double (std::pair<double,double>) >;

int main(int argc, char** argv) {

  std::vector<std::pair<double, double> > persistence;

  persistence.push_back(std::make_pair(1, 2));
  persistence.push_back(std::make_pair(6, 8));
  persistence.push_back(std::make_pair(0, 4));
  persistence.push_back(std::make_pair(3, 8));

  double min_x = 0.0; double max_x = 10.0; int res_x = 100; double min_y = 0.0; double max_y = 10.0; int res_y = 100; double sigma = 1.0; Weight weight = Gudhi::Persistence_representations::linear_weight;

  PI pim(persistence, min_x, max_x, res_x, min_y, max_y, res_y, weight, sigma);
  std::vector<std::vector<double> > P = pim.vectorize();

  for(int i = 0; i < res_y; i++){
    for(int j = 0; j < res_x; j++) std::cout << P[i][j] << " ";
    std::cout << std::endl;
  }

  return 0;
}
