/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2015  INRIA (France)
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
#include <gudhi/Persistence_heat_maps.h>

#include <iostream>
#include <vector>

using namespace Gudhi;
using namespace Gudhi::Persistence_representations;

double epsilon = 0.0000005;

int main(int argc, char** argv) {
  // create two simple vectors with birth--death pairs:

  std::vector<std::pair<double, double> > persistence1;
  std::vector<std::pair<double, double> > persistence2;

  persistence1.push_back(std::make_pair(1, 2));
  persistence1.push_back(std::make_pair(6, 8));
  persistence1.push_back(std::make_pair(0, 4));
  persistence1.push_back(std::make_pair(3, 8));

  persistence2.push_back(std::make_pair(2, 9));
  persistence2.push_back(std::make_pair(1, 6));
  persistence2.push_back(std::make_pair(3, 5));
  persistence2.push_back(std::make_pair(6, 10));

  // over here we define a function we sill put on a top on every birth--death pair in the persistence interval. It can
  // be anything. Over here we will use standard Gaussian
  std::vector<std::vector<double> > filter = create_Gaussian_filter(5, 1);

  // creating two heat maps.
  Persistence_heat_maps<constant_scaling_function> hm1(persistence1, filter, false, 20, 0, 11);
  Persistence_heat_maps<constant_scaling_function> hm2(persistence2, filter, false, 20, 0, 11);

  std::vector<Persistence_heat_maps<constant_scaling_function>*> vector_of_maps;
  vector_of_maps.push_back(&hm1);
  vector_of_maps.push_back(&hm2);

  // compute median/mean of a vector of heat maps:
  Persistence_heat_maps<constant_scaling_function> mean;
  mean.compute_mean(vector_of_maps);
  Persistence_heat_maps<constant_scaling_function> median;
  median.compute_median(vector_of_maps);

  // to compute L^1 distance between hm1 and hm2:
  std::cout << "The L^1 distance is : " << hm1.distance(hm2, 1) << std::endl;

  // to average of hm1 and hm2:
  std::vector<Persistence_heat_maps<constant_scaling_function>*> to_average;
  to_average.push_back(&hm1);
  to_average.push_back(&hm2);
  Persistence_heat_maps<constant_scaling_function> av;
  av.compute_average(to_average);

  // to compute scalar product of hm1 and hm2:
  std::cout << "Scalar product is : " << hm1.compute_scalar_product(hm2) << std::endl;

  return 0;
}
