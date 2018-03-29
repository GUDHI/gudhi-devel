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

#include <gudhi/Sliced_Wasserstein.h>

#include <iostream>
#include <vector>
#include <utility>

using SW = Gudhi::Persistence_representations::Sliced_Wasserstein;

int main(int argc, char** argv) {

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


  SW sw1(persistence1, 1, 100);
  SW sw2(persistence2, 1, 100);

  SW swex1(persistence1, 1, -1);
  SW swex2(persistence2, 1, -1);

  std::cout << "Approx SW distance: " << sw1.compute_sliced_wasserstein_distance(sw2) << std::endl;
  std::cout << "Exact SW distance: " << swex1.compute_sliced_wasserstein_distance(swex2) << std::endl;
  std::cout << "Approx SW kernel: " << sw1.compute_scalar_product(sw2) << std::endl;
  std::cout << "Exact  SW kernel: " << swex1.compute_scalar_product(swex2) << std::endl;
  std::cout << "Distance induced by approx SW kernel: " << sw1.distance(sw2) << std::endl;
  std::cout << "Distance induced by exact  SW kernel: " << swex1.distance(swex2) << std::endl;

  return 0;
}
