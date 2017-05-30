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

#include <gudhi/Persistence_landscape.h>

using namespace Gudhi;
using namespace Gudhi::Persistence_representations;

#include <iostream>

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

  // create two persistence landscapes based on persistence1 and persistence2:
  Persistence_landscape l1(persistence1);
  Persistence_landscape l2(persistence2);

  // This is how to compute integral of landscapes:
  std::cout << "Integral of the first landscape : " << l1.compute_integral_of_landscape() << std::endl;
  std::cout << "Integral of the second landscape : " << l2.compute_integral_of_landscape() << std::endl;

  // And here how to write landscapes to stream:
  std::cout << "l1 : " << l1 << std::endl;
  std::cout << "l2 : " << l2 << std::endl;

  // Arithmetic operations on landscapes:
  Persistence_landscape sum = l1 + l2;
  std::cout << "sum : " << sum << std::endl;

  // here are the maxima of the functions:
  std::cout << "Maximum of l1 : " << l1.compute_maximum() << std::endl;
  std::cout << "Maximum of l2 : " << l2.compute_maximum() << std::endl;

  // here are the norms of landscapes:
  std::cout << "L^1 Norm of l1 : " << l1.compute_norm_of_landscape(1.) << std::endl;
  std::cout << "L^1 Norm of l2 : " << l2.compute_norm_of_landscape(1.) << std::endl;

  // here is the average of landscapes:
  Persistence_landscape average;
  average.compute_average({&l1, &l2});
  std::cout << "average : " << average << std::endl;

  // here is the distance of landscapes:
  std::cout << "Distance : " << l1.distance(l2) << std::endl;

  // here is the scalar product of landscapes:
  std::cout << "Scalar product : " << l1.compute_scalar_product(l2) << std::endl;

  // here is how to create a file which is suitable for visualization via gnuplot:
  average.plot("average_landscape");

  return 0;
}
