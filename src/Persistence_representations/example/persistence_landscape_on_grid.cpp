/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/Persistence_landscape_on_grid.h>

#include <iostream>
#include <utility>
#include <vector>

using Persistence_landscape_on_grid = Gudhi::Persistence_representations::Persistence_landscape_on_grid;

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
  Persistence_landscape_on_grid l1(persistence1, 0, 11, 20);
  Persistence_landscape_on_grid l2(persistence2, 0, 11, 20);

  // This is how to compute integral of landscapes:
  std::clog << "Integral of the first landscape : " << l1.compute_integral_of_landscape() << std::endl;
  std::clog << "Integral of the second landscape : " << l2.compute_integral_of_landscape() << std::endl;

  // And here how to write landscapes to stream:
  std::clog << "l1 : " << l1 << std::endl;
  std::clog << "l2 : " << l2 << std::endl;

  // here are the maxima of the functions:
  std::clog << "Maximum of l1 : " << l1.compute_maximum() << std::endl;
  std::clog << "Maximum of l2 : " << l2.compute_maximum() << std::endl;

  // here are the norms of landscapes:
  std::clog << "L^1 Norm of l1 : " << l1.compute_norm_of_landscape(1.) << std::endl;
  std::clog << "L^1 Norm of l2 : " << l2.compute_norm_of_landscape(1.) << std::endl;

  // here is the average of landscapes:
  Persistence_landscape_on_grid average;
  average.compute_average({&l1, &l2});
  std::clog << "average : " << average << std::endl;

  // here is the distance of landscapes:
  std::clog << "Distance : " << l1.distance(l2) << std::endl;

  // here is the scalar product of landscapes:
  std::clog << "Scalar product : " << l1.compute_scalar_product(l2) << std::endl;

  // here is how to create a file which is suitable for visualization via gnuplot:
  average.plot("average_landscape");

  return 0;
}
