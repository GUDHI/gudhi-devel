/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/Persistence_vectors.h>

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <limits>
#include <utility>

using Vector_distances_in_diagram =
    Gudhi::Persistence_representations::Vector_distances_in_diagram<Gudhi::Euclidean_distance>;

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

  // create two persistence vectors based on persistence1 and persistence2:
  Vector_distances_in_diagram v1(persistence1, std::numeric_limits<size_t>::max());
  Vector_distances_in_diagram v2(persistence2, std::numeric_limits<size_t>::max());

  // writing to a stream:
  std::cout << "v1 : " << v1 << std::endl;
  std::cout << "v2 : " << v2 << std::endl;

  // averages:
  Vector_distances_in_diagram average;
  average.compute_average({&v1, &v2});
  std::cout << "Average : " << average << std::endl;

  // computations of distances:
  std::cout << "l^1 distance : " << v1.distance(v2) << std::endl;

  // computations of scalar product:
  std::cout << "Scalar product of l1 and l2 : " << v1.compute_scalar_product(v2) << std::endl;

  // create a file with a gnuplot script:
  v1.plot("plot_of_vector_representation");

  return 0;
}
