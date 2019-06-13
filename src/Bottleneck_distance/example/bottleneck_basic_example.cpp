/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Authors:       Francois Godi, small modifications by Pawel Dlotko
 *
 *    Copyright (C) 2015 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/Bottleneck.h>

#include <iostream>
#include <vector>
#include <utility>  // for pair
#include <limits>  // for numeric_limits

int main() {
  std::vector< std::pair<double, double> > v1, v2;

  v1.emplace_back(2.7, 3.7);
  v1.emplace_back(9.6, 14.);
  v1.emplace_back(34.2, 34.974);
  v1.emplace_back(3., std::numeric_limits<double>::infinity());

  v2.emplace_back(2.8, 4.45);
  v2.emplace_back(9.5, 14.1);
  v2.emplace_back(3.2, std::numeric_limits<double>::infinity());


  double b = Gudhi::persistence_diagram::bottleneck_distance(v1, v2);

  std::cout << "Bottleneck distance = " << b << std::endl;

  b = Gudhi::persistence_diagram::bottleneck_distance(v1, v2, 0.1);

  std::cout << "Approx bottleneck distance = " << b << std::endl;
}
