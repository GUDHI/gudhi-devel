/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Mathieu Carriere
 *
 *    Copyright (C) 2018  INRIA (France)
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/Sliced_Wasserstein.h>

#include <iostream>
#include <vector>
#include <utility>

using Persistence_diagram = Gudhi::Persistence_representations::Persistence_diagram;
using SW = Gudhi::Persistence_representations::Sliced_Wasserstein;

int main(int argc, char** argv) {

  Persistence_diagram persistence1, persistence2;

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

  std::cout << "Approx SW kernel: " << sw1.compute_scalar_product(sw2) << std::endl;
  std::cout << "Exact  SW kernel: " << swex1.compute_scalar_product(swex2) << std::endl;
  std::cout << "Distance induced by approx SW kernel: " << sw1.distance(sw2) << std::endl;
  std::cout << "Distance induced by exact  SW kernel: " << swex1.distance(swex2) << std::endl;

  return 0;
}
