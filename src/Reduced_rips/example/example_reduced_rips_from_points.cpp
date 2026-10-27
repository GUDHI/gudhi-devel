/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Thomas Burnett
 *
 *    Copyright (C) 2026 Thomas Burnett, Musashi Koyama
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/Reduced_rips.h>

#include <iostream>
#include <vector>

int main() {
  // A hand-typed point cloud: the four corners of a unit square in the plane. Its degree-1 persistent
  // homology has a single loop, born when the four unit-length sides close the square and dying when a
  // diagonal (length sqrt 2) fills it in. The barcode is the single bar (1, sqrt 2).
  std::vector<std::vector<double>> points = {{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}};

  auto rr = Gudhi::reduced_rips::Reduced_rips<>::from_points(points);

  std::cout << "Degree-1 barcode (birth death):\n";
  for (const auto& bar : rr.persistence()) std::cout << bar[0] << ' ' << bar[1] << '\n';
  return 0;
}
