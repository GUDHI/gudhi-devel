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

#include <cmath>
#include <iostream>
#include <vector>

int main() {
  // A hand-typed distance matrix in lower-triangular form: row i lists the distances from point i to points
  // 0..i-1 (the diagonal is taken as zero). These are the pairwise distances of a unit square unit sides.
  // The barcode here is the single loop (1, sqrt 2).
  const double square_diag = std::sqrt(2.0);
  std::vector<std::vector<double>> distances = {
      {},                      // point 0
      {1.0},                   // d(1,0)
      {square_diag, 1.0},      // d(2,0), d(2,1)
      {1.0, square_diag, 1.0}  // d(3,0), d(3,1), d(3,2)
  };

  auto rr = Gudhi::reduced_rips::Reduced_rips<>::from_distance_matrix(distances);

  std::cout << "Degree-1 barcode (birth death):\n";
  for (const auto& bar : rr.persistence()) std::cout << bar[0] << ' ' << bar[1] << '\n';
  return 0;
}
