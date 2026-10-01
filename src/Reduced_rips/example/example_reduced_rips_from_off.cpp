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
#include <gudhi/Points_off_io.h>

#include <iostream>
#include <string>
#include <vector>

int main(int argc, char** argv) {
  if (argc != 2 && argc != 3) {
    std::cerr << "Usage: " << argv[0] << " <input.off> [num_neighbors=0]\n";
    return 1;
  }
  unsigned int num_neighbors = argc == 3 ? static_cast<unsigned int>(std::stoul(argv[2])) : 0;

  Gudhi::Points_off_reader<std::vector<double>> off_reader(argv[1]);
  if (!off_reader.is_valid()) {
    std::cerr << "Could not read OFF file: " << argv[1] << '\n';
    return 1;
  }

  auto ph1 = Gudhi::reduced_rips::Reduced_rips<>::from_points(off_reader.get_point_cloud(), num_neighbors);
  std::cout << "Degree-1 barcode (birth death):\n";
  for (const auto& bar : ph1.persistence()) std::cout << bar[0] << ' ' << bar[1] << '\n';
  return 0;
}
