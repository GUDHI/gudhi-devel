/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2026 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <iostream>  // std::cout
#include <vector>

#include <gudhi/Matrix.h>
#include <gudhi/persistence_matrix_options.h>

using Gudhi::persistence_matrix::Default_options;

struct RU_swap_options : Default_options<> {
  static const bool has_vine_update = true;      // for swap
  static const bool has_column_pairings = true;  // for barcode
};

using Matrix = Gudhi::persistence_matrix::Matrix<RU_swap_options>;
using Bar = typename Matrix::Bar;

void print_current_barcode(const Matrix& M, const std::vector<double>& filtrationValues, int step)
{
  std::cout << "Barcode " << step << ":\n\n";
  const auto& barcode = M.get_current_barcode();
  for (const auto& bar : barcode) {
    // bar with indices
    std::cout << bar << "\n";
    // bar with filtration values
    std::cout << "[" << bar.dim << "] ";
    std::cout << filtrationValues[bar.birth] << " - ";
    if (bar.death == Bar::inf)
      std::cout << "inf";
    else
      std::cout << filtrationValues[bar.death];
    std::cout << "\n\n";
  }
  std::cout << "\n";
}

int main()
{
  Matrix mp({{}, {}, {}, {}, {0, 1}, {0, 3}, {0, 2}, {1, 2}, {2, 3}, {4, 6, 7}});
  std::vector<double> filtrationValues = {0, 0.1, 0.2, 0.3, 1, 1.1, 2, 2.1, 3, 4};

  std::cout << "Initialized state:\n\n";
  print_current_barcode(mp, filtrationValues, 0);

  std::cout << "Filtration value of the fourth vertex changes to 1.2: ";
  std::cout << "swap index 3 and 4 as well as 4 and 5:\n\n";
  filtrationValues = {0, 0.1, 0.2, 1, 1.1, 1.2, 2, 2.1, 3, 4};
  mp.vine_swap(3);
  mp.vine_swap(4);
  print_current_barcode(mp, filtrationValues, 1);

  std::cout << "Filtration value of the triangle changes to 2.4: ";
  std::cout << "swap index 8 and 9:\n\n";
  filtrationValues = {0, 0.1, 0.2, 1, 1.1, 1.2, 2, 2.1, 2.4, 3};
  mp.vine_swap(8);
  print_current_barcode(mp, filtrationValues, 2);

  std::cout << "Filtration value of the fourth edge changes to 1.5: ";
  std::cout << "swap index 6 and 7:\n\n";
  filtrationValues = {0, 0.1, 0.2, 1, 1.1, 1.2, 1.5, 2, 2.4, 3};
  mp.vine_swap(6);
  print_current_barcode(mp, filtrationValues, 3);
}
