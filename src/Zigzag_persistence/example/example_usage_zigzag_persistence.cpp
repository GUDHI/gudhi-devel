/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2024 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <iostream>

#include <gudhi/zigzag_persistence.h>

using Zigzag_persistence = Gudhi::zigzag_persistence::Zigzag_persistence<>;
using Dimension = Zigzag_persistence::Dimension;
using Index = Zigzag_persistence::Index;

int main() {
  std::clog << "************* Minimalistic example of usage of the Zigzag_persistence class *************" << std::endl;

  // Zigzag_persistence(callback) with for example callback method as a anonymous lambda
  Zigzag_persistence zp([](Dimension dim, Index birth, Index death) {
    std::cout << "[" << dim << "] " << birth << " - " << death << std::endl;
  });

  // It is important that the operations of insertions and removals are made **in the same order** as in the zigzag
  // filtration one wants to compute the barcode from.
  // A cell has to be identified in the boundaries by the operation number the cell was inserted with in the sequence.

  // inserts vertex 0 -> birth at 0 of 0-cycle
  zp.insert_cell({}, 0);
  // inserts vertex 1 -> birth at 1 of 0-cycle
  zp.insert_cell({}, 0);
  // inserts edge 2 = (0,1) -> death at 2 -> outputs (0, 1, 2)
  zp.insert_cell({0, 1}, 1);
  // inserts vertex 3 -> birth at 3 of 0-cycle
  zp.insert_cell({}, 0);
  // inserts edge 4 = (0,3) -> death at 4 -> outputs (0, 3, 4)
  zp.insert_cell({0, 3}, 1);
  // inserts edge 5 = (1,3) -> birth at 5 of 1-cycle
  zp.insert_cell({1, 3}, 1);
  // removes edge 4 -> death at 6 -> outputs (1, 5, 6)
  zp.remove_cell(4);
  // removes edge 2 -> birth at 7 of 0-cycle
  zp.remove_cell(2);

  // Only the closed bars were output so far, so the open/infinite bars still need to be retrieved.

  // in this example, computes (0, 0) and (0, 7)
  zp.get_current_infinite_intervals(
      [](Dimension dim, Index birth) { std::cout << "[" << dim << "] " << birth << " - inf" << std::endl; });

  return 0;
}
