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

#include <gudhi/filtered_zigzag_persistence.h>

using Zigzag_persistence = Gudhi::zigzag_persistence::Filtered_zigzag_persistence<>;
using Dimension = Zigzag_persistence::Dimension;
using Filtration_value = Zigzag_persistence::Filtration_value;

int main() {
  std::clog << "********* Minimalistic example of usage of the Filtered_zigzag_persistence class ********" << std::endl;

  // Filtered_zigzag_persistence(callback) with for example callback method as a anonymous lambda
  Zigzag_persistence zp([](Dimension dim, Filtration_value birth, Filtration_value death) {
    std::cout << "[" << dim << "] " << birth << " - " << death << std::endl;
  });

  // It is important that the operations of insertions and removals are made **in the same order** as in the zigzag
  // filtration ones wants to compute the barcode from.
  // A face can be identified in the boundaries by any given numerical label, it is just important that the given
  // filtration values are monotonous (ie., either only increasing or only decreasing).

  // inserts vertex 2 at filtration value 0.1 -> birth at 0.1 of 0-cycle
  zp.insert_cell(2, {}, 0, 0.1);
  // inserts vertex 4 at filtration value 0.1 -> birth at 0.1 of 0-cycle
  zp.insert_cell(4, {}, 0, 0.1);
  // inserts edge 5 = (2,4) at filtration value 0.3 -> death at 0.3 -> outputs (0, 0.1, 0.3)
  zp.insert_cell(5, {2, 4}, 1, 0.3);
  // inserts vertex 3 at filtration value 0.4 -> birth at 0.4 of 0-cycle
  zp.insert_cell(3, {}, 0, 0.4);
  // inserts edge 6 = (2,3) at filtration value 0.4 -> death at 0.4 of the cycle born at 0.4 -> outputs nothing
  zp.insert_cell(6, {2, 3}, 1, 0.4);
  // inserts edge 9 = (3,4) at filtration value 1.2 -> birth at 1.2 of 1-cycle
  zp.insert_cell(9, {4, 3}, 1, 1.2);
  // removes edge 6 at filtration value 1.5 -> death at 1.5 -> outputs (1, 1.2, 1.5)
  zp.remove_cell(6, 1.5);
  // removes edge 5 at filtration value 2.0 -> birth at 2.0 of 0-cycle
  zp.remove_cell(5, 2.0);

  // Only the closed bars where output so far, so the open/infinite bars still need to be retrieved.

  // in this example, computes (0, 0.1) and (0, 2.0)
  zp.get_current_infinite_intervals([](Dimension dim, Filtration_value birth) {
    std::cout << "[" << dim << "] " << birth << " - inf" << std::endl;
  });

  return 0;
}
