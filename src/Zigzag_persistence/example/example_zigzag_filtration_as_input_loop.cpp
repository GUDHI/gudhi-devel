/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <iostream>
#include <vector>

#include <gudhi/filtered_zigzag_persistence.h>

using ZP = Gudhi::zigzag_persistence::Filtered_zigzag_persistence_with_storage<>;
using Cell_handle = ZP::Cell_key;
using Filtration_value = ZP::Filtration_value;
using Interval_filtration = ZP::Filtration_value_interval;

void print_barcode(ZP& zp) {
  std::cout << std::endl << "Current barcode:" << std::endl;
  for (Interval_filtration& bar : zp.get_persistence_diagram(0, true)) {
    //stream out content of bar
    std::cout << bar << std::endl;
    //to access the content of the bar, it can either be used as a struct:
    //  bar.birth
    //  bar.death
    //  bar.dim
    //or as a tuple
    //  std::get<0>(bar) <- birth
    //  std::get<1>(bar) <- death
    //  std::get<2>(bar) <- dim
  }
}

void print_indices(ZP& zp) {
  std::cout << std::endl << "Current pairs:" << std::endl;
  for (auto& bar : zp.get_index_persistence_diagram()) {
    //stream out content of bar
    std::cout << bar << std::endl;
    //to access the content of the bar, it can either be used as a struct:
    //  bar.birth
    //  bar.death
    //  bar.dim
    //or as a tuple:
    //  std::get<0>(bar) <- birth
    //  std::get<1>(bar) <- death
    //  std::get<2>(bar) <- dim
  }
}

std::vector<std::vector<Cell_handle> > get_boundaries() {
  return {{},
          {},
          {},
          {0, 1},
          {0, 2},
          {},
          {1, 2},
          {},
          {5, 7},
          {},
          {3, 4, 6},
          {7, 9},
          {5, 9},
          {8, 11, 12},
          {10},                         // remove
          {13},                         // remove
          {1, 7},
          {3, 4, 6},
          {2, 7},
          {8, 11, 12},
          {0, 7},
          {4, 18, 20},
          {6, 16, 18},
          {3, 16, 20},
          {19},                         // remove
          {8},                          // remove
          {12},                         // remove
          {17, 21, 22, 23},
          {27}};                        // remove
}

std::vector<Filtration_value> get_filtration_values() {
  return {0, 0, 0, 
          1, 1, 1, 
          2, 2, 2, 
          3, 3, 3, 3, 
          4, 
          5, 
          6, 6, 6, 
          7, 7, 7, 7, 7, 7, 
          8, 
          9, 9, 9, 
          10};
}

std::vector<bool> get_directions() {
  return {true,  true, true, true, true, true, true, true, true, true,  true,  true,  true, true, 
          false, false, 
          true, true, true, true, true, true, true, true, 
          false, false, false, 
          true, 
          false};
}

std::vector<unsigned int> get_batch_sizes() {
  return {14, 2, 8, 3, 1, 1};
}

int main(int argc, char* const argv[]) {
  std::clog << "********** Example **********" << std::endl;
  
  ZP zp;

  std::vector<std::vector<Cell_handle> > simplices = get_boundaries();
  std::vector<Filtration_value> fils = get_filtration_values();
  std::vector<bool> dirs = get_directions();

  for (unsigned int i = 0; i < simplices.size(); ++i) {
    if (i > 0 && dirs[i] != dirs[i - 1]) {
      print_barcode(zp);
      print_indices(zp);
    }
    if (dirs[i]) {
      int dim = simplices[i].size() == 0 ? 0 : simplices[i].size() - 1;
      zp.insert_cell(i, simplices[i], dim, fils[i]);
    } else {
      auto id = simplices[i][0];
      zp.remove_cell(id, fils[i]);
    }
  }

  print_barcode(zp);
  print_indices(zp);

  return 0;
}
