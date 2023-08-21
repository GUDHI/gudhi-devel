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
#include <utility>
#include <vector>

#include <gudhi/Zigzag_persistence.h>
#include <gudhi/Simplex_tree.h>

using ST = Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_zigzag_persistence>;
using ZP = Gudhi::zigzag_persistence::Zigzag_persistence<ST>;
using Vertex_handle = ST::Vertex_handle;
using Filtration_value = ST::Filtration_value;
using interval_filtration = ZP::filtration_value_interval;

// void print_complex(ZP& zp) {
//   std::clog << std::endl << "Current complex:" << std::endl;
//   const auto& cpx = zp.get_complex();
//   for (const auto& sh : cpx.complex_simplex_range()) {
//     for (auto v : cpx.simplex_vertex_range(sh)) {
//       std::cout << v << " ";
//     }
//     std::cout << " - " << cpx.filtration(sh) << "" << std::endl;
//   }
// }

void print_barcode(ZP& zp) {
  std::clog << std::endl << "Current barcode:" << std::endl;
  for (auto& bar : zp.get_persistence_diagram(0, true)) {
    std::clog << std::floor(bar.birth()) << " - ";
    if (bar.death() == std::numeric_limits<Filtration_value>::infinity()) {
      std::clog << "inf";
    } else {
      std::clog << std::floor(bar.death());
    }
    std::clog << " (" << bar.dim() << ")" << std::endl;
  }
}

void print_indices(ZP& zp) {
  std::clog << std::endl << "Current pairs:" << std::endl;
  for (auto& bar : zp.get_index_persistence_diagram()) {
    std::clog << bar.birth() << " - ";
    std::clog << bar.death();
    std::clog << " (" << bar.dim() << ")" << std::endl;
  }
}

std::vector<std::vector<Vertex_handle> > get_simplices() {
  return {{0},
          {1},
          {2},
          {0, 1},
          {0, 2},
          {3},
          {1, 2},
          {4},
          {3, 4},
          {5},
          {0, 1, 2},
          {4, 5},
          {3, 5},
          {3, 4, 5},
          {0, 1, 2},            // remove
          {3, 4, 5},            // remove
          {1, 4},
          {0, 1, 2},
          {2, 4},
          {3, 4, 5},
          {0, 4},
          {0, 2, 4},
          {1, 2, 4},
          {0, 1, 4},
          {3, 4, 5},            // remove
          {3, 4},                   // remove
          {3, 5},                   // remove
          {0, 1, 2, 4},
          {0, 1, 2, 4}};    // remove
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

void one_by_one() {
  ZP zp;

  std::vector<std::vector<Vertex_handle> > simplices = get_simplices();
  std::vector<Filtration_value> fils = get_filtration_values();
  std::vector<bool> dirs = get_directions();

  for (unsigned int i = 0; i < simplices.size(); ++i) {
    if (i > 0 && dirs[i] != dirs[i - 1]) {
      // print_complex(zp);
      print_barcode(zp);
      print_indices(zp);
    }
    if (dirs[i]) {
      zp.insert_simplex(simplices[i], fils[i]);
    } else {
      zp.remove_simplex(simplices[i], fils[i]);
    }
  }
  // print_complex(zp);
  print_barcode(zp);
  print_indices(zp);
}

void in_batches() {
  ZP zp;

  std::vector<std::vector<Vertex_handle> > simplices = get_simplices();
  std::vector<Filtration_value> fils = get_filtration_values();
    std::vector<unsigned int> sizes = get_batch_sizes();

  unsigned int start;
  unsigned int end = 0;
  bool dir = true;  //first operation has to be an insertion
  for (auto s : sizes){
    start = end;
    end += s;
    if (dir){
        zp.insert_simplices_contiguously(simplices.begin() + start, 
                                       simplices.begin() + end, 
                                       fils.begin() + start);
    } else {
        zp.remove_simplices_contiguously(simplices.begin() + start, 
                                       simplices.begin() + end, 
                                       fils.begin() + start);
    }
    dir = !dir;
    // print_complex(zp);
    print_barcode(zp);
    print_indices(zp);
  }
}

int main(int argc, char* const argv[]) {
  std::clog << "********** Example one_by_one **********" << std::endl;
  one_by_one();

  std::clog << std::endl << "********** Example in_batches **********" << std::endl;
  in_batches();

  return 0;
}
