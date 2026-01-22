/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Jānis Lazovskis
 *
 *    Copyright (C) 2025 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <iostream>

#include <gudhi/Matrix.h>
#include <gudhi/persistence_matrix_options.h>

using Gudhi::persistence_matrix::Column_types;
using Gudhi::persistence_matrix::Default_options;

struct RU_rep_cycles_options : Default_options<Column_types::INTRUSIVE_LIST, true> {
  static const bool can_retrieve_representative_cycles = true;
  static const bool has_vine_update = true;        // for insert, remove, swap
  static const bool has_removable_columns = true;  // for remove
};

using RU_matrix = Gudhi::persistence_matrix::Matrix<RU_rep_cycles_options>;

void print_current_representative_cycles(const RU_matrix& M) {
  const auto& rc = M.get_all_representative_cycles();
  for (const auto& cycle : rc) {
    std::cout << M.get_column_dimension(cycle[0]);
    std::cout << "-cycle: ";
    for (auto index : cycle) {
      std::cout << index << ", ";
    }
    std::cout << "\n";
  }
}

template <class Matrix>
void remove_insert_swap_with_cycles() {
  Matrix mp({{}, {}, {}, {}, {0, 1}, {0, 3}, {0, 2}, {1, 2}, {2, 3}});

  std::cout << "Representative cycles at input:\n";
  print_current_representative_cycles(mp);

  std::cout << "Representative cycles after swapping 6 and 7:\n";
  mp.vine_swap(6);
  mp.update_all_representative_cycles();
  print_current_representative_cycles(mp);

  std::cout << "Representative cycles after inserting 0-cell at position 4:\n";
  mp.insert_maximal_cell(4, {});
  mp.update_all_representative_cycles();
  print_current_representative_cycles(mp);

  std::cout << "Representative cycles after swapping 8 and 9:\n";
  mp.vine_swap(8);
  mp.update_all_representative_cycles();
  print_current_representative_cycles(mp);

  std::cout << "Representative cycles after removing 5:\n";
  mp.remove_maximal_cell(5);
  mp.update_all_representative_cycles();
  print_current_representative_cycles(mp);
}

int main() { remove_insert_swap_with_cycles<RU_matrix>(); }
