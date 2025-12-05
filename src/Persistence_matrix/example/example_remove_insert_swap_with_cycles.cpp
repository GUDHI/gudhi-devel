/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2024 Inria
 *
 *    Modification(s):
 *      - 2025/11 Jānis Lazovskis: Modified from example_representative_cycles_from_matrix
 */

#include <iostream>

#include <gudhi/Matrix.h>
#include <gudhi/persistence_matrix_options.h>

using Gudhi::persistence_matrix::Default_options;
using Gudhi::persistence_matrix::Column_types;

struct RU_rep_cycles_options : Default_options<Column_types::INTRUSIVE_LIST, true> 
{
  static const bool can_retrieve_representative_cycles = true;
  static const bool has_vine_update = true; // for insert, remove, swap
  static const bool has_removable_columns = true; // for remove
};

using RU_matrix = Gudhi::persistence_matrix::Matrix<RU_rep_cycles_options>;

template <class Matrix>
void print_representative_cycles_example()
{
  Matrix mp({ { },
              { },
              { },
              { },
              { 0 ,1 },
              { 0, 3 },
              { 0, 2 },
              { 1, 2 },
              { 2, 3 } 
            });

  std::cout << "Representative cycles at input:\n";
  auto rc = mp.get_representative_cycles();
  for (auto cycle : rc) {
    std::cout << mp.get_column_dimension(cycle[0]);
    std::cout << "-cycle: ";
    for (auto index : cycle) { std::cout << index << ", "; }
    std::cout << "\n";
  }

  std::cout << "Representative cycles after swapping 6 and 7:\n";
  mp.vine_swap(6);
  mp.update_representative_cycles();
  rc = mp.get_representative_cycles();
  for (auto cycle : rc) {
    std::cout << mp.get_column_dimension(cycle[0]);
    std::cout << "-cycle: ";
    for (auto index : cycle) { std::cout << index << ", "; }
    std::cout << "\n";
  }

  std::cout << "Representative cycles after inserting {1,3} at position 7:\n";
  mp.insert_maximal_cell( 7, {1,3}, 1);
  mp.update_representative_cycles();
  rc = mp.get_representative_cycles();
  for (auto cycle : rc) {
    std::cout << mp.get_column_dimension(cycle[0]);
    std::cout << "-cycle: ";
    for (auto index : cycle) { std::cout << index << ", "; }
    std::cout << "\n";
  }

  std::cout << "Representative cycles after swapping 8 and 9:\n";
  mp.vine_swap(8);
  mp.update_representative_cycles();
  rc = mp.get_representative_cycles();
  for (auto cycle : rc) {
    std::cout << mp.get_column_dimension(cycle[0]);
    std::cout << "-cycle: ";
    for (auto index : cycle) { std::cout << index << ", "; }
    std::cout << "\n";
  }

  std::cout << "Representative cycles after removing 5:\n";
  mp.remove_maximal_cell(5);
  mp.update_representative_cycles();
  rc = mp.get_representative_cycles();
  for (auto cycle : rc) {
    std::cout << mp.get_column_dimension(cycle[0]);
    std::cout << "-cycle: ";
    for (auto index : cycle) { std::cout << index << ", "; }
    std::cout << "\n";
  }

}

int main() {
  print_representative_cycles_example<RU_matrix>();
}
