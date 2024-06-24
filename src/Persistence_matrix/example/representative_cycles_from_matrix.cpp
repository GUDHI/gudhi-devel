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

#include <gudhi/matrix.h>
#include <gudhi/persistence_matrix_options.h>

using Gudhi::persistence_matrix::Default_options;
using Gudhi::persistence_matrix::Column_types;

struct RU_rep_cycles_options : Default_options<Column_types::INTRUSIVE_LIST, true> 
{
  static const bool can_retrieve_representative_cycles = true;
};

struct Chain_rep_cycles_options : Default_options<Column_types::INTRUSIVE_LIST, true> 
{
  static const bool can_retrieve_representative_cycles = true;
  static const bool is_of_boundary_type = false;
};

using RU_matrix = Gudhi::persistence_matrix::Matrix<RU_rep_cycles_options>;
using Chain_matrix = Gudhi::persistence_matrix::Matrix<Chain_rep_cycles_options>;

template <class Matrix>
void print_representative_cycles_example()
{
  Matrix mp({ { },
              { },
              { },
              { },
              { },
              { },
              { },
              { 2, 3 },
              { 4, 5 },
              { 0, 2 },
              { 0, 1 },
              { 1, 3 },
              { 1, 2 },
              { 7, 11, 12 },
              { 9, 10, 12 },
              { 5, 6 },
              { 2, 4 },
              { 4, 6 },
              { 8, 15, 17 },
              { 3, 6 } 
            });

  auto rc = mp.get_representative_cycles();
  for (auto cycle : rc) {
    // cycle[0] gives the row index of a simplex in the cycle
    // because the simplices where indexed from 0 continously, the simplex represented by the row index cycle[0] is
    // the same simplex represented by the column at position cycle[0] in RU
    // that is why `mp.get_column_dimension(cycle[0])` gives us the dimension of the simplex for RU_matrix
    //
    // for the chain matrix, the row index will always represent a simplex ID. So,
    // `mp.get_column_dimension(mp.get_column_with_pivot(cycle[0]))` will always work to get the dimension
    // of the simplex. But in this particlar case, because of the simplex indexation and the fact that no swap
    // occured, mp.get_column_with_pivot(cycle[0]) == cycle[0] and so `mp.get_column_dimension(cycle[0])` also works.
    std::cout << mp.get_column_dimension(cycle[0]);
    std::cout << "-cycle: ";
    for (auto index : cycle) {
      std::cout << index << ", ";
    }
    std::cout << "\n";
  }
}

int main() {
  std::cout << "RU_matrix:\n";
  print_representative_cycles_example<RU_matrix>();
  std::cout << "\n";
  std::cout << "Chain_matrix:\n";
  print_representative_cycles_example<Chain_matrix>();
}
