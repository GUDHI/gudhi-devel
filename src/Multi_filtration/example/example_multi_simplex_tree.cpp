/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       David Loiseaux
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <initializer_list>

#include <gudhi/Simplex_tree.h>
#include <gudhi/Multi_parameter_filtration.h>

// Options for the simplex tree
template <typename MultiFiltrationValue>
struct Simplex_tree_options_multidimensional_filtration : Gudhi::Simplex_tree_options_default {
  using Filtration_value = MultiFiltrationValue;
};

using Multi_filtration_value = Gudhi::multi_filtration::Multi_parameter_filtration<double>;
using ST = Gudhi::Simplex_tree<Simplex_tree_options_multidimensional_filtration<Multi_filtration_value> >;

int main()
{
  ST st;

  /* Complex to build. */
  /*    1              */
  /*    o              */
  /*   /X\             */
  /*  o---o---o        */
  /*  2   0   3        */

  auto triangle012 = {0, 1, 2};
  auto edge03 = {0, 3};
  // Inserts the triangle with multi filtration value {1,2,3}
  st.insert_simplex_and_subfaces(triangle012, Multi_filtration_value({1, 2, 3}));
  // The filtration value can also be given as an initializer list
  st.insert_simplex_and_subfaces(edge03, std::initializer_list<double>{4, 5, 6});

  auto edge02 = {0, 2};
  ST::Simplex_handle e = st.find(edge02);
  // Multi_filtration_value has an operator<<
  std::cout << st.filtration(e) << std::endl;
  GUDHI_CHECK(st.filtration(st.find(edge03)) == ST::Filtration_value({4, 5, 6}),
              "edge03 does not have the right value.");
}
