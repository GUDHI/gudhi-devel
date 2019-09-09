/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Marc Glisse
 *
 *    Copyright (C) 2015 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/Simplex_tree.h>
#include <iostream>
#include <initializer_list>

struct MyOptions : Gudhi::Simplex_tree_options_full_featured {
  // Not doing persistence, so we don't need those
  static const bool store_key = false;
  static const bool store_filtration = false;
  // I have few vertices
  typedef short Vertex_handle;
};

using ST = Gudhi::Simplex_tree<MyOptions>;

// Dictionary should be private, but for now this is the easiest way.
static_assert(sizeof(ST::Dictionary::value_type) < sizeof(Gudhi::Simplex_tree<>::Dictionary::value_type),
    "Not storing the filtration and key should save some space");

int main() {
  ST st;

  /* Complex to build. */
  /*    1              */
  /*    o              */
  /*   /X\             */
  /*  o---o---o        */
  /*  2   0   3        */

  auto triangle012 = {0, 1, 2};
  auto edge03 = {0, 3};
  st.insert_simplex_and_subfaces(triangle012);
  st.insert_simplex_and_subfaces(edge03);

  auto edge02 = {0, 2};
  ST::Simplex_handle e = st.find(edge02);
  // We are not using filtrations so everything has value 0
  assert(st.filtration(e) == 0);
  for (ST::Simplex_handle t : st.cofaces_simplex_range(e, 1)) {
    // Only coface is 012
    for (ST::Vertex_handle v : st.simplex_vertex_range(t))  // v in { 0, 1, 2 }
        std::cout << v;
      std::cout << '\n';
  }
}
