/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Marc Glisse
 *
 *    Copyright (C) 2015  Inria
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
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
