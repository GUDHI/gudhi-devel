/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Marc Glisse
 *
 *    Copyright (C) 2015  INRIA Saclay - Ile-de-France (France)
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

// FIXME: remove the first include
#include <gudhi/graph_simplicial_complex.h>
#include <gudhi/Simplex_tree.h>
#include <iostream>
#include <initializer_list>

using namespace Gudhi;

struct MyOptions : Simplex_tree_options_full_featured {
  // Not doing persistence, so we don't need those
  static constexpr bool store_key = false;
  static constexpr bool store_filtration = false;
  // I have few vertices
  typedef short Vertex_handle;
};
typedef Simplex_tree<MyOptions> ST;

// Dictionary should be private, but for now this is the easiest way.
static_assert(sizeof(ST::Dictionary::value_type) < sizeof(Simplex_tree<>::Dictionary::value_type),
    "Not storing the filtration and key should save some space");

int main() {
  ST st;

  /* Complex to build. */
  /*    1              */
  /*    o              */
  /*   /X\             */
  /*  o---o---o        */
  /*  2   0   3        */

  // FIXME: Replace std::vector<short> with auto
  std::vector<short> triangle012 = {0, 1, 2};
  std::vector<short> edge03 = {0, 3};
  st.insert_simplex_and_subfaces(triangle012);
  st.insert_simplex_and_subfaces(edge03);
  // FIXME: Remove this line
  st.set_dimension(2);

  std::vector<short> edge02 = {0, 2};
  ST::Simplex_handle e = st.find(edge02);
  assert(st.filtration(e) == 0); // We are not using filtrations so everything has value 0
  for(ST::Simplex_handle t : st.cofaces_simplex_range(e, 1)) // Only coface is 012
    {
      for(ST::Vertex_handle v : st.simplex_vertex_range(t)) // v in { 0, 1, 2 }
        std::cout << v;
      std::cout << '\n';
    }
}
