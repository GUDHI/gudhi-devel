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
#include <gudhi/Simplex_tree/Simplex_tree_multi.h>

#include <iostream>
#include <initializer_list>

struct ST_MULTI {
public:
	typedef Gudhi::linear_indexing_tag Indexing_tag;
	typedef int Vertex_handle;
	typedef float value_type;
	using Filtration_value = Gudhi::multiparameter::multi_filtrations::Finitely_critical_multi_filtration<value_type>;
	typedef std::uint32_t Simplex_key;
	static const bool store_key = true;
	static const bool store_filtration = true;
	static const bool contiguous_vertices = false;
	static const bool link_nodes_by_label = true;
	static const bool stable_simplex_handles = false;
	static const bool is_multi_parameter = true;
};

using ST = Gudhi::Simplex_tree<ST_MULTI>;


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
  st.insert_simplex_and_subfaces(triangle012, {1,2,3}); // {1,2,3} can be any array-like vector-like
  st.insert_simplex_and_subfaces(edge03, {4,5,6});

  auto edge02 = {0, 2};
  ST::Simplex_handle e = st.find(edge02);
  // Finitely_critical_multi_filtration has an operator<<
  std::cout << st.filtration(e) << std::endl;
  assert(st.filtration(st.find(edge03)) == std::vector<float>({4,5,6}));
  
}
