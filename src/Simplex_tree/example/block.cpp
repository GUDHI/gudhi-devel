/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2014
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

#include <gudhi/graph_simplicial_complex.h>
#include <gudhi/Simplex_tree.h>

#include <iostream>
#include <utility>  // for pair
#include <vector>

using Simplex_tree = Gudhi::Simplex_tree<>;
using Vertex_handle = Simplex_tree::Vertex_handle;
using Filtration_value = Simplex_tree::Filtration_value;
using typeVectorVertex = std::vector< Vertex_handle >;
using typePairSimplexBool = std::pair< Simplex_tree::Simplex_handle, bool >;

int main(int argc, char * const argv[]) {

  // Construct the Simplex Tree
  Simplex_tree simplexTree;

  simplexTree.insert_simplex({0, 1});
  simplexTree.insert_simplex({0, 2});
  simplexTree.insert_simplex({0, 3});
  simplexTree.insert_simplex({1, 2});
  simplexTree.insert_simplex({1, 3});
  simplexTree.insert_simplex({2, 3});
  simplexTree.insert_simplex({2, 4});
  simplexTree.insert_simplex({3, 6});
  simplexTree.insert_simplex({4, 5});
  simplexTree.insert_simplex({4, 6});
  simplexTree.insert_simplex({5, 6});
  simplexTree.insert_simplex({6});

  std::cout << "********************************************************************\n";
  // Display the Simplex_tree - Can not be done in the middle of 2 inserts
  std::cout << "* The complex contains " << simplexTree.num_simplices() << " simplices\n";
  std::cout << "   - dimension " << simplexTree.dimension() << "   - filtration " << simplexTree.filtration() << "\n";
  std::cout << "* Iterator on Simplices in the filtration, with [filtration value]:\n";
  for (auto f_simplex : simplexTree.filtration_simplex_range()) {
    std::cout << "   " << "[" << simplexTree.filtration(f_simplex) << "] ";
    for (auto vertex : simplexTree.simplex_vertex_range(f_simplex))
      std::cout << "(" << vertex << ")";
    std::cout << std::endl;
  }

  simplexTree.expansion_with_blockers(3, [](){return true;});

  simplexTree.initialize_filtration();
  std::cout << "********************************************************************\n";
  // Display the Simplex_tree - Can not be done in the middle of 2 inserts
  std::cout << "* The complex contains " << simplexTree.num_simplices() << " simplices\n";
  std::cout << "   - dimension " << simplexTree.dimension() << "   - filtration " << simplexTree.filtration() << "\n";
  std::cout << "* Iterator on Simplices in the filtration, with [filtration value]:\n";
  for (auto f_simplex : simplexTree.filtration_simplex_range()) {
    std::cout << "   " << "[" << simplexTree.filtration(f_simplex) << "] ";
    for (auto vertex : simplexTree.simplex_vertex_range(f_simplex))
      std::cout << "(" << vertex << ")";
    std::cout << std::endl;
  }

  return 0;
}
