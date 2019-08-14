/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/Simplex_tree.h>

#include <iostream>

using Simplex_tree = Gudhi::Simplex_tree<>;
using Simplex_handle = Simplex_tree::Simplex_handle;

int main(int argc, char* const argv[]) {
  // Construct the Simplex Tree with a 1-skeleton graph example
  Simplex_tree stree;

  stree.insert_simplex({0, 1}, 0.);
  stree.insert_simplex({0, 2}, 1.);
  stree.insert_simplex({0, 3}, 2.);
  stree.insert_simplex({1, 2}, 3.);
  stree.insert_simplex({1, 3}, 4.);
  stree.insert_simplex({2, 3}, 5.);
  stree.insert_simplex({2, 4}, 6.);
  stree.insert_simplex({3, 6}, 7.);
  stree.insert_simplex({4, 5}, 8.);
  stree.insert_simplex({4, 6}, 9.);
  stree.insert_simplex({5, 6}, 10.);
  stree.insert_simplex({6}, 10.);

  stree.expansion_with_blockers(3, [&](Simplex_handle sh) {
    bool result = false;
    std::cout << "Blocker on [";
    // User can loop on the vertices from the given simplex_handle i.e.
    for (auto vertex : stree.simplex_vertex_range(sh)) {
      // We block the expansion, if the vertex '6' is in the given list of vertices
      if (vertex == 6) result = true;
      std::cout << vertex << ", ";
    }
    std::cout << "] ( " << stree.filtration(sh);
    // User can re-assign a new filtration value directly in the blocker (default is the maximal value of boudaries)
    stree.assign_filtration(sh, stree.filtration(sh) + 1.);

    std::cout << " + 1. ) = " << result << std::endl;

    return result;
  });

  std::cout << "********************************************************************\n";
  std::cout << "* The complex contains " << stree.num_simplices() << " simplices";
  std::cout << "   - dimension " << stree.dimension() << "\n";
  std::cout << "* Iterator on Simplices in the filtration, with [filtration value]:\n";
  for (auto f_simplex : stree.filtration_simplex_range()) {
    std::cout << "   "
              << "[" << stree.filtration(f_simplex) << "] ";
    for (auto vertex : stree.simplex_vertex_range(f_simplex)) std::cout << "(" << vertex << ")";
    std::cout << std::endl;
  }

  return 0;
}
