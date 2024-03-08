/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2024 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef RANDOM_GRAPH_GENERATORS_H_
#define RANDOM_GRAPH_GENERATORS_H_

#include <algorithm>  // for std::prev_permutation
#include <vector>
#include <utility>  // for std::pair
#include <random>
#ifdef DEBUG_TRACES
#include <iostream>
#endif  // DEBUG_TRACES

std::random_device rd;

namespace Gudhi {

template <typename Vertex_handle>
std::vector<std::pair<Vertex_handle, Vertex_handle>> random_edges(Vertex_handle nb_vertices, double density = 0.15) {
  std::vector<std::pair<Vertex_handle, Vertex_handle>> permutations;
  if (nb_vertices < 2)
    return permutations;

  std::uniform_real_distribution<double> unif(0., 1.);
  std::mt19937 rand_engine(rd());

  std::vector<bool> to_permute(nb_vertices);
  std::fill(to_permute.begin(), to_permute.begin() + 2, true);
  
  do {
    // Keep only X% of the possible edges
    if (unif(rand_engine) > density)
      continue;

    std::pair<Vertex_handle, Vertex_handle> permutation = {-1, -1};
    for (Vertex_handle idx = 0; idx < nb_vertices; ++idx) {
      if (to_permute[idx]) {
        if (permutation.first == -1) {
          permutation.first = idx;
        } else {
          permutation.second = idx;
          // No need to go further, only 2 'true' values
          break;
        }
      }
    }
#ifdef DEBUG_TRACES
    std::cout << permutation.first << ", " << permutation.second << std::endl;
#endif  // DEBUG_TRACES
    permutations.push_back(permutation);
    //std::cout << permutation.first << ", " << permutation.second << "\n";
  } while (std::prev_permutation(to_permute.begin(), to_permute.end()));

  return permutations;
}

}  // namespace Gudhi

#endif  // RANDOM_GRAPH_GENERATORS_H_
