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

#include <gudhi/Simplex_tree.h>


#include <algorithm>  // for std::prev_permutation
#include <vector>
#include <array>
#include <random>
#include <numeric>  // for std::iota
#include <cstddef>  // for std::size_t
#ifdef DEBUG_TRACES
#include <iostream>
#endif  // DEBUG_TRACES

std::random_device rd;

namespace Gudhi {

template <typename Floating_type>
std::vector<Floating_type> N_random_values(std::size_t nb_values, Floating_type min = 0., Floating_type max = 1.) {
  std::uniform_real_distribution<Floating_type> unif(0., 1.);
  std::mt19937 rand_engine(rd());
  std::vector<Floating_type> random_filtrations(nb_values);
  std::generate(random_filtrations.begin(), random_filtrations.end(),
                [&unif, &rand_engine]() { return unif(rand_engine); });
  return random_filtrations;
}

template <typename Vertex_handle>
std::vector<std::array<Vertex_handle, 2>> random_edges(Vertex_handle nb_vertices, double density = 0.15) {
  std::vector<std::array<Vertex_handle, 2>> permutations;
  if (nb_vertices < 2)
    return permutations;

  std::vector<bool> to_permute(nb_vertices);
  std::fill(to_permute.begin(), to_permute.begin() + 2, true);
  
  std::size_t nb_permutations = (nb_vertices * (nb_vertices - 1)) / 2;
  auto random_values = N_random_values<double>(nb_permutations);
  std::size_t idx = 0;
  do {
    // Keep only X% of the possible edges
    if (random_values[idx] > density) {
      idx++;
      continue;
    }
    idx++;

    std::array<Vertex_handle, 2> permutation = {-1, -1};
    for (Vertex_handle idx = 0; idx < nb_vertices; ++idx) {
      if (to_permute[idx]) {
        if (permutation[0] == -1) {
          permutation[0] = idx;
        } else {
          permutation[1] = idx;
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

template<typename Simplex_tree>
void simplex_tree_random_flag_complex(
    Simplex_tree& st,
    typename Simplex_tree::Vertex_handle nb_vertices,
    double density = 0.15) {
  using Vertex_handle = typename Simplex_tree::Vertex_handle;
  std::vector<Vertex_handle> vertices(nb_vertices);
  std::iota(vertices.begin(), vertices.end(), 0); // vertices is { 0, 1, 2, ..., 99 } when nb_vertices is 100
  st.insert_batch_vertices(vertices);

  auto edges = random_edges(nb_vertices, density);

  auto random_filtrations = N_random_values<typename Simplex_tree::Filtration_value>(edges.size());

  std::size_t idx = 0;
  for (auto edge : edges) {
    st.insert_simplex(edge, random_filtrations[idx]);
    idx++;
  }

}

}  // namespace Gudhi

#endif  // RANDOM_GRAPH_GENERATORS_H_
