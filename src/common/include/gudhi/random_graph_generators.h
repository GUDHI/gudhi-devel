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

#include <gudhi/random.h>


#include <algorithm>  // for std::shuffle
#include <vector>
#include <array>
#include <cmath>  // for std::round
#include <cstddef>  // for std::size_t

#include <boost/range/irange.hpp>

namespace Gudhi {

namespace random {

/* @brief Returns exactly \f$ density * nb_vertices(nb_vertices - 1) / 2 \f$ edges, chosen randomly.
 * 
 * @param[in] nb_vertices The number of vertices.
 * @param[in] density The percentage of edges to keep.
 * @return Some random edges.
 */
template <typename Vertex_handle>
std::vector<std::array<Vertex_handle, 2>> random_edges(Vertex_handle nb_vertices, double density = 0.15) {
  std::vector<std::array<Vertex_handle, 2>> edges;
  if (nb_vertices < 2)
    return edges;

  std::size_t nb_permutations = (nb_vertices * (nb_vertices - 1)) / 2;
  edges.reserve(nb_permutations);
  
  for (Vertex_handle u = 0; u < nb_vertices; u++) {
    for (Vertex_handle v = u + 1; v < nb_vertices; v++) {
      edges.push_back({u, v});
    }
  }
  
  std::shuffle(edges.begin(), edges.end(), Gudhi::random::get_default_random());
  edges.resize(std::round(nb_permutations * density));
  return edges;
}

/* @brief Constructs the Simplex_tree with nb_vertices and @ref Gudhi::random::random_edges. Each edge is initialized
 * with a random filtration value in the interval [filt_min, filt_max].
 * 
 * @param[in] st The Simplex_tree to initialize.
 * @param[in] nb_vertices The number of vertices.
 * @param[in] density The percentage of edges to keep.
 * @param[in] filt_min The minimal filtration value for the random filtration value of the edges. Default value is `0.`.
 * @param[in] filt_max The maximal filtration value for the random filtration value of the edges. Default value is `1.`.
 */
template<typename Simplex_tree>
void simplex_tree_random_graph(
    Simplex_tree& st,
    typename Simplex_tree::Vertex_handle nb_vertices,
    double density = 0.15,
    typename Simplex_tree::Filtration_value filt_min = 0.,
    typename Simplex_tree::Filtration_value filt_max = 1.) {
  st.insert_batch_vertices(boost::irange(nb_vertices));

  auto edges = random_edges(nb_vertices, density);

  using Filtration_value = typename Simplex_tree::Filtration_value;
  std::vector<Filtration_value> random_filtrations;
  // Shortcut when filt_min == filt_max
  if (filt_min == filt_max)
    random_filtrations.resize(edges.size(), filt_min);
  else
    random_filtrations = Gudhi::random::get_uniform_range<Filtration_value>(edges.size(), filt_min, filt_max);

  std::size_t idx = 0;
  for (auto edge : edges) {
    st.insert_simplex(edge, random_filtrations[idx]);
    idx++;
  }

}

}  // namespace random

}  // namespace Gudhi

#endif  // RANDOM_GRAPH_GENERATORS_H_
