/*    This file is part of the Gudhi Library. The Gudhi library 
 *    (Geometric Understanding in Higher Dimensions) is a generic C++ 
 *    library for computational topology.
 *
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2014 Inria
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

#ifndef GRAPH_SIMPLICIAL_COMPLEX_H_
#define GRAPH_SIMPLICIAL_COMPLEX_H_

#include <boost/graph/adjacency_list.hpp>

#include <utility>  // for pair<>
#include <vector>
#include <map>
#include <tuple>  // for std::tie

namespace Gudhi {

/* Edge tag for Boost PropertyGraph. */
struct edge_filtration_t {
  typedef boost::edge_property_tag kind;
};

/* Vertex tag for Boost PropertyGraph. */
struct vertex_filtration_t {
  typedef boost::vertex_property_tag kind;
};

/** \brief Proximity_graph contains the vertices and edges with their filtration values in order to store the result
 * of `Gudhi::compute_proximity_graph` function.
 *
 * \tparam SimplicialComplexForProximityGraph furnishes `Filtration_value` type definition.
 *
 */
template <typename SimplicialComplexForProximityGraph>
using Proximity_graph = typename boost::adjacency_list < boost::vecS, boost::vecS, boost::directedS
, boost::property < vertex_filtration_t, typename SimplicialComplexForProximityGraph::Filtration_value >
, boost::property < edge_filtration_t, typename SimplicialComplexForProximityGraph::Filtration_value >>;

// setS returns pointer-like vertices (don't know why):
// source = 0x55ef4a232310 - target = 0x55ef4a232670 - filtration = 6.08276
// instead of :
// source = 0 - target = 1 - filtration = 6.08276
// template <typename SimplicialComplexForProximityGraph>
// using Proximity_graph = typename boost::adjacency_list < boost::setS, boost::setS, boost::undirectedS
//     , boost::property < vertex_filtration_t, typename SimplicialComplexForProximityGraph::Filtration_value >
//     , boost::property < edge_filtration_t, typename SimplicialComplexForProximityGraph::Filtration_value >>;

/** \brief Computes the proximity graph of the points.
 *
 * If points contains n elements, the proximity graph is the graph with n vertices, and an edge [u,v] iff the
 * distance function between points u and v is smaller than threshold.
 *
 * \tparam ForwardPointRange furnishes `.begin()` and `.end()` methods.
 *
 * \tparam Distance furnishes `operator()(const Point& p1, const Point& p2)`, where
 * `Point` is a point from the `ForwardPointRange`, and that returns a `Filtration_value`.
 */
template< typename SimplicialComplexForProximityGraph
          , typename ForwardPointRange
          , typename Distance >
Proximity_graph<SimplicialComplexForProximityGraph> compute_proximity_graph(
    const ForwardPointRange& points,
    typename SimplicialComplexForProximityGraph::Filtration_value threshold,
    Distance distance) {
  using Vertex_handle = typename SimplicialComplexForProximityGraph::Vertex_handle;
  using Filtration_value = typename SimplicialComplexForProximityGraph::Filtration_value;

  std::vector<std::pair< Vertex_handle, Vertex_handle >> edges;
  std::vector< Filtration_value > edges_fil;
  std::map< Vertex_handle, Filtration_value > vertices;

  Vertex_handle idx_u, idx_v;
  Filtration_value fil;
  idx_u = 0;
  for (auto it_u = points.begin(); it_u != points.end(); ++it_u) {
    idx_v = idx_u + 1;
    for (auto it_v = it_u + 1; it_v != points.end(); ++it_v, ++idx_v) {
      fil = distance(*it_u, *it_v);
      if (fil <= threshold) {
        auto edges_iter = edges.begin();
        auto edges_fil_iter = edges_fil.begin();
        while (edges_iter < edges.end() && edges_fil_iter != edges_fil.end() && fil > *edges_fil_iter) {
          ++edges_iter;
          ++edges_fil_iter;
        }

        edges.insert(edges_iter, std::make_pair(idx_u, idx_v));
        edges_fil.insert(edges_fil_iter, fil);
      }
    }
    ++idx_u;
  }

  // Points are labeled from 0 to idx_u-1
  Proximity_graph<SimplicialComplexForProximityGraph> skel_graph(edges.begin(), edges.end(), edges_fil.begin(), idx_u);

  auto vertex_prop = boost::get(vertex_filtration_t(), skel_graph);

  typename boost::graph_traits<Proximity_graph<SimplicialComplexForProximityGraph>>::vertex_iterator vi, vi_end;
  for (std::tie(vi, vi_end) = boost::vertices(skel_graph);
       vi != vi_end; ++vi) {
    boost::put(vertex_prop, *vi, 0.);
  }

  return skel_graph;
}

}  // namespace Gudhi

#endif  // GRAPH_SIMPLICIAL_COMPLEX_H_
