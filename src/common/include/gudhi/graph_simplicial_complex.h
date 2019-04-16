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
#include <tuple>      // for std::tie
#include <algorithm>  // for std::sort

namespace Gudhi {

/* Edge tag for Boost PropertyGraph. */
struct edge_filtration_t {
  typedef boost::edge_property_tag kind;
};

/* Vertex tag for Boost PropertyGraph. */
struct vertex_filtration_t {
  typedef boost::vertex_property_tag kind;
};

/** \brief Filtered_edges_vector contains the edges with their filtration values, sorted by ascending filtration values,
 * in order to store the result of `Gudhi::compute_edge_graph` function.
 *
 * \tparam SimplicialComplexForProximityGraph furnishes `Filtration_value` and `Vertex_handle` types definition.
 *
 */
template <typename SimplicialComplexForFilteredEdges>
class Filtered_edges_vector {
 public:
  /** \brief Type of the `Gudhi::Filtered_edges_vector` data structure.
   */
  using Filtered_edge = std::tuple<typename SimplicialComplexForFilteredEdges::Filtration_value,
                                   typename SimplicialComplexForFilteredEdges::Vertex_handle,
                                   typename SimplicialComplexForFilteredEdges::Vertex_handle>;
  using Filtered_edge_set = std::vector<Filtered_edge>;

  template <class EdgeIterator, class EdgePropertyIterator>
  Filtered_edges_vector(EdgeIterator first, EdgeIterator last, EdgePropertyIterator ep_iter,
                        typename SimplicialComplexForFilteredEdges::Vertex_handle n) {
    for (; first < last; ++first) {
      edges_.push_back(std::make_tuple(*ep_iter, first->first, first->second));
      ++ep_iter;
    }
    // By default the sort is done on the first element (Filtration_value) in ascending order.
    std::sort(edges_.begin(), edges_.end());
  }

  Filtered_edge_set edges() const { return edges_; }

  /** \brief Returns the edges vector at a specific index filtration value.
   */
  typename SimplicialComplexForFilteredEdges::Filtration_value get_filtration_at(std::size_t new_index) const {
    return std::get<0>(*(edges_.begin() + new_index));
  }

  /** \brief Returns the edges vector min filtration value.
   */
  typename SimplicialComplexForFilteredEdges::Filtration_value get_filtration_min() const {
    return std::get<0>(*(edges_.begin()));
  }

  /** \brief Returns the edges vector max filtration value.
   */
  typename SimplicialComplexForFilteredEdges::Filtration_value get_filtration_max() const {
    return std::get<0>(*(edges_.rbegin()));
  }

  /** \brief Returns the edges vector.
   */
  std::size_t size() const { return edges_.size(); }

  /** \brief Returns the sub-filtered edges vector from a `new_threshold` filtration value.
   */
  Filtered_edge_set sub_filter_edges_by_filtration(typename SimplicialComplexForFilteredEdges::Filtration_value new_threshold) const {
    auto edge_it = edges_.begin();
    while ((std::get<0>(*edge_it) < new_threshold) && (edge_it < edges_.end())) {
      ++edge_it;
    }
#ifdef DEBUG_TRACES
    if (edges_.begin() != edge_it) {
      Filtered_edge_set output(edges_.begin(), edge_it);
      Filtered_edge back = output.back();
      std::cout << "Filtered_edges_vector::sub_filter_edges_by_filtration threshold = "
                << std::get<0>(back) << " - size = " << output.size() << " - back = "
                << std::get<0>(back) << " - u = " << std::get<1>(back)
                << " - v = " << std::get<2>(back) << std::endl;
    } else {
      std::cout << "Filtered_edges_vector::sub_filter_edges_by_filtration is empty " << std::endl;
    }
#endif  // DEBUG_TRACES
    return Filtered_edge_set(edges_.begin(), edge_it);
  }

  /** \brief Returns the sub-filtered edges vector from a `new_index` index value.
   */
  Filtered_edge_set sub_filter_edges_by_index(std::size_t new_index) const {
    if (new_index >= size())
      return edges_;

#ifdef DEBUG_TRACES
    if (edges_.begin() != edges_.begin() + new_index + 1) {
      Filtered_edge_set output(edges_.begin(), edges_.begin() + new_index + 1);
      Filtered_edge back(output.back());
      std::cout << "Filtered_edges_vector::sub_filter_edges_by_filtration threshold = "
                << std::get<0>(back) << " - size = " << output.size() << " - back = "
                << std::get<0>(back) << " - u = " << std::get<1>(back)
                << " - v = " << std::get<2>(back) << std::endl;
    } else {
      std::cout << "Filtered_edges_vector::sub_filter_edges_by_index is empty " << std::endl;
    }
#endif  // DEBUG_TRACES
    return Filtered_edge_set(edges_.begin(), edges_.begin() + new_index + 1);
  }

private:
  Filtered_edge_set edges_;
};

/** \brief Computes the edge graph of the points.
 *
 * If points contains n elements, the edge graph is the graph with all edges [u,v] iff the
 * distance function between points u and v is smaller than threshold.
 *
 * \tparam ForwardPointRange furnishes `.begin()` and `.end()` methods.
 *
 * \tparam Distance furnishes `operator()(const Point& p1, const Point& p2)`, where
 * `Point` is a point from the `ForwardPointRange`, and that returns a `Filtration_value`.
 *
 * \tparam SimplicialComplexForEdgeGraph furnishes `Vertex_handle` and `Filtration_value` types declarations.
 */
template <template <class> class Edge_graph, class SimplicialComplexForEdgeGraph, typename ForwardPointRange,
          typename Distance>
Edge_graph<SimplicialComplexForEdgeGraph> compute_edge_graph(
    const ForwardPointRange& points, typename SimplicialComplexForEdgeGraph::Filtration_value threshold,
    Distance distance) {
  using Vertex_handle = typename SimplicialComplexForEdgeGraph::Vertex_handle;
  using Filtration_value = typename SimplicialComplexForEdgeGraph::Filtration_value;

  std::vector<std::pair<Vertex_handle, Vertex_handle>> edges;
  std::vector<Filtration_value> edges_fil;
  std::map<Vertex_handle, Filtration_value> vertices;

  Vertex_handle idx_u, idx_v;
  Filtration_value fil;
  idx_u = 0;
  for (auto it_u = points.begin(); it_u != points.end(); ++it_u) {
    idx_v = idx_u + 1;
    for (auto it_v = it_u + 1; it_v != points.end(); ++it_v, ++idx_v) {
      fil = distance(*it_u, *it_v);
      if (fil <= threshold) {
        edges.emplace_back(idx_u, idx_v);
        edges_fil.push_back(fil);
      }
    }
    ++idx_u;
  }

  // Points are labeled from 0 to idx_u-1
  Edge_graph<SimplicialComplexForEdgeGraph> graph_edges(edges.begin(), edges.end(), edges_fil.begin(), idx_u);
  return graph_edges;
}

/** \brief Proximity_graph contains the vertices and edges with their filtration values in order to store the result
 * of `Gudhi::compute_proximity_graph` function.
 *
 * \tparam SimplicialComplexForProximityGraph furnishes `Filtration_value` type definition.
 *
 */
template <typename SimplicialComplexForProximityGraph>
using Proximity_graph = typename boost::adjacency_list<
    boost::vecS, boost::vecS, boost::directedS,
    boost::property<vertex_filtration_t, typename SimplicialComplexForProximityGraph::Filtration_value>,
    boost::property<edge_filtration_t, typename SimplicialComplexForProximityGraph::Filtration_value>>;

/** \brief Computes the proximity graph of the points.
 *
 * If points contains n elements, the proximity graph is the graph with n vertices, and an edge [u,v] iff the
 * distance function between points u and v is smaller than threshold.
 *
 * \tparam ForwardPointRange furnishes `.begin()` and `.end()` methods.
 *
 * \tparam Distance furnishes `operator()(const Point& p1, const Point& p2)`, where
 * `Point` is a point from the `ForwardPointRange`, and that returns a `Filtration_value`.
 *
 * \tparam SimplicialComplexForProximityGraph furnishes `Vertex_handle` and `Filtration_value` types declarations.
 */
template <typename SimplicialComplexForProximityGraph, typename ForwardPointRange, typename Distance>
Proximity_graph<SimplicialComplexForProximityGraph> compute_proximity_graph(
    const ForwardPointRange& points, typename SimplicialComplexForProximityGraph::Filtration_value threshold,
    Distance distance) {
  // Points are labeled from 0 to idx_u-1
  Proximity_graph<SimplicialComplexForProximityGraph> skel_graph =
      compute_edge_graph<Proximity_graph, SimplicialComplexForProximityGraph>(points, threshold, distance);

  // Insert all vertices with a 0. filtration value
  auto vertex_prop = boost::get(vertex_filtration_t(), skel_graph);

  typename boost::graph_traits<Proximity_graph<SimplicialComplexForProximityGraph>>::vertex_iterator vi, vi_end;
  for (std::tie(vi, vi_end) = boost::vertices(skel_graph); vi != vi_end; ++vi) {
    boost::put(vertex_prop, *vi, 0.);
  }

  return skel_graph;
}

}  // namespace Gudhi

#endif  // GRAPH_SIMPLICIAL_COMPLEX_H_
