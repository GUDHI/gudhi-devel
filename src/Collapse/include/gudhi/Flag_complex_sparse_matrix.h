/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Siddharth Pritam
 *
 *    Copyright (C) 2020 Inria
 *
 *    Modification(s):
 *      - 2020/03 Vincent Rouvreau: integration to the gudhi library
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef FLAG_COMPLEX_SPARSE_MATRIX_H_
#define FLAG_COMPLEX_SPARSE_MATRIX_H_

#include <gudhi/graph_simplicial_complex.h>
#include <gudhi/Debug_utils.h>

#include <boost/functional/hash.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <Eigen/Sparse>

#ifdef GUDHI_USE_TBB
#include <tbb/parallel_sort.h>
#endif

#include <iostream>
#include <utility>  // for std::pair
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <tuple>  // for std::tie
#include <algorithm>  // for std::includes
#include <iterator>  // for std::inserter

namespace Gudhi {

namespace collapse {

/**
 * \class Flag_complex_sparse_matrix
 * \brief Flag complex sparse matrix data structure.
 *
 * \ingroup collapse
 *
 * \details
 * This class stores a <a target="_blank" href="https://en.wikipedia.org/wiki/Clique_complex">Flag complex</a>
 * in an <a target="_blank" href="https://eigen.tuxfamily.org/dox/group__TutorialSparse.html">Eigen sparse matrix</a>.
 *
 * \tparam Vertex type must be a signed integer type. It admits a total order <.
 * \tparam Filtration type for the value of the filtration function. Must be comparable with <.
 */
template<typename Vertex, typename Filtration>
class Flag_complex_sparse_matrix {
 public:
  /** \brief Re-define Vertex as Vertex_handle type to ease the interface with compute_proximity_graph. */
  using Vertex_handle = Vertex;
  /** \brief Re-define Filtration as Filtration_value type to ease the interface with compute_proximity_graph. */
  using Filtration_value = Filtration;
  /** \brief This is an ordered pair, An edge is stored with convention of the first element being the smaller i.e
   * {2,3} not {3,2}. However this is at the level of row indices on actual vertex lables.
   */
  using Edge = std::pair<Vertex_handle, Vertex_handle>;

 private:
  // internal numbering of vertices and edges
  using IVertex = std::size_t;
  using Edge_index = std::size_t;

  // The sparse matrix data type
  // (Eigen::SparseMatrix<Edge_index, Eigen::RowMajor> has slow insertions)
  using Sparse_vector = Eigen::SparseVector<Edge_index>;
  using Sparse_row_matrix = std::vector<Sparse_vector>;

  // A range of row indices
  using IVertex_vector = std::vector<IVertex>;

 public:
  /** \brief Filtered_edge is a type to store an edge with its filtration value. */
  using Filtered_edge = std::pair<Edge, Filtration_value>;
  /** \brief Proximity_graph is a type that can be used to construct easily a Flag_complex_sparse_matrix. */
  using Proximity_graph = Gudhi::Proximity_graph<Flag_complex_sparse_matrix>;

 private:
  // Map from row index to its vertex handle
  std::vector<Vertex_handle> row_to_vertex_;

  // Index of the current edge in the backwards walk. Edges <= current_backward are part of the temporary graph,
  // while edges > current_backward are removed unless critical_edge_indicator_.
  Edge_index current_backward = -1;

  // Map from edge to its index
  std::unordered_map<Edge, Edge_index, boost::hash<Edge>> edge_to_index_map_;

  // Boolean vector to indicate if the edge is critical.
  std::vector<bool> critical_edge_indicator_;

  // Map from vertex handle to its row index
  std::unordered_map<Vertex_handle, IVertex> vertex_to_row_;

  // Stores the Sparse matrix of Filtration values representing the original graph.
  // The matrix rows and columns are indexed by IVertex.
  Sparse_row_matrix sparse_row_adjacency_matrix_;

  // The input, a vector of filtered edges.
  std::vector<Filtered_edge> f_edge_vector_;

  // Edge e is the actual edge (u,v), with Vertex_handle u and v, not IVertex.
  bool edge_is_dominated(const Edge& edge) const
  {
    Vertex_handle u = std::get<0>(edge);
    Vertex_handle v = std::get<1>(edge);

    const IVertex rw_u = vertex_to_row_.at(u);
    const IVertex rw_v = vertex_to_row_.at(v);
#ifdef DEBUG_TRACES
    std::cout << "The edge {" << u << ", " << v <<  "} is going for domination check." << std::endl;
#endif  // DEBUG_TRACES
    auto common_neighbours = open_common_neighbours_row_index(rw_u, rw_v);
#ifdef DEBUG_TRACES
    std::cout << "And its common neighbours are." << std::endl;
    for (auto neighbour : common_neighbours) {
      std::cout << row_to_vertex_[neighbour] << ", " ;
    }
    std::cout<< std::endl;
#endif  // DEBUG_TRACES
    if (common_neighbours.size() == 1)
      return true;
    else
      for (auto rw_c : common_neighbours) {
        auto neighbours_c = neighbours_row_index(rw_c, true);
        // If neighbours_c contains the common neighbours.
        if (std::includes(neighbours_c.begin(), neighbours_c.end(),
                          common_neighbours.begin(), common_neighbours.end()))
          return true;
      }
    return false;
  }

  // Returns the edges connecting u and v (extremities of crit) to their common neighbors (not themselves)
  std::set<Edge_index> three_clique_indices(Edge_index crit) {
    std::set<Edge_index> edge_indices;

    Edge edge = std::get<0>(f_edge_vector_[crit]);
    Vertex_handle u = std::get<0>(edge);
    Vertex_handle v = std::get<1>(edge);

#ifdef DEBUG_TRACES
    std::cout << "The  current critical edge to re-check criticality with filt value is : f {" << u << "," << v
              << "} = " << std::get<1>(f_edge_vector_[crit]) << std::endl;
#endif  // DEBUG_TRACES
    auto rw_u = vertex_to_row_[u];
    auto rw_v = vertex_to_row_[v];

    IVertex_vector common_neighbours = open_common_neighbours_row_index(rw_u, rw_v);

    for (auto rw_c : common_neighbours) {
      auto e_with_new_nbhr_v = std::minmax(u, row_to_vertex_[rw_c]);
      auto e_with_new_nbhr_u = std::minmax(v, row_to_vertex_[rw_c]);
      edge_indices.emplace(edge_to_index_map_[e_with_new_nbhr_v]);
      edge_indices.emplace(edge_to_index_map_[e_with_new_nbhr_u]);
    }
    return edge_indices;
  }

  // Detect and set all edges that are becoming critical
  template<typename FilteredEdgeOutput>
  void set_edge_critical(Edge_index indx, Filtration_value filt, FilteredEdgeOutput filtered_edge_output) {
#ifdef DEBUG_TRACES
    std::cout << "The curent index  with filtration value " << indx << ", " << filt << " is primary critical" <<
    std::endl;
#endif  // DEBUG_TRACES
    std::set<Edge_index> effected_indices = three_clique_indices(indx);
    // Cannot use boost::adaptors::reverse in such dynamic cases apparently
    for (auto it = effected_indices.rbegin(); it != effected_indices.rend(); ++it) {
      current_backward = *it;
      Edge edge = std::get<0>(f_edge_vector_[current_backward]);
      Vertex_handle u = std::get<0>(edge);
      Vertex_handle v = std::get<1>(edge);
      // If current_backward is not critical so it should be processed, otherwise it stays in the graph
      if (!critical_edge_indicator_[current_backward]) {
        if (!edge_is_dominated(edge)) {
#ifdef DEBUG_TRACES
          std::cout << "The curent index became critical " << current_backward  << std::endl;
#endif  // DEBUG_TRACES
          critical_edge_indicator_[current_backward] = true;
          filtered_edge_output({u, v}, filt);
          std::set<Edge_index> inner_effected_indcs = three_clique_indices(current_backward);
          for (auto inr_idx : inner_effected_indcs) {
            if(inr_idx < current_backward) // && !critical_edge_indicator_[inr_idx]
              effected_indices.emplace(inr_idx);
          }
#ifdef DEBUG_TRACES
          std::cout << "The following edge is critical with filt value: {" << u << "," << v << "}; "
            << filt << std::endl;
#endif  // DEBUG_TRACES
        }
      }
    }
    // Clear the implicit "removed from graph" data structure
    current_backward = -1;
  }

  // Returns list of neighbors of a particular vertex.
  IVertex_vector neighbours_row_index(IVertex rw_u, bool closed) const
  {
    IVertex_vector neighbors;
    neighbors.reserve(sparse_row_adjacency_matrix_[rw_u].nonZeros()); // too much, but who cares
#ifdef DEBUG_TRACES
    std::cout << "The neighbours of the vertex: " << row_to_vertex_[rw_u] << " are. " << std::endl;
#endif  // DEBUG_TRACES
    // Iterate over the neighbors
    for (typename Sparse_vector::InnerIterator it(sparse_row_adjacency_matrix_[rw_u]); it; ++it) {
      IVertex rw_v = it.index();
      if (!closed && rw_u == rw_v) continue;
      Edge_index ei;
      // If the vertex v is not dominated and the edge {u,v} is still in the matrix
      if ((closed && rw_u == rw_v) ||
          (ei = it.value()) <= current_backward ||
          critical_edge_indicator_[ei]) {
        neighbors.push_back(rw_v);
#ifdef DEBUG_TRACES
        std::cout << row_to_vertex_[rw_v] << ", " ;
#endif  // DEBUG_TRACES
      }
    }
#ifdef DEBUG_TRACES
    std::cout << std::endl;
#endif  // DEBUG_TRACES
    return neighbors;
  }

  // Returns the list of open neighbours of the edge :{u,v}.
  IVertex_vector open_common_neighbours_row_index(IVertex rw_u, IVertex rw_v) const
  {
    IVertex_vector non_zero_indices_u = neighbours_row_index(rw_u, false);
    IVertex_vector non_zero_indices_v = neighbours_row_index(rw_v, false);
    IVertex_vector common;
    common.reserve(std::min(non_zero_indices_u.size(), non_zero_indices_v.size()));
    std::set_intersection(non_zero_indices_u.begin(), non_zero_indices_u.end(), non_zero_indices_v.begin(),
                          non_zero_indices_v.end(), std::back_inserter(common));

    return common;
  }

  // Insert a vertex in the data structure
  IVertex insert_vertex(Vertex_handle vertex) {
    auto n = row_to_vertex_.size();
    auto result = vertex_to_row_.emplace(vertex, n);
    // If it was not already inserted - Value won't be updated by emplace if it is already present
    if (result.second) {
      // Expand the matrix. The size of rows is irrelevant.
      sparse_row_adjacency_matrix_.emplace_back((std::numeric_limits<Eigen::Index>::max)());
      // Initializing the diagonal element of the adjency matrix corresponding to rw_b.
      sparse_row_adjacency_matrix_[n].insert(n) = -1; // not an edge
      // Must be done after reading its size()
      row_to_vertex_.push_back(vertex);
    }
    return result.first->second;
  }

  // Insert an edge in the data structure
  // @exception std::invalid_argument In debug mode, if u == v
  void insert_new_edge(Vertex_handle u, Vertex_handle v, Edge_index idx)
  {
    GUDHI_CHECK((u != v), std::invalid_argument("Flag_complex_sparse_matrix::insert_new_edge with u == v"));
    // The edge must not be added before, it should be a new edge.
    IVertex rw_u = insert_vertex(u);
    IVertex rw_v = insert_vertex(v);
#ifdef DEBUG_TRACES
    std::cout << "Inserting the edge " << u <<", " << v << std::endl;
#endif  // DEBUG_TRACES
    sparse_row_adjacency_matrix_[rw_u].insert(rw_v) = idx;
    sparse_row_adjacency_matrix_[rw_v].insert(rw_u) = idx;
  }

 public:
  /** \brief Flag_complex_sparse_matrix constructor from a range of filtered edges.
   *
   * @param[in] filtered_edge_range Range of filtered edges. Filtered edges must be in
   * `Flag_complex_sparse_matrix::Filtered_edge`.
   *
   * There is no need the range to be sorted, as it will be performed in
   * `Flag_complex_sparse_matrix::filtered_edge_collapse`.
   */
  template<typename Filtered_edge_range>
  Flag_complex_sparse_matrix(const Filtered_edge_range& filtered_edge_range)
  : f_edge_vector_(filtered_edge_range.begin(), filtered_edge_range.end()) { }

  /** \brief Flag_complex_sparse_matrix constructor from a proximity graph, cf. `Gudhi::compute_proximity_graph`.
   *
   * @param[in] one_skeleton_graph The one skeleton graph. The graph must be in
   * `Flag_complex_sparse_matrix::Proximity_graph`.
   *
   * The constructor is computing and filling a vector of `Flag_complex_sparse_matrix::Filtered_edge`
   */
  Flag_complex_sparse_matrix(const Proximity_graph& one_skeleton_graph) {
    // Insert all edges
    for (auto edge_it = boost::edges(one_skeleton_graph);
         edge_it.first != edge_it.second; ++edge_it.first) {
      auto edge = *(edge_it.first);
      Vertex_handle u = source(edge, one_skeleton_graph);
      Vertex_handle v = target(edge, one_skeleton_graph);
      f_edge_vector_.push_back({{u, v}, boost::get(Gudhi::edge_filtration_t(), one_skeleton_graph, edge)});
    }
  }

  /** \brief Performs edge collapse in a increasing sequence of the filtration value.
   *
   * \tparam FilteredEdgeOutput is a functor that furnishes `({Vertex_handle u, Vertex_handle v}, Filtration_value f)`
   * that will get called on the output edges, in non-decreasing order of filtration.
   */
  template<typename FilteredEdgeOutput>
  void filtered_edge_collapse(FilteredEdgeOutput filtered_edge_output) {
    // Sort edges
    auto sort_by_filtration = [](const Filtered_edge& edge_a, const Filtered_edge& edge_b) -> bool
    {
      return (get<1>(edge_a) < get<1>(edge_b)); 
    };

#ifdef GUDHI_USE_TBB
    tbb::parallel_sort(f_edge_vector_.begin(), f_edge_vector_.end(), sort_by_filtration);
#else
    std::sort(f_edge_vector_.begin(), f_edge_vector_.end(), sort_by_filtration);
#endif

    for (Edge_index endIdx = 0; endIdx < f_edge_vector_.size(); endIdx++) {
      Filtered_edge fec = f_edge_vector_[endIdx];
      Edge edge = std::get<0>(fec);
      Vertex_handle u = std::get<0>(edge);
      Vertex_handle v = std::get<1>(edge);
      Filtration_value filt = std::get<1>(fec);

      // Inserts the edge in the sparse matrix to update the graph (G_i)
      insert_new_edge(u, v, endIdx);

      edge_to_index_map_.emplace(std::minmax(u, v), endIdx);
      critical_edge_indicator_.push_back(false);

      if (!edge_is_dominated(edge)) {
        critical_edge_indicator_[endIdx] = true;
        filtered_edge_output({u, v}, filt);
        if (endIdx > 1)
          set_edge_critical(endIdx, filt, filtered_edge_output);
      }
    }
  }

};

}  // namespace collapse

}  // namespace Gudhi

#endif  // FLAG_COMPLEX_SPARSE_MATRIX_H_
