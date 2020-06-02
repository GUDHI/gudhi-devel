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
  // Row_index type in the sparse matrix
  using Row_index = std::size_t;

  // The sparse matrix data type
  using Sparse_row_matrix = Eigen::SparseMatrix<Filtration_value, Eigen::RowMajor>;

  // A range of row indices
  using Row_indices_vector = std::vector<Row_index>;

 public:
  /** \brief Filtered_edge is a type to store an edge with its filtration value. */
  using Filtered_edge = std::pair<Edge, Filtration_value>;
  /** \brief Proximity_graph is a type that can be used to construct easily a Flag_complex_sparse_matrix. */
  using Proximity_graph = Gudhi::Proximity_graph<Flag_complex_sparse_matrix>;

 private:
  // Map from row index to its vertex handle
  std::vector<Vertex_handle> row_to_vertex_;

  // Unordered set of removed edges. (to enforce removal from the matrix)
  std::unordered_set<Edge, boost::hash<Edge>> u_set_removed_edges_;

  // Unordered set of dominated edges. (to inforce removal from the matrix)
  std::unordered_set<Edge, boost::hash<Edge>> u_set_dominated_edges_;

  // Map from edge to its index
  std::unordered_map<Edge, Row_index, boost::hash<Edge>> edge_to_index_map_;

  // Boolean vector to indicate if the index is critical or not.
  std::vector<bool> critical_edge_indicator_;

  // Map from vertex handle to its row index
  std::unordered_map<Vertex_handle, Row_index> vertex_to_row_;

  // Stores the Sparse matrix of Filtration values representing the original graph.
  // This is row-major version of the same sparse-matrix, to facilitate easy access
  // to elements when traversing the matrix row-wise.
  Sparse_row_matrix sparse_row_adjacency_matrix_;

  // Vector of filtered edges, for edge-collapse, the indices of the edges are the row-indices.
  std::vector<Filtered_edge> f_edge_vector_;

  // Edge e is the actual edge (u,v). Not the row ids in the matrixs
  bool edge_is_dominated(const Edge& edge) const
  {
    Vertex_handle u = std::get<0>(edge);
    Vertex_handle v = std::get<1>(edge);

    const Row_index rw_u = vertex_to_row_.at(u);
    const Row_index rw_v = vertex_to_row_.at(v);
#ifdef DEBUG_TRACES
    std::cout << "The edge {" << u << ", " << v <<  "} is going for domination check." << std::endl;
#endif  // DEBUG_TRACES
    auto common_neighbours = closed_common_neighbours_row_index(rw_u, rw_v);
#ifdef DEBUG_TRACES
    std::cout << "And its common neighbours are." << std::endl;
    for (auto neighbour : common_neighbours) {
      std::cout << row_to_vertex_[neighbour] << ", " ;
    }
    std::cout<< std::endl;
#endif  // DEBUG_TRACES
    if (common_neighbours.size() > 2) {
      if (common_neighbours.size() == 3)
        return true;
      else
        for (auto rw_c : common_neighbours) {
          if (rw_c != rw_u && rw_c != rw_v) {
            auto neighbours_c = closed_neighbours_row_index(rw_c);
            // If neighbours_c contains the common neighbours.
            if (std::includes(neighbours_c.begin(), neighbours_c.end(), common_neighbours.begin(),
                              common_neighbours.end()))
              return true;
          }
        }
    }
    return false;
  }

  // Returns the edges connecting u and v (extremities of crit) to their common neighbors (not themselves)
  std::set<Row_index> three_clique_indices(Row_index crit) {
    std::set<Row_index> edge_indices;

    Edge edge = std::get<0>(f_edge_vector_[crit]);
    Vertex_handle u = std::get<0>(edge);
    Vertex_handle v = std::get<1>(edge);

#ifdef DEBUG_TRACES
    std::cout << "The  current critical edge to re-check criticality with filt value is : f {" << u << "," << v
              << "} = " << std::get<1>(f_edge_vector_[crit]) << std::endl;
#endif  // DEBUG_TRACES
    auto rw_u = vertex_to_row_[u];
    auto rw_v = vertex_to_row_[v];

    Row_indices_vector common_neighbours = closed_common_neighbours_row_index(rw_u, rw_v);

    if (common_neighbours.size() > 2) {
      for (auto rw_c : common_neighbours) {
        if (rw_c != rw_u && rw_c != rw_v) {
          auto e_with_new_nbhr_v = std::minmax(u, row_to_vertex_[rw_c]);
          auto e_with_new_nbhr_u = std::minmax(v, row_to_vertex_[rw_c]);
          edge_indices.emplace(edge_to_index_map_[e_with_new_nbhr_v]);
          edge_indices.emplace(edge_to_index_map_[e_with_new_nbhr_u]);
        }
      }
    }
    return edge_indices;
  }

  // Detect and set all indices that are becoming critical
  template<typename FilteredEdgeOutput>
  void set_edge_critical(Row_index indx, Filtration_value filt, FilteredEdgeOutput filtered_edge_output) {
#ifdef DEBUG_TRACES
    std::cout << "The curent index  with filtration value " << indx << ", " << filt << " is primary critical" <<
    std::endl;
#endif  // DEBUG_TRACES
    std::set<Row_index> effected_indices = three_clique_indices(indx);
    if (effected_indices.size() > 0) {
      for (auto idx = indx - 1; idx > 0; idx--) {
        Edge edge = std::get<0>(f_edge_vector_[idx]);
        Vertex_handle u = std::get<0>(edge);
        Vertex_handle v = std::get<1>(edge);
        // If idx is not critical so it should be processed, otherwise it stays in the graph
        if (!critical_edge_indicator_[idx]) {
          // If idx is affected
          if (effected_indices.find(idx) != effected_indices.end()) {
            if (!edge_is_dominated(edge)) {
#ifdef DEBUG_TRACES
              std::cout << "The curent index became critical " << idx  << std::endl;
#endif  // DEBUG_TRACES
              critical_edge_indicator_[idx] = true;
              filtered_edge_output({u, v}, filt);
              std::set<Row_index> inner_effected_indcs = three_clique_indices(idx);
              for (auto inr_idx = inner_effected_indcs.rbegin(); inr_idx != inner_effected_indcs.rend(); inr_idx++) {
                if (*inr_idx < idx) effected_indices.emplace(*inr_idx);
              }
#ifdef DEBUG_TRACES
              std::cout << "The following edge is critical with filt value: {" << u << "," << v << "}; "
                        << filt << std::endl;
#endif  // DEBUG_TRACES
            } else
              u_set_dominated_edges_.emplace(std::minmax(vertex_to_row_[u], vertex_to_row_[v]));
          } else
            // Idx is not affected hence dominated.
            u_set_dominated_edges_.emplace(std::minmax(vertex_to_row_[u], vertex_to_row_[v]));
        }
      }
    }
    u_set_dominated_edges_.clear();
  }

  // Returns list of non-zero columns of a particular indx.
  Row_indices_vector closed_neighbours_row_index(Row_index rw_u) const
  {
    Row_indices_vector non_zero_indices;
#ifdef DEBUG_TRACES
    std::cout << "The neighbours of the vertex: " << row_to_vertex_[rw_u] << " are. " << std::endl;
#endif  // DEBUG_TRACES
    // Iterate over the non-zero columns
    for (typename Sparse_row_matrix::InnerIterator it(sparse_row_adjacency_matrix_, rw_u); it; ++it) {
      Row_index rw_v = it.index();
      // If the vertex v is not dominated and the edge {u,v} is still in the matrix
      if (u_set_removed_edges_.find(std::minmax(rw_u, rw_v)) == u_set_removed_edges_.end() &&
          u_set_dominated_edges_.find(std::minmax(rw_u, rw_v)) == u_set_dominated_edges_.end()) {
        // inner index, here it is equal to it.columns()
        non_zero_indices.push_back(rw_v);
#ifdef DEBUG_TRACES
        std::cout << row_to_vertex_[rw_v] << ", " ;
#endif  // DEBUG_TRACES
      }
    }
#ifdef DEBUG_TRACES
    std::cout << std::endl;
#endif  // DEBUG_TRACES
    return non_zero_indices;
  }

  // Returns the list of closed neighbours of the edge :{u,v}.
  Row_indices_vector closed_common_neighbours_row_index(Row_index rw_u, Row_index rw_v) const
  {
    Row_indices_vector non_zero_indices_u = closed_neighbours_row_index(rw_u);
    Row_indices_vector non_zero_indices_v = closed_neighbours_row_index(rw_v);
    Row_indices_vector common;
    std::set_intersection(non_zero_indices_u.begin(), non_zero_indices_u.end(), non_zero_indices_v.begin(),
                          non_zero_indices_v.end(), std::inserter(common, common.begin()));

    return common;
  }

  // Insert a vertex in the data structure
  void insert_vertex(Vertex_handle vertex, Filtration_value filt_val) {
    auto result = vertex_to_row_.emplace(vertex, row_to_vertex_.size());
    // If it was not already inserted - Value won't be updated by emplace if it is already present
    if (result.second) {
      // Initializing the diagonal element of the adjency matrix corresponding to rw_b.
      sparse_row_adjacency_matrix_.insert(row_to_vertex_.size(), row_to_vertex_.size()) = filt_val;
      // Must be done after sparse_row_adjacency_matrix_ insertion
      row_to_vertex_.push_back(vertex);
    }
  }

  // Insert an edge in the data structure
  void insert_new_edge(Vertex_handle u, Vertex_handle v, Filtration_value filt_val)
  {
    // The edge must not be added before, it should be a new edge.
    insert_vertex(u, filt_val);
    if (u != v) {
      insert_vertex(v, filt_val);
#ifdef DEBUG_TRACES
      std::cout << "Insertion of the edge begins " << u <<", " << v << std::endl;
#endif  // DEBUG_TRACES

      auto rw_u = vertex_to_row_.find(u);
      auto rw_v = vertex_to_row_.find(v);
#ifdef DEBUG_TRACES
      std::cout << "Inserting the edge " << u <<", " << v << std::endl;
#endif  // DEBUG_TRACES
      sparse_row_adjacency_matrix_.insert(rw_u->second, rw_v->second) = filt_val;
      sparse_row_adjacency_matrix_.insert(rw_v->second, rw_u->second) = filt_val;
    }
#ifdef DEBUG_TRACES
    else {
     	std::cout << "Already a member simplex,  skipping..." << std::endl;
    }
#endif  // DEBUG_TRACES
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
  : f_edge_vector_(filtered_edge_range.begin(), filtered_edge_range.end()) {
    // To get the number of vertices
    std::unordered_set<Vertex_handle> vertices;
    for (Filtered_edge filtered_edge : filtered_edge_range) {
      Vertex_handle u;
      Vertex_handle v;
      std::tie(u,v) = std::get<0>(filtered_edge);
      vertices.emplace(u);
      vertices.emplace(v);
    }
    // Initializing sparse_row_adjacency_matrix_, This is a row-major sparse matrix.
    sparse_row_adjacency_matrix_ = Sparse_row_matrix(vertices.size(), vertices.size());
  }

  /** \brief Flag_complex_sparse_matrix constructor from a proximity graph, cf. `Gudhi::compute_proximity_graph`.
   *
   * @param[in] one_skeleton_graph The one skeleton graph. The graph must be in
   * `Flag_complex_sparse_matrix::Proximity_graph`.
   *
   * The constructor is computing and filling a vector of `Flag_complex_sparse_matrix::Filtered_edge`
   */
  Flag_complex_sparse_matrix(const Proximity_graph& one_skeleton_graph)
  : sparse_row_adjacency_matrix_(boost::num_vertices(one_skeleton_graph), boost::num_vertices(one_skeleton_graph)) {
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
    Row_index endIdx = 0;
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

    while (endIdx < f_edge_vector_.size()) {
      Filtered_edge fec = f_edge_vector_[endIdx];
      Edge edge = std::get<0>(fec);
      Vertex_handle u = std::get<0>(edge);
      Vertex_handle v = std::get<1>(edge);
      Filtration_value filt = std::get<1>(fec);

      // Inserts the edge in the sparse matrix to update the graph (G_i)
      insert_new_edge(u, v, filt);

      edge_to_index_map_.emplace(std::minmax(u, v), endIdx);
      critical_edge_indicator_.push_back(false);

      if (!edge_is_dominated(edge)) {
        critical_edge_indicator_[endIdx] = true;
        filtered_edge_output({u, v}, filt);
        if (endIdx > 1)
          set_edge_critical(endIdx, filt, filtered_edge_output);
      }
      endIdx++;
    }
  }

};

}  // namespace collapse

}  // namespace Gudhi

#endif  // FLAG_COMPLEX_SPARSE_MATRIX_H_
