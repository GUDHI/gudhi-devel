/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Siddharth Pritam
 *
 *    Copyright (C) 2018 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef FLAG_COMPLEX_SPARSE_MATRIX_H_
#define FLAG_COMPLEX_SPARSE_MATRIX_H_

#include <gudhi/Rips_edge_list.h>
#include <gudhi/graph_simplicial_complex.h>

#include <boost/functional/hash.hpp>

#include <Eigen/Sparse>

#ifdef GUDHI_USE_TBB
#include <tbb/parallel_sort.h>
#endif

#include <iostream>
#include <utility>
#include <vector>
#include <queue>
#include <unordered_map>
#include <tuple>
#include <list>
#include <algorithm>
#include <chrono>

#include <ctime>
#include <fstream>

namespace Gudhi {

namespace collapse {


typedef std::size_t Vertex;
using Edge = std::pair<Vertex, Vertex>;  // This is an ordered pair, An edge is stored with convention of the first
                                         // element being the smaller i.e {2,3} not {3,2}. However this is at the level
                                         // of row indices on actual vertex lables
using Filtered_edge = std::pair<Edge, double>;

using Map_vertex_to_index = std::unordered_map<Vertex, std::size_t>;

using Sparse_row_matrix = Eigen::SparseMatrix<double, Eigen::RowMajor>;

using doubleVector = std::vector<double>;

using Filtered_sorted_edge_list = std::vector<std::tuple<double, Vertex, Vertex>>;

//!  Class SparseMsMatrix
/*!
  The class for storing the Vertices v/s MaxSimplices Sparse Matrix and performing collapses operations using the N^2()
  Algorithm.
*/
class Flag_complex_sparse_matrix {
 private:
  std::unordered_map<int, Vertex> row_to_vertex_;

  // Vertices stored as an unordered_set
  std::unordered_set<Vertex> vertices_;

  // Unordered set of  removed edges. (to enforce removal from the matrix)
  std::unordered_set<Edge, boost::hash<Edge>> u_set_removed_redges_;

  // Unordered set of  dominated edges. (to inforce removal from the matrix)
  std::unordered_set<Edge, boost::hash<Edge>> u_set_dominated_redges_;

  // Map from egde to its index
  std::unordered_map<Edge, std::size_t, boost::hash<Edge>> edge_to_index_map_;
  // Boolean vector to indicate if the index is critical or not.
  std::vector<bool> critical_edge_indicator_;  // critical indicator

  // Boolean vector to indicate if the index is critical or not.
  std::vector<bool> dominated_edge_indicator_;  // domination indicator

  //! Stores the Map between vertices_<B>row_to_vertex_  and row indices <B>row_to_vertex_ -> row-index</B>.
  /*!
    So, if the original simplex tree had vertices_ 0,1,4,5 <br>
    <B>row_to_vertex_</B> would store : <br>
    \verbatim
    Values =  | 0 | 1 | 4 | 5 |
    Indices =   0   1   2   3
    \endverbatim
    And <B>vertex_to_row_</B> would be a map like the following : <br>
    \verbatim
    0 -> 0
    1 -> 1
    4 -> 2
    5 -> 3
    \endverbatim
  */
  std::unordered_map<Vertex, std::size_t> vertex_to_row_;

  //! Stores the Sparse matrix of double values representing the Original Simplicial Complex.
  /*!
    \code
    Sparse_row_matrix   = Eigen::SparseMatrix<double, Eigen::RowMajor> ;
    \endcode
    ;
      */

  Sparse_row_matrix sparse_row_adjacency_matrix_;  // This is row-major version of the same sparse-matrix, to facilitate easy access
                                       // to elements when traversing the matrix row-wise.

  //! Stores <I>true</I> for dominated rows and  <I>false</I> for undominated rows.
  /*!
    Initialised to a vector of length equal to the value of the variable <B>rows</B> with all <I>false</I> values.
    Subsequent removal of dominated vertices_ is reflected by concerned entries changing to <I>true</I> in this vector.
  */
  std::vector<bool> domination_indicator_;  //(domination indicator)

  // Vector of filtered edges, for edge-collapse, the indices of the edges are the row-indices.
  std::vector<Filtered_edge> f_edge_vector_;

  // Stores the indices from the sorted filtered edge vector.
  // std::set<std::size_t> recurCriticalCoreIndcs;

  //! Stores the number of vertices_ in the original Simplicial Complex.
  /*!
    This stores the count of vertices_ (which is also the number of rows in the Matrix).
  */
  std::size_t rows;

  // Edge e is the actual edge (u,v). Not the row ids in the matrixs
  bool check_edge_domination(Edge e)
  {
    auto u = std::get<0>(e);
    auto v = std::get<1>(e);

    auto rw_u = vertex_to_row_[u];
    auto rw_v = vertex_to_row_[v];
    auto rw_e = std::make_pair(rw_u, rw_v);
#ifdef DEBUG_TRACES
    std::cout << "The edge {" << u << ", " << v <<  "} is going for domination check." << std::endl;
#endif  // DEBUG_TRACES
    auto common_neighbours = closed_common_neighbours_row_index(rw_e);
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
          if (rw_c != rw_u and rw_c != rw_v) {
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

  // The edge should be sorted by the indices and indices are original
  bool check_domination_indicator(Edge e)
  {
    return dominated_edge_indicator_[edge_to_index_map_[e]];
  }

  std::set<std::size_t> three_clique_indices(std::size_t crit) {
    std::set<std::size_t> edge_indices;

    Edge e = std::get<0>(f_edge_vector_[crit]);
    Vertex u = std::get<0>(e);
    Vertex v = std::get<1>(e);

#ifdef DEBUG_TRACES
    std::cout << "The  current critical edge to re-check criticality with filt value is : f {" << u << "," << v
              << "} = " << std::get<1>(f_edge_vector_[crit]) << std::endl;
#endif  // DEBUG_TRACES
    auto rw_u = vertex_to_row_[u];
    auto rw_v = vertex_to_row_[v];
    auto rw_critical_edge = std::make_pair(rw_u, rw_v);

    doubleVector common_neighbours = closed_common_neighbours_row_index(rw_critical_edge);

    if (common_neighbours.size() > 2) {
      for (auto rw_c : common_neighbours) {
        if (rw_c != rw_u and rw_c != rw_v) {
          auto e_with_new_nbhr_v = std::minmax(u, row_to_vertex_[rw_c]);
          auto e_with_new_nbhr_u = std::minmax(v, row_to_vertex_[rw_c]);
          edge_indices.emplace(edge_to_index_map_[e_with_new_nbhr_v]);
          edge_indices.emplace(edge_to_index_map_[e_with_new_nbhr_u]);
        }
      }
    }
    return edge_indices;
  }

  template<typename FilteredEdgeInsertion>
  void set_edge_critical(std::size_t indx, double filt, FilteredEdgeInsertion filtered_edge_insert) {
#ifdef DEBUG_TRACES
    std::cout << "The curent index  with filtration value " << indx << ", " << filt << " is primary critical" <<
    std::endl;
#endif  // DEBUG_TRACES
    std::set<std::size_t> effectedIndcs = three_clique_indices(indx);
    if (effectedIndcs.size() > 0) {
      for (auto idx = indx - 1; idx > 0; idx--) {
        Edge e = std::get<0>(f_edge_vector_[idx]);
        Vertex u = std::get<0>(e);
        Vertex v = std::get<1>(e);
        // If idx is not critical so it should be proceses, otherwise it stays in the graph // prev
        // code : recurCriticalCoreIndcs.find(idx) == recurCriticalCoreIndcs.end()
        if (not critical_edge_indicator_[idx]) {
          // If idx is affected
          if (effectedIndcs.find(idx) != effectedIndcs.end()) {
            if (not check_edge_domination(e)) {
#ifdef DEBUG_TRACES
              std::cout << "The curent index became critical " << idx  << std::endl;
#endif  // DEBUG_TRACES
              critical_edge_indicator_[idx] = true;
              filtered_edge_insert({u, v}, filt);
              std::set<std::size_t> inner_effected_indcs = three_clique_indices(idx);
              for (auto inr_idx = inner_effected_indcs.rbegin(); inr_idx != inner_effected_indcs.rend(); inr_idx++) {
                if (*inr_idx < idx) effectedIndcs.emplace(*inr_idx);
              }
              inner_effected_indcs.clear();
#ifdef DEBUG_TRACES
              std::cout << "The following edge is critical with filt value: {" << std::get<0>(e) << "," <<
              std::get<1>(e) << "}; " << filt << std::endl;
#endif  // DEBUG_TRACES
            } else
              u_set_dominated_redges_.emplace(std::minmax(vertex_to_row_[u], vertex_to_row_[v]));
          } else
            // Idx is not affected hence dominated.
            u_set_dominated_redges_.emplace(std::minmax(vertex_to_row_[u], vertex_to_row_[v]));
        }
      }
    }
    effectedIndcs.clear();
    u_set_dominated_redges_.clear();
  }

  // Returns list of non-zero columns of the particular indx.
  doubleVector closed_neighbours_row_index(double indx)
  {
    doubleVector non_zero_indices;
    Vertex u = indx;
    Vertex v;
#ifdef DEBUG_TRACES
    std::cout << "The neighbours of the vertex: " << row_to_vertex_[u] << " are. " << std::endl;
#endif  // DEBUG_TRACES
    if (not domination_indicator_[indx]) {
      // Iterate over the non-zero columns
      for (Sparse_row_matrix::InnerIterator it(sparse_row_adjacency_matrix_, indx); it; ++it) {
        v = it.index();
        // If the vertex v is not dominated and the edge {u,v} is still in the matrix
        if (not domination_indicator_[v] and u_set_removed_redges_.find(std::minmax(u, v)) == u_set_removed_redges_.end() and
            u_set_dominated_redges_.find(std::minmax(u, v)) == u_set_dominated_redges_.end()) {
          // inner index, here it is equal to it.columns()
          non_zero_indices.push_back(it.index());
#ifdef DEBUG_TRACES
          std::cout << row_to_vertex_[it.index()] << ", " ;
#endif  // DEBUG_TRACES
        }
      }
#ifdef DEBUG_TRACES
      std::cout << std::endl;
#endif  // DEBUG_TRACES
    }
    return non_zero_indices;
  }

  doubleVector closed_common_neighbours_row_index(Edge e)  // Returns the list of closed neighbours of the edge :{u,v}.
  {
    doubleVector common;
    doubleVector non_zero_indices_u;
    doubleVector non_zero_indices_v;
    double u = std::get<0>(e);
    double v = std::get<1>(e);

    non_zero_indices_u = closed_neighbours_row_index(u);
    non_zero_indices_v = closed_neighbours_row_index(v);
    std::set_intersection(non_zero_indices_u.begin(), non_zero_indices_u.end(), non_zero_indices_v.begin(),
                          non_zero_indices_v.end(), std::inserter(common, common.begin()));

    return common;
  }
  
  void insert_vertex(const Vertex& vertex, double filt_val) {
    auto rw = vertex_to_row_.find(vertex);
    if (rw == vertex_to_row_.end()) {
      // Initializing the diagonal element of the adjency matrix corresponding to rw_b.
      sparse_row_adjacency_matrix_.insert(rows, rows) = filt_val;
      domination_indicator_.push_back(false);
      vertex_to_row_.insert(std::make_pair(vertex, rows));
      row_to_vertex_.insert(std::make_pair(rows, vertex));
      rows++;
    }
  }

  void insert_new_edges(const Vertex& u, const Vertex& v, double filt_val)
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
  //! Main Constructor
  /*!
    Argument is an instance of Filtered_sorted_edge_list. <br>
    This is THE function that initialises all data members to appropriate values. <br>
    <B>row_to_vertex_</B>, <B>vertex_to_row_</B>, <B>rows</B>, <B>cols</B>, <B>sparse_row_adjacency_matrix_</B> are initialised here.
    <B>domination_indicator_</B> are initialised by init() function which is
    called at the begining of this. <br>
  */
  Flag_complex_sparse_matrix(const Filtered_sorted_edge_list& edge_t)
  : rows(0) {
    for (size_t bgn_idx = 0; bgn_idx < edge_t.size(); bgn_idx++) {
      Vertex u = std::get<1>(edge_t[bgn_idx]);
      Vertex v = std::get<2>(edge_t[bgn_idx]);
      f_edge_vector_.push_back({{u, v}, std::get<0>(edge_t[bgn_idx])});
      vertices_.emplace(u);
      vertices_.emplace(v);
    }
  }

  template<typename OneSkeletonGraph>
  Flag_complex_sparse_matrix(const OneSkeletonGraph& one_skeleton_graph)
  : rows(0) {
    // Insert all vertices_
    for (auto v_it = boost::vertices_(one_skeleton_graph); v_it.first != v_it.second; ++v_it.first) {
      vertices_.emplace(*(v_it.first));
    }
    // Insert all edges
    for (auto edge_it = boost::edges(one_skeleton_graph);
         edge_it.first != edge_it.second; ++edge_it.first) {
      auto edge = *(edge_it.first);
      Vertex u = source(edge, one_skeleton_graph);
      Vertex v = target(edge, one_skeleton_graph);
      f_edge_vector_.push_back({{u, v}, boost::get(Gudhi::edge_filtration_t(), one_skeleton_graph, edge)});
    }
    // Sort edges
    auto sort_by_filtration = [](const Filtered_edge& edge_a, const Filtered_edge& edge_b) -> bool
    {
      return (get<1>(edge_a) < get<1>(edge_b)); 
    };

#ifdef GUDHI_USE_TBB
    tbb::parallel_sort(f_edge_vector_.begin(), f_edge_vector_.end(), sort_by_filtration);
#else
    std::stable_sort(f_edge_vector_.begin(), f_edge_vector_.end(), sort_by_filtration);
#endif
  }

  // Performs edge collapse in a decreasing sequence of the filtration value.
  template<typename FilteredEdgeInsertion>
  void filtered_edge_collapse(FilteredEdgeInsertion filtered_edge_insert) {
    std::size_t endIdx = 0;

    u_set_removed_redges_.clear();
    u_set_dominated_redges_.clear();
    critical_edge_indicator_.clear();

    // Initializing sparse_row_adjacency_matrix_, This is a row-major sparse matrix.
    sparse_row_adjacency_matrix_ = Sparse_row_matrix(vertices_.size(), vertices_.size());

    while (endIdx < f_edge_vector_.size()) {
      Filtered_edge fec = f_edge_vector_[endIdx];
      Edge e = std::get<0>(fec);
      Vertex u = std::get<0>(e);
      Vertex v = std::get<1>(e);
      double filt = std::get<1>(fec);

      // Inserts the edge in the sparse matrix to update the graph (G_i)
      insert_new_edges(u, v, filt);

      edge_to_index_map_.emplace(std::minmax(u, v), endIdx);
      critical_edge_indicator_.push_back(false);
      dominated_edge_indicator_.push_back(false);

      if (not check_edge_domination(e)) {
        critical_edge_indicator_[endIdx] = true;
        dominated_edge_indicator_[endIdx] = false;
        filtered_edge_insert({u, v}, filt);
        if (endIdx > 1)
          set_edge_critical(endIdx, filt, filtered_edge_insert);
      } else
        dominated_edge_indicator_[endIdx] = true;
      endIdx++;
    }
  }

  std::size_t num_vertices() const { return vertices_.size(); }

};

}  // namespace collapse

}  // namespace Gudhi

#endif  // FLAG_COMPLEX_SPARSE_MATRIX_H_
