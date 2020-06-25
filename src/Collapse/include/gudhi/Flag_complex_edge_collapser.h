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

#ifndef FLAG_COMPLEX_EDGE_COLLAPSER_H_
#define FLAG_COMPLEX_EDGE_COLLAPSER_H_

#include <gudhi/graph_simplicial_complex.h>
#include <gudhi/Debug_utils.h>

#include <boost/functional/hash.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/iterator/iterator_facade.hpp>

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
 * \class Flag_complex_edge_collapser
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
class Flag_complex_edge_collapser {
 public:
  /** \brief Re-define Vertex as Vertex_handle type to ease the interface with `Gudhi::Proximity_graph`. */
  using Vertex_handle = Vertex;
  /** \brief Re-define Filtration as Filtration_value type to ease the interface with `Gudhi::Proximity_graph`. */
  using Filtration_value = Filtration;

 private:
  // internal numbering of vertices and edges
  using IVertex = std::size_t;
  using Edge_index = std::size_t;
  using IEdge = std::pair<IVertex, IVertex>;

  // The sparse matrix data type
  // (Eigen::SparseMatrix<Edge_index, Eigen::RowMajor> has slow insertions)
  using Sparse_vector = Eigen::SparseVector<Edge_index>;
  using Sparse_row_matrix = std::vector<Sparse_vector>;

  // Range of neighbors of a vertex
  template<bool closed>
  struct Neighbours {
    class iterator : public boost::iterator_facade<iterator,
      IVertex, /* value_type */
      std::input_iterator_tag, // or boost::single_pass_traversal_tag
      IVertex /* reference */ >
    {
      public:
        iterator():ptr(nullptr){}
        iterator(Neighbours const*p):ptr(p){find_valid();}
      private:
        friend class boost::iterator_core_access;
        Neighbours const*ptr;
        void increment(){
          ++ptr->it;
          find_valid();
        }
        void find_valid(){
          auto& it = ptr->it;
          do {
            if(!it) { ptr=nullptr; break; }
            if(IVertex(it.index()) == ptr->u) {
              if(closed) break;
              else continue;
            }
            Edge_index e = it.value();
            if(e <= ptr->ec->current_backward || ptr->ec->critical_edge_indicator_[e]) break;
          } while(++it, true);
        }
        bool equal(iterator const& other) const { return ptr == other.ptr; }
        IVertex dereference() const { return ptr->it.index(); }
    };
    typedef iterator const_iterator;
    mutable typename Sparse_vector::InnerIterator it;
    Flag_complex_edge_collapser const*ec;
    IVertex u;
    iterator begin() const { return this; }
    iterator end() const { return {}; }
    explicit Neighbours(Flag_complex_edge_collapser const*p,IVertex u):it(p->sparse_row_adjacency_matrix_[u]),ec(p),u(u){}
  };

  // A range of row indices
  using IVertex_vector = std::vector<IVertex>;

 public:
  /** \brief Filtered_edge is a type to store an edge with its filtration value. */
  using Filtered_edge = std::tuple<Vertex_handle, Vertex_handle, Filtration_value>;

 private:
  // Map from row index to its vertex handle
  std::vector<Vertex_handle> row_to_vertex_;

  // Index of the current edge in the backwards walk. Edges <= current_backward are part of the temporary graph,
  // while edges > current_backward are removed unless critical_edge_indicator_.
  Edge_index current_backward = -1;

  // Map from IEdge to its index
  std::unordered_map<IEdge, Edge_index, boost::hash<IEdge>> iedge_to_index_map_;

  // Boolean vector to indicate if the edge is critical.
  std::vector<bool> critical_edge_indicator_;

  // Map from vertex handle to its row index
  std::unordered_map<Vertex_handle, IVertex> vertex_to_row_;

  // Stores the Sparse matrix of Filtration values representing the original graph.
  // The matrix rows and columns are indexed by IVertex.
  Sparse_row_matrix sparse_row_adjacency_matrix_;

  // The input, a vector of filtered edges.
  std::vector<Filtered_edge> f_edge_vector_;

  // Edge is the actual edge (u,v), with Vertex_handle u and v, not IVertex.
  bool edge_is_dominated(Vertex_handle u, Vertex_handle v) const
  {
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
        auto neighbours_c = neighbours_row_index<true>(rw_c);
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

    Vertex_handle u = std::get<0>(f_edge_vector_[crit]);
    Vertex_handle v = std::get<1>(f_edge_vector_[crit]);

#ifdef DEBUG_TRACES
    std::cout << "The  current critical edge to re-check criticality with filt value is : f {" << u << "," << v
              << "} = " << std::get<2>(f_edge_vector_[crit]) << std::endl;
#endif  // DEBUG_TRACES
    auto rw_u = vertex_to_row_[u];
    auto rw_v = vertex_to_row_[v];

    IVertex_vector common_neighbours = open_common_neighbours_row_index(rw_u, rw_v);

    for (auto rw_c : common_neighbours) {
      IEdge e_with_new_nbhr_v = std::minmax(rw_u, rw_c);
      IEdge e_with_new_nbhr_u = std::minmax(rw_v, rw_c);
      edge_indices.emplace(iedge_to_index_map_[e_with_new_nbhr_v]);
      edge_indices.emplace(iedge_to_index_map_[e_with_new_nbhr_u]);
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
      Vertex_handle u = std::get<0>(f_edge_vector_[current_backward]);
      Vertex_handle v = std::get<1>(f_edge_vector_[current_backward]);
      // If current_backward is not critical so it should be processed, otherwise it stays in the graph
      if (!critical_edge_indicator_[current_backward]) {
        if (!edge_is_dominated(u, v)) {
#ifdef DEBUG_TRACES
          std::cout << "The curent index became critical " << current_backward  << std::endl;
#endif  // DEBUG_TRACES
          critical_edge_indicator_[current_backward] = true;
          filtered_edge_output(u, v, filt);
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
  template<bool closed>
  auto neighbours_row_index(IVertex rw_u) const
  {
    return Neighbours<closed>(this, rw_u);
  }

  // Returns the list of open neighbours of the edge :{u,v}.
  IVertex_vector open_common_neighbours_row_index(IVertex rw_u, IVertex rw_v) const
  {
    auto non_zero_indices_u = neighbours_row_index<false>(rw_u);
    auto non_zero_indices_v = neighbours_row_index<false>(rw_v);
    IVertex_vector common;
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
  IEdge insert_new_edge(Vertex_handle u, Vertex_handle v, Edge_index idx)
  {
    GUDHI_CHECK((u != v), std::invalid_argument("Flag_complex_edge_collapser::insert_new_edge with u == v"));
    // The edge must not be added before, it should be a new edge.
    IVertex rw_u = insert_vertex(u);
    IVertex rw_v = insert_vertex(v);
#ifdef DEBUG_TRACES
    std::cout << "Inserting the edge " << u <<", " << v << std::endl;
#endif  // DEBUG_TRACES
    sparse_row_adjacency_matrix_[rw_u].insert(rw_v) = idx;
    sparse_row_adjacency_matrix_[rw_v].insert(rw_u) = idx;
    return std::minmax(rw_u, rw_v);
  }

 public:
  /** \brief Flag_complex_edge_collapser constructor from a range of filtered edges.
   *
   * @param[in] begin Iterator on the first element of a filtered edges range, aka. `std::begin`. Filtered edges must
   * be in `Flag_complex_edge_collapser::Filtered_edge`.
   *
   * @param[in] end Iterator on the final element of a filtered edges range, aka. `std::end`. Filtered edges must be
   * in `Flag_complex_edge_collapser::Filtered_edge`.
   *
   * There is no need the range to be sorted, as it will be performed in
   * `Flag_complex_edge_collapser::process_edges`.
   */
  template<typename Filtered_edge_iterator>
  Flag_complex_edge_collapser(Filtered_edge_iterator begin, Filtered_edge_iterator end)
  : f_edge_vector_(begin, end) { }

  /** \brief Inserts all edges given by a OneSkeletonGraph into a vector of
   * `Flag_complex_edge_collapser::Filtered_edge`.
   * OneSkeletonGraph must be a model of
   * <a href="http://www.boost.org/doc/libs/1_65_1/libs/graph/doc/EdgeListGraph.html">boost::EdgeListGraph</a>
   * and <a href="http://www.boost.org/doc/libs/1_65_1/libs/graph/doc/PropertyGraph.html">boost::PropertyGraph</a>.
   *
   * The edge filtration value is accessible through the property tag
   * edge_filtration_t.
   *
   * boost::graph_traits<OneSkeletonGraph>::vertex_descriptor
   *                                    must be Vertex_handle.
   * boost::graph_traits<OneSkeletonGraph>::directed_category
   *                                    can be directed_tag (the fastest, the least RAM use), undirected_tag or even
   *                                    bidirected_tag.
   *
   * If an edge appears with multiplicity, the function will arbitrarily pick one representative to read the filtration
   * value.
   * 
   * `Gudhi::Proximity_graph<Flag_complex_edge_collapser>` is a good candidate for OneSkeletonGraph.
   */
  template<class OneSkeletonGraph>
  Flag_complex_edge_collapser(const OneSkeletonGraph& one_skeleton_graph) {
    // Insert all edges
    for (auto edge_it = edges(one_skeleton_graph);
         edge_it.first != edge_it.second; ++edge_it.first) {
      auto edge = *(edge_it.first);
      Vertex_handle u = source(edge, one_skeleton_graph);
      Vertex_handle v = target(edge, one_skeleton_graph);
      f_edge_vector_.emplace_back(u, v, get(Gudhi::edge_filtration_t(), one_skeleton_graph, edge));
    }
  }

  /** \brief Performs edge collapse in a increasing sequence of the filtration value.
   *
   * \tparam filtered_edge_output is a functor that is called on the output edges, in non-decreasing order of
   * filtration, as filtered_edge_output(u, v, f) where u and v are Vertex_handle representing the extremities of the
   * edge, and f is its new Filtration_value.
   */
  template<typename FilteredEdgeOutput>
  void process_edges(FilteredEdgeOutput filtered_edge_output) {
    // Sort edges
    auto sort_by_filtration = [](const Filtered_edge& edge_a, const Filtered_edge& edge_b) -> bool
    {
      return (get<2>(edge_a) < get<2>(edge_b)); 
    };

#ifdef GUDHI_USE_TBB
    tbb::parallel_sort(f_edge_vector_.begin(), f_edge_vector_.end(), sort_by_filtration);
#else
    std::sort(f_edge_vector_.begin(), f_edge_vector_.end(), sort_by_filtration);
#endif

    for (Edge_index endIdx = 0; endIdx < f_edge_vector_.size(); endIdx++) {
      Filtered_edge fec = f_edge_vector_[endIdx];
      Vertex_handle u = std::get<0>(fec);
      Vertex_handle v = std::get<1>(fec);
      Filtration_value filt = std::get<2>(fec);

      // Inserts the edge in the sparse matrix to update the graph (G_i)
      IEdge ie = insert_new_edge(u, v, endIdx);

      iedge_to_index_map_.emplace(ie, endIdx);
      critical_edge_indicator_.push_back(false);

      if (!edge_is_dominated(u, v)) {
        critical_edge_indicator_[endIdx] = true;
        filtered_edge_output(u, v, filt);
        if (endIdx > 1)
          set_edge_critical(endIdx, filt, filtered_edge_output);
      }
    }
  }

};

}  // namespace collapse

}  // namespace Gudhi

#endif  // FLAG_COMPLEX_EDGE_COLLAPSER_H_
