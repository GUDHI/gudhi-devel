/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Siddharth Pritam, Marc Glisse
 *
 *    Copyright (C) 2020 Inria
 *
 *    Modification(s):
 *      - 2020/03 Vincent Rouvreau: integration to the gudhi library
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef FLAG_COMPLEX_EDGE_COLLAPSER_H_
#define FLAG_COMPLEX_EDGE_COLLAPSER_H_

#include <gudhi/Debug_utils.h>

#include <boost/functional/hash.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/container/flat_map.hpp>
#include <boost/container/flat_set.hpp>

#include <Eigen/Sparse>
#include <Eigen/src/Core/util/Macros.h>  // for EIGEN_VERSION_AT_LEAST

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
#include <type_traits>  // for std::decay

// Make compilation fail - required for external projects - https://github.com/GUDHI/gudhi-devel/issues/10
#if !EIGEN_VERSION_AT_LEAST(3,1,0)
# error Edge Collapse is only available for Eigen3 >= 3.1.0
#endif

namespace Gudhi {

namespace collapse {

/** \private
 * 
 * \brief Flag complex sparse matrix data structure.
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
   * @param[in] edges Range of Filtered edges range.There is no need the range to be sorted, as it will be performed in
   * `Flag_complex_edge_collapser::process_edges`.
   *
   * \tparam FilteredEdgeRange must be a range for which std::begin and std::end return iterators on a
   * `Flag_complex_edge_collapser::Filtered_edge`.
   */
  template<typename FilteredEdgeRange>
  Flag_complex_edge_collapser(const FilteredEdgeRange& edges)
  : f_edge_vector_(std::begin(edges), std::end(edges)) { }

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
      return (std::get<2>(edge_a) < std::get<2>(edge_b)); 
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

template<typename Vertex, typename Filtration_value>
struct Flag_complex_edge_collapser2 {
  using Filtered_edge = std::tuple<Vertex, Vertex, Filtration_value>;
  typedef std::pair<Vertex,Vertex> Edge;
  struct Cmpi { template<class T, class U> bool operator()(T const&a, U const&b)const{return b<a; } };
  typedef std::list<Edge> Edge_group; // vector? on erase elements...
  typedef std::map<Filtration_value, Edge_group, Cmpi> Events; // decreasing order
  Events events;
  typedef boost::container::flat_map<Vertex, Filtration_value> Ngb_list;
  //typedef std::unordered_map<Vertex, Ngb_list> Neighbors;
  typedef std::vector<Ngb_list> Neighbors; // we still need to find the edge in the group then...
  Neighbors neighbors;
  Vertex num_vertices;

  template<class FilteredEdgeRange>
  void read_edges(FilteredEdgeRange const&r){
    neighbors.resize(num_vertices);
    std::vector<typename Ngb_list::sequence_type> neighbors2(num_vertices);
    for(auto&&e : r){
      using std::get;
      Vertex u = get<0>(e);
      Vertex v = get<1>(e);
      Filtration_value f = get<2>(e);
      auto i = events.emplace_hint(events.end(), f, 0);
      neighbors2[u].emplace_back(v, f);
      neighbors2[v].emplace_back(u, f);
      i->second.push_back(std::minmax({u, v}));
    }
    for(int i=0;i<neighbors2.size();++i){
      neighbors[i].adopt_sequence(std::move(neighbors2[i])); // calls sort
    }
    neighbors2.clear();
  }

  void common_neighbors(boost::container::flat_set<Vertex>& e_ngb, std::vector<std::pair<Filtration_value, Vertex>>& e_ngb_later, Ngb_list const&nu, Ngb_list const&nv, Filtration_value f_event){
    auto ui = nu.begin();
    auto ue = nu.end();
    auto vi = nv.begin();
    auto ve = nv.end();
    assert(ui != ue && vi != ve);
    while(ui != ue && vi != ve){
      Vertex w = ui->first;
      if(w < vi->first) { ++ui; continue; }
      if(w > vi->first) { ++vi; continue; }
      Filtration_value f = std::max(ui->second, vi->second);
      if(f > f_event)
        e_ngb_later.emplace_back(f, w);
      else
        e_ngb.insert(e_ngb.end(), w);
      ++ui; ++vi;
    }
  }

  // Test if the neighborhood of e is included in the closed neighborhood of c
  template<class Ngb>
  bool is_dominated_by(Ngb const& e_ngb, Vertex c, Filtration_value f){
    Ngb_list const&nc = neighbors[c];
    // The best strategy probably depends on the distribution, how sparse / dense the adjacency matrix is, how (un)balanced the sizes of e_ngb and nc are.
    // Some efficient operations on sets work best with bitsets, although the need for a map complicates things.
#if 0
    // if few neighbors, use dichotomy? Seems slower.
    auto ci = nc.begin();
    auto ce = nc.end();
    for(auto v : e_ngb) {
      if(v==c)continue;
      ci = std::lower_bound(ci, ce, v, [](auto a, auto b){return a.first < b;});
      if(ci == nc.end() || ci->first != v || ci->second > f) return false;
    }
    return true;
#elif 0
    // I tried storing a copy of neighbors as a vector<absl::flat_hash_map> and using it for nc, but it was a bit slower here. It did help noticably with neighbors[dominator].find(w) in the main function though.
    for(auto v : e_ngb) {
      if(v==c)continue;
      auto it = nc.find(v);
      if(it == nc.end() || it->second > f) return false;
    }
    return true;
#else
    auto ci = nc.begin();
    auto ce = nc.end();
    auto eni = e_ngb.begin();
    auto ene = e_ngb.end();
    assert(eni != ene);
    assert(ci != ce);
    if(*eni == c && ++eni == ene) return true;
    for(;;){
      Vertex ve = *eni;
      Vertex vc = ci->first;
      while(ve > vc) {
        // try a gallop strategy (exponential search)? Seems slower
        if(++ci == ce) return false;
        vc = ci->first;
      }
      if(ve < vc) return false;
      // ve == vc
      if(ci->second > f) return false;
      if(++eni == ene)return true;
      // Should we store a closed neighborhood of c (including c) so we can avoid testing for c at each iteration?
      if(*eni == c && ++eni == ene)return true;
      if(++ci == ce) return false;
    }
#endif
  }

  // delay = [](double d){return d*1.05;}
  template<class FilteredEdgeRange, class Delay>
  void process_edges(FilteredEdgeRange const& edges, Delay&& delay) {
    {
      Vertex maxi = 0, maxj = 0;
      for(auto& fe : edges) {
        Vertex i = std::get<0>(fe);
        Vertex j = std::get<1>(fe);
        if (i > maxi) maxi = i;
        if (j > maxj) maxj = j;
      }
      num_vertices = std::max(maxi, maxj) + 1;
    }

    read_edges(edges);

    typename Events::iterator evi;
    boost::container::flat_set<Vertex> e_ngb;
    e_ngb.reserve(num_vertices);
    //std::multimap<Filtration_value, Vertex> e_ngb_later;
    std::vector<std::pair<Filtration_value, Vertex>> e_ngb_later;
    // do not use reverse_iterator, it behaves badly with erase.
    for(auto next_evi = events.begin(); (evi = next_evi) != events.end(); ){
      ++next_evi; // do it now, in case we remove evi
      Edge_group& eg = evi->second;
      //Filtration_value f_event = evi->first;
      auto next_ei = eg.begin();
      auto ei = next_ei;
      for(; (ei = next_ei) != eg.end();){
        next_ei = std::next(ei); // do it now, in case we remove ei
        Vertex u = ei->first;
        Vertex v = ei->second;
        auto time = evi;
        Filtration_value a_bit_later = delay(evi->first);
        if(a_bit_later != evi->first) {
#if 0
        auto next_time = time;
        while (next_time != events.begin() && (--next_time)->first <= a_bit_later)
          time = next_time;
#else
        // too bad there isn't a map::lower_bound_after
        time = events.lower_bound(a_bit_later);
#endif
        }
        auto start_time = time;
        e_ngb.clear();
        e_ngb_later.clear();
        common_neighbors(e_ngb, e_ngb_later, neighbors[u], neighbors[v], time->first);
        auto cmp1=[](auto const&a, auto const&b){return a.first > b.first;};
        auto e_ngb_later_begin=e_ngb_later.begin();
        auto e_ngb_later_end=e_ngb_later.end();
        bool heapified = false;

        bool dead = false;
        while(true) {
          Vertex dominator = -1;
          // special case for size 1
          // if(e_ngb.size()==1){dominator=*e_ngb.begin();}else
          // TODO: try testing dominators in order of filtration value
          for(auto c : e_ngb){
            if(is_dominated_by(e_ngb, c, time->first)){
              dominator = c;
              break;
            }
          }
          if(dominator==-1) break;
          bool still_dominated;
          do {
            if(e_ngb_later_begin == e_ngb_later_end) {
              dead = true; goto end_move;
            }
            if(!heapified) {
              // Eagerly sorting can be slow
              std::make_heap(e_ngb_later_begin, e_ngb_later_end, cmp1);
              heapified=true;
            }
            // go back to storing an iterator in e_ngb_later?
            time = events.lower_bound(e_ngb_later_begin->first);
            still_dominated = true;
            while (e_ngb_later_begin != e_ngb_later_end && e_ngb_later_begin->first <= time->first) {
              Vertex w = e_ngb_later_begin->second;
              auto wit = neighbors[dominator].find(w);
              if (wit == neighbors[dominator].end() || wit->second > e_ngb_later_begin->first)
                still_dominated = false;
              e_ngb.insert(w);
              std::pop_heap(e_ngb_later_begin, e_ngb_later_end--, cmp1);
            }
          } while (still_dominated); // this doesn't seem to help that much...
        }
end_move:
        if(dead) {
          neighbors[u].erase(v);
          neighbors[v].erase(u);
          eg.erase(ei);
        } else if(start_time != time){
          time->second.push_back(*ei);
          neighbors[u][v]=time->first;
          neighbors[v][u]=time->first;
          eg.erase(ei);
        }
      }
      // check if event is dead (all edges have moved)
      if (evi->second.empty()) {
        events.erase(evi);
      }
    }
  }

  std::vector<Filtered_edge> output() {
    std::vector<std::tuple<Vertex, Vertex, Filtration_value>> r;
    for(auto& ev : events)
      for(auto& e : ev.second)
        r.emplace_back(e.first, e.second, ev.first);
    return r;
  }

};

/** \brief Implicitly constructs a flag complex from edges as an input, collapses edges while preserving the persistent
 * homology and returns the remaining edges as a range.
 *
 * \param[in] edges Range of Filtered edges.There is no need the range to be sorted, as it will be performed.
 *
 * \tparam FilteredEdgeRange furnishes `std::begin` and `std::end` methods and returns an iterator on a
 * FilteredEdge of type `std::tuple<Vertex_handle, Vertex_handle, Filtration_value>` where `Vertex_handle` is the type
 * of a vertex index and `Filtration_value` is the type of an edge filtration value.
 *
 * \return Remaining edges after collapse as a range of
 * `std::tuple<Vertex_handle, Vertex_handle, Filtration_value>`.
 * 
 * \ingroup edge_collapse
 * 
 */
template<class FilteredEdgeRange> auto flag_complex_collapse_edges(const FilteredEdgeRange& edges) {
  auto first_edge_itr = std::begin(edges);
  using Vertex_handle = std::decay_t<decltype(std::get<0>(*first_edge_itr))>;
  using Filtration_value = std::decay_t<decltype(std::get<2>(*first_edge_itr))>;
  using Edge_collapser = Flag_complex_edge_collapser<Vertex_handle, Filtration_value>;
  std::vector<typename Edge_collapser::Filtered_edge> remaining_edges;
  if (first_edge_itr != std::end(edges)) {
    Edge_collapser edge_collapser(edges);
    edge_collapser.process_edges(
      [&remaining_edges](Vertex_handle u, Vertex_handle v, Filtration_value filtration) {
          // insert the edge
          remaining_edges.emplace_back(u, v, filtration);
        });
  }
  return remaining_edges;
}

// Would it help to label the points according to some spatial sorting?
template<class FilteredEdgeRange, class Delay> auto flag_complex_collapse_edges2(FilteredEdgeRange&& edges, Delay&&delay) {
  auto first_edge_itr = std::begin(edges);
  using Vertex = std::decay_t<decltype(std::get<0>(*first_edge_itr))>;
  using Filtration_value = std::decay_t<decltype(std::get<2>(*first_edge_itr))>;
  using Edge_collapser = Flag_complex_edge_collapser2<Vertex, Filtration_value>;
  if (first_edge_itr != std::end(edges)) {
    std::vector<typename Edge_collapser::Filtered_edge> edges2;
    if(std::is_same<FilteredEdgeRange, std::vector<typename Edge_collapser::Filtered_edge>>::value)
      edges2 = std::move(edges);
    else
      edges2.insert(edges2.end(), edges.begin(), edges.end());
    // The sorting is not necessary, but inserting in the map is faster if done in order
    std::sort(edges2.begin(), edges2.end(), [](auto const&a, auto const&b){return std::get<2>(a)>std::get<2>(b);});
    Edge_collapser edge_collapser;
    edge_collapser.process_edges(edges2, std::forward<Delay>(delay));
    return edge_collapser.output();
  }
  return std::vector<typename Edge_collapser::Filtered_edge>();
}
template<class FilteredEdgeRange> auto flag_complex_collapse_edges2(const FilteredEdgeRange& edges) {
  return flag_complex_collapse_edges2(edges, [](auto const&d){return d;});
}

}  // namespace collapse

}  // namespace Gudhi

#endif  // FLAG_COMPLEX_EDGE_COLLAPSER_H_
