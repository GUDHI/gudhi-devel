/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Siddharth Pritam, Marc Glisse
 *
 *    Copyright (C) 2020 Inria
 *
 *    Modification(s):
 *      - 2020/03 Vincent Rouvreau: integration to the gudhi library
 *      - 2021 Marc Glisse: complete rewrite
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef FLAG_COMPLEX_EDGE_COLLAPSER_H_
#define FLAG_COMPLEX_EDGE_COLLAPSER_H_

#include <gudhi/Debug_utils.h>

#include <boost/container/flat_map.hpp>
#include <boost/container/flat_set.hpp>

#ifdef GUDHI_USE_TBB
#include <tbb/parallel_sort.h>
#endif

#include <utility>
#include <vector>
#include <tuple>
#include <algorithm>
#include <limits>

namespace Gudhi {

namespace collapse {

/** \private
 *
 * \brief Flag complex sparse matrix data structure.
 *
 * \tparam Vertex type must be an integer type.
 * \tparam Filtration type for the value of the filtration function.
 */
template<typename Vertex, typename Filtration_value>
struct Flag_complex_edge_collapser {
  using Filtered_edge = std::tuple<Vertex, Vertex, Filtration_value>;
  typedef std::pair<Vertex,Vertex> Edge;
  struct Cmpi { template<class T, class U> bool operator()(T const&a, U const&b)const{return b<a; } };
  typedef boost::container::flat_map<Vertex, Filtration_value> Ngb_list;
  typedef std::vector<Ngb_list> Neighbors;
  Neighbors neighbors; // closed neighborhood
  std::size_t num_vertices;
  std::vector<std::tuple<Vertex, Vertex, Filtration_value>> res;

#ifdef GUDHI_COLLAPSE_USE_DENSE_ARRAY
  // Minimal matrix interface
  // Using this matrix generally helps performance, but the memory use may be excessive for a very sparse graph
  // (and in extreme cases the constant initialization of the matrix may start to dominate the running time).
  // Are there cases where the matrix is too big but a hash table would help?
  std::vector<Filtration_value> neighbors_data;
  void init_neighbors_dense(){
    neighbors_data.clear();
    neighbors_data.resize(num_vertices*num_vertices, std::numeric_limits<Filtration_value>::infinity());
  }
  Filtration_value& neighbors_dense(Vertex i, Vertex j){return neighbors_data[num_vertices*j+i];}
#endif

  // This does not touch the events list, only the adjacency matrix(es)
  void delay_neighbor(Vertex u, Vertex v, Filtration_value f) {
    neighbors[u][v]=f;
    neighbors[v][u]=f;
#ifdef GUDHI_COLLAPSE_USE_DENSE_ARRAY
    neighbors_dense(u,v)=f;
    neighbors_dense(v,u)=f;
#endif
  }
  void remove_neighbor(Vertex u, Vertex v) {
    neighbors[u].erase(v);
    neighbors[v].erase(u);
#ifdef GUDHI_COLLAPSE_USE_DENSE_ARRAY
    neighbors_dense(u,v)=std::numeric_limits<Filtration_value>::infinity();
    neighbors_dense(v,u)=std::numeric_limits<Filtration_value>::infinity();
#endif
  }

  template<class FilteredEdgeRange>
  void read_edges(FilteredEdgeRange const&r){
    neighbors.resize(num_vertices);
#ifdef GUDHI_COLLAPSE_USE_DENSE_ARRAY
    init_neighbors_dense();
#endif
    // Use the raw sequence to avoid maintaining the order
    std::vector<typename Ngb_list::sequence_type> neighbors_seq(num_vertices);
    for(auto&&e : r){
      using std::get;
      Vertex u = get<0>(e);
      Vertex v = get<1>(e);
      Filtration_value f = get<2>(e);
      neighbors_seq[u].emplace_back(v, f);
      neighbors_seq[v].emplace_back(u, f);
#ifdef GUDHI_COLLAPSE_USE_DENSE_ARRAY
      neighbors_dense(u,v)=f;
      neighbors_dense(v,u)=f;
#endif
    }
    for(std::size_t i=0;i<neighbors_seq.size();++i){
      neighbors_seq[i].emplace_back(i, -std::numeric_limits<Filtration_value>::infinity());
      neighbors[i].adopt_sequence(std::move(neighbors_seq[i])); // calls sort
#ifdef GUDHI_COLLAPSE_USE_DENSE_ARRAY
      neighbors_dense(i,i)=-std::numeric_limits<Filtration_value>::infinity();
#endif
    }
  }

  // Open neighborhood
  // At some point it helped gcc to add __attribute__((noinline)) here, otherwise we had +50% on the running time
  // on one example. It looks ok now, or I forgot which example that was.
  void common_neighbors(boost::container::flat_set<Vertex>& e_ngb,
      std::vector<std::pair<Filtration_value, Vertex>>& e_ngb_later,
      Vertex u, Vertex v, Filtration_value f_event){
    // Using neighbors_dense here seems to hurt, even if we loop on the smaller of nu and nv.
    Ngb_list const&nu = neighbors[u];
    Ngb_list const&nv = neighbors[v];
    auto ui = nu.begin();
    auto ue = nu.end();
    auto vi = nv.begin();
    auto ve = nv.end();
    assert(ui != ue && vi != ve);
    while(ui != ue && vi != ve){
      Vertex w = ui->first;
      if(w < vi->first) { ++ui; continue; }
      if(w > vi->first) { ++vi; continue; }
      // nu and nv are closed, so we need to exclude e here.
      if(w != u && w != v) {
        Filtration_value f = std::max(ui->second, vi->second);
        if(f > f_event)
          e_ngb_later.emplace_back(f, w);
        else
          e_ngb.insert(e_ngb.end(), w);
      }
      ++ui; ++vi;
    }
  }

  // Test if the neighborhood of e is included in the closed neighborhood of c
  template<class Ngb>
  bool is_dominated_by(Ngb const& e_ngb, Vertex c, Filtration_value f){
    // The best strategy probably depends on the distribution, how sparse / dense the adjacency matrix is,
    // how (un)balanced the sizes of e_ngb and nc are.
    // Some efficient operations on sets work best with bitsets, although the need for a map complicates things.
#ifdef GUDHI_COLLAPSE_USE_DENSE_ARRAY
    for(auto v : e_ngb) {
      // if(v==c)continue;
      if(neighbors_dense(v,c) > f) return false;
    }
    return true;
#else
    auto&&nc = neighbors[c];
    // if few neighbors, use dichotomy? Seems slower.
    // I tried storing a copy of neighbors as a vector<absl::flat_hash_map> and using it for nc, but it was
    // a bit slower here. It did help with neighbors[dominator].find(w) in the main function though,
    // sometimes enough, sometimes not.
    auto ci = nc.begin();
    auto ce = nc.end();
    auto eni = e_ngb.begin();
    auto ene = e_ngb.end();
    assert(eni != ene);
    assert(ci != ce);
    // if(*eni == c && ++eni == ene) return true;
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
      // If we stored an open neighborhood of c (excluding c), we would need to test for c here and before the loop
      // if(*eni == c && ++eni == ene)return true;
      if(++ci == ce) return false;
    }
#endif
  }

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

    boost::container::flat_set<Vertex> e_ngb;
    e_ngb.reserve(num_vertices);
    std::vector<std::pair<Filtration_value, Vertex>> e_ngb_later;
    for(auto&e:edges) {
      {
        Vertex u = std::get<0>(e);
        Vertex v = std::get<1>(e);
        Filtration_value input_time = std::get<2>(e);
        auto time = delay(input_time);
        auto start_time = time;
        e_ngb.clear();
        e_ngb_later.clear();
        common_neighbors(e_ngb, e_ngb_later, u, v, time);
        // If we identify a good candidate (the first common neighbor) for being a dominator of e until infinity,
        // we could check that a bit more cheaply. It does not seem to help though.
        auto cmp1=[](auto const&a, auto const&b){return a.first > b.first;};
        auto e_ngb_later_begin=e_ngb_later.begin();
        auto e_ngb_later_end=e_ngb_later.end();
        bool heapified = false;

        bool dead = false;
        while(true) {
          Vertex dominator = -1;
          // special case for size 1
          // if(e_ngb.size()==1){dominator=*e_ngb.begin();}else
          // It is tempting to test the dominators in increasing order of filtration value, which is likely to reduce
          // the number of calls to is_dominated_by before finding a dominator, but sorting, even partially / lazily,
          // is very expensive.
          for(auto c : e_ngb){
            if(is_dominated_by(e_ngb, c, time)){
              dominator = c;
              break;
            }
          }
          if(dominator==-1) break;
          // Push as long as dominator remains a dominator.
          // Iterate on times where at least one neighbor appears.
          for (bool still_dominated = true; still_dominated; ) {
            if(e_ngb_later_begin == e_ngb_later_end) {
              dead = true; goto end_move;
            }
            if(!heapified) {
              // Eagerly sorting can be slow
              std::make_heap(e_ngb_later_begin, e_ngb_later_end, cmp1);
              heapified=true;
            }
            time = e_ngb_later_begin->first; // first place it may become critical
            // Update the neighborhood for this new time, while checking if any of the new neighbors break domination.
            while (e_ngb_later_begin != e_ngb_later_end && e_ngb_later_begin->first <= time) {
              Vertex w = e_ngb_later_begin->second;
#ifdef GUDHI_COLLAPSE_USE_DENSE_ARRAY
              if (neighbors_dense(dominator,w) > e_ngb_later_begin->first)
                still_dominated = false;
#else
              auto& ngb_dom = neighbors[dominator];
              auto wit = ngb_dom.find(w); // neighborhood may be open or closed, it does not matter
              if (wit == ngb_dom.end() || wit->second > e_ngb_later_begin->first)
                still_dominated = false;
#endif
              e_ngb.insert(w);
              std::pop_heap(e_ngb_later_begin, e_ngb_later_end--, cmp1);
            }
          } // this doesn't seem to help that much...
        }
end_move:
        if(dead) {
          remove_neighbor(u, v);
        } else if(start_time != time) {
          delay_neighbor(u, v, time);
          res.emplace_back(u, v, time);
        } else {
          res.emplace_back(u, v, input_time);
        }
      }
    }
  }

  std::vector<Filtered_edge> output() {
    return std::move(res);
  }

};

template<class R> R to_range(R&& r) { return std::move(r); }
template<class R, class T> R to_range(T const& t) { R r; r.insert(r.end(), t.begin(), t.end()); return r; }

template<class FilteredEdgeRange, class Delay>
auto flag_complex_collapse_edges(FilteredEdgeRange&& edges, Delay&&delay) {
  // Would it help to label the points according to some spatial sorting?
  auto first_edge_itr = std::begin(edges);
  using Vertex = std::decay_t<decltype(std::get<0>(*first_edge_itr))>;
  using Filtration_value = std::decay_t<decltype(std::get<2>(*first_edge_itr))>;
  using Edge_collapser = Flag_complex_edge_collapser<Vertex, Filtration_value>;
  if (first_edge_itr != std::end(edges)) {
    auto edges2 = to_range<std::vector<typename Edge_collapser::Filtered_edge>>(std::forward<FilteredEdgeRange>(edges));
#ifdef GUDHI_USE_TBB
    // I think this sorting is always negligible compared to the collapse, but parallelizing it shouldn't hurt.
    tbb::parallel_sort(edges2.begin(), edges2.end(),
        [](auto const&a, auto const&b){return std::get<2>(a)>std::get<2>(b);});
#else
    std::sort(edges2.begin(), edges2.end(), [](auto const&a, auto const&b){return std::get<2>(a)>std::get<2>(b);});
#endif
    Edge_collapser edge_collapser;
    edge_collapser.process_edges(edges2, std::forward<Delay>(delay));
    return edge_collapser.output();
  }
  return std::vector<typename Edge_collapser::Filtered_edge>();
}

/** \brief Implicitly constructs a flag complex from edges as an input, collapses edges while preserving the persistent
 * homology and returns the remaining edges as a range. The filtration value of vertices is irrelevant to this function.
 *
 * \param[in] edges Range of Filtered edges. There is no need for the range to be sorted, as it will be done internally.
 *
 * \tparam FilteredEdgeRange Range of `std::tuple<Vertex_handle, Vertex_handle, Filtration_value>`
 * where `Vertex_handle` is the type of a vertex index.
 *
 * \return Remaining edges after collapse as a range of
 * `std::tuple<Vertex_handle, Vertex_handle, Filtration_value>`.
 *
 * \ingroup edge_collapse
 *
 * \note
 * Advanced: Defining the macro GUDHI_COLLAPSE_USE_DENSE_ARRAY tells gudhi to allocate a square table of size the
 * maximum vertex index. This usually speeds up the computation for dense graphs. However, for sparse graphs, the memory
 * use may be problematic and initializing this large table may be slow.
 */
template<class FilteredEdgeRange> auto flag_complex_collapse_edges(const FilteredEdgeRange& edges) {
  return flag_complex_collapse_edges(edges, [](auto const&d){return d;});
}

}  // namespace collapse

}  // namespace Gudhi

#endif  // FLAG_COMPLEX_EDGE_COLLAPSER_H_
