/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s): Marc Glisse
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

//  ds_find_set_ is inspired from code in Boost.Graph that is
//
//  (C) Copyright Jeremy Siek 2004
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)


#ifndef PERSISTENCE_ON_RECTANGLE_H
#define PERSISTENCE_ON_RECTANGLE_H

#include <gudhi/Debug_utils.h>
#ifdef GUDHI_DETAILED_TIMES
 #include <gudhi/Clock.h>
#endif

#include <boost/range/adaptor/reversed.hpp>

#ifdef GUDHI_USE_TBB
 #include <tbb/parallel_sort.h>
#endif

#ifdef DEBUG_TRACES
 #include <iostream>
#endif
#include <vector>
#include <memory>
#include <algorithm>
#include <stdexcept>
#include <cstddef>

namespace Gudhi::cubical_complex {

// When building a cubical complex from top-dimensional cells, there are
// normally more vertices than input top cells ((x+1)*(y+1) instead of x*y).
// However, the top cells along the boundary turn out to be collapsible, and we
// only need to work with (x-1)*(y-1) vertices.
template <class Filtration_value, class Index = std::size_t, bool output_index = false>
struct Persistence_on_rectangle {
  // If we want to save space, we don't have to store the redundant 'first'
  // field in T_with_index. However, it would slow down the primal pass.
  struct T_with_index {
    Filtration_value first; Index second;
    T_with_index() = default;
    T_with_index(Filtration_value f, Index i) : first(f), second(i) {}
    bool operator<(T_with_index const& other) const {
      return std::tie(first, second) < std::tie(other.first, other.second);
    }
    Index out() const { return second; }
  };
  // Don't store the index if we don't want to output it.
  struct T_no_index {
    Filtration_value first;
    T_no_index() = default;
    T_no_index(Filtration_value f, Index) : first(f) {}
    bool operator<(T_no_index const& other) const { return first < other.first; }
    Filtration_value out() const { return first; }
  };
  typedef std::conditional_t<output_index, T_with_index, T_no_index> T;

  Filtration_value const* input_p;
  Filtration_value input(Index i) const { return input_p[i]; }

  // size_* counts the number of vertices in each direction.
  Index size_x, size_y, input_size;
  // The square i + dy is right above i.
  Index dy;

  // Squares keep their index from the input.
  // Vertices have the same index as the square at their bottom left (smaller x and y)
  // Store the filtration value of vertices that could be critical. We could store them
  // in some map, or recompute them on demand, but this strongly affects performance.
  std::unique_ptr<T[]> data_v_;
  T& data_vertex(Index i){ return data_v_[i]; }
  T data_vertex(Index i) const { return data_v_[i]; }

  std::conditional_t<output_index, Index, Filtration_value> global_min;

  // Information on a cluster
  // We do not use the rank/size heuristics, they do not go well with the pre-pairing and end up slowing things down.
  // We thus use the same representative for disjoint-sets and persistence (the minimum).
  std::unique_ptr<Index[]> ds_parent_v_;
  std::vector<Index> ds_parent_s_;
  Index& ds_parent_vertex(Index n) { return ds_parent_v_[n]; }
  Index& ds_parent_square(Index n) { return ds_parent_s_[n]; }

  template<class Parent>
  Index ds_find_set_(Index v, Parent&&ds_parent) {
    // Experimentally, path halving is currently the fastest. Note that with a
    // different algorithm, full compression was faster, so make sure to check
    // again if the algorithm changes.
    // (the setting is unusual because we start from a forest with broken ranks)
#if 0
    // Full compression
    Index old = v;
    Index ancestor = ds_parent(v);
    while (ancestor != v)
    {
      v = ancestor;
      ancestor = ds_parent(v);
    }
    v = ds_parent(old);
    while (ancestor != v)
    {
      ds_parent(old) = ancestor;
      old = v;
      v = ds_parent(old);
    }
    return ancestor;
#elif 1
    // Path halving
    Index parent = ds_parent(v);
    Index grandparent = ds_parent(parent);
    while (parent != grandparent)
    {
      ds_parent(v) = grandparent;
      v = grandparent;
      parent = ds_parent(v);
      grandparent = ds_parent(parent);
    }
    return parent;
#elif 1
    // Path splitting
    Index parent = ds_parent(v);
    Index grandparent = ds_parent(parent);
    while (parent != grandparent)
    {
      ds_parent(v) = grandparent;
      v = parent;
      parent = grandparent;
      grandparent = ds_parent(parent);
    }
    return parent;
#elif 1
    // No compression (just for reference)
    Index parent;
    while (v != (parent = ds_parent(v)))
      v = parent;
    return v;
#endif
  }
  Index ds_find_set_vertex(Index v) {
    return ds_find_set_(v, [this](Index i) -> Index& { return ds_parent_vertex(i); });
  }
  Index ds_find_set_square(Index v) {
    return ds_find_set_(v, [this](Index i) -> Index& { return ds_parent_square(i); });
  }

  struct Edge {
    T f;
    Index v1, v2; // v1 < v2
    Edge() = default;
    Edge(T f, Index v1, Index v2) : f(f), v1(v1), v2(v2) {}
    Filtration_value filt() const { return f.first; }
    bool operator<(Edge const& other) const { return filt() < other.filt(); }
  };
  void dualize_edge(Edge& e) const {
    Index new_v2 = e.v1 + (dy + 1);
    e.v1 = e.v2;
    e.v2 = new_v2;
  };
  std::vector<Edge> edges;

  void init(const Filtration_value* input_, Index n_rows, Index n_cols) {
    input_size = n_rows * n_cols;
    input_p = input_;
#ifdef DEBUG_TRACES
    std::clog << "Input\n";
    for(Index i = 0; i < input_size; ++i) {
      std::clog << i << '\t' << input(i) << '\n';
    }
#endif
    dy = n_cols;
    size_x = dy - 1;
    size_y = n_rows - 1;
    // The unique_ptr could be std::vector, but the initialization is useless.
    data_v_.reset(new T[input_size - dy - 1]); // 1 row/column less for vertices than squares
    ds_parent_v_.reset(new Index[input_size - dy - 1]);
    // Initializing the boundary squares to 0 is important, it represents the infinite exterior cell.
    ds_parent_s_.resize(input_size);
    // What is a good estimate here? For a random 1000x1000 input, we get ~311k edges. For a checkerboard, ~498k.
    edges.reserve(input_size / 2);
  }

  bool has_larger_input(Index a, Index b, Filtration_value fb) const {
    // Is passing fb useful, or would the compiler notice that it already has it available?
    GUDHI_CHECK(a != b, std::logic_error("Bug in Gudhi: comparing a cell to itself"));
    Filtration_value fa = input(a);
    if (fb < fa) return true;
    if (fa < fb) return false;
    return a > b; // Arbitrary, but has to be consistent
  }
  void set_parent_vertex(Index child, Index parent) {
    GUDHI_CHECK(child != parent, std::logic_error("Bug in Gudhi: use mark_*_critical instead of set_parent"));
    ds_parent_vertex(child) = parent;
  }
  void set_parent_square(Index child, Index parent) {
    GUDHI_CHECK(child != parent, std::logic_error("Bug in Gudhi: use mark_*_critical instead of set_parent"));
    ds_parent_square(child) = parent;
  }

  // Locally pair simplices around each square.
  // Work implicitly from input, only store the filtration value of critical vertices (squares are already in input).
  // Store critical edges for later processing.
  void fill_and_pair() {
    Index i; // Index of the current square
    Filtration_value f; // input(i)
    auto mark_vertex_critical = [&](Index c) {
      ds_parent_vertex(c) = c;
      data_vertex(c) = T(f, i);
    };
    auto mark_square_critical = [&]() {
      ds_parent_square(i) = i;
    };
    auto mark_edge_critical = [&](Index v1, Index v2) {
      edges.emplace_back(T(f, i), v1, v2);
    };
    auto v_up_left    = [&](){ return i - 1; };
    auto v_up_right   = [&](){ return i; };
    auto v_down_left  = [&](){ return i - dy - 1; };
    auto v_down_right = [&](){ return i - dy; };
    auto pair_square_up    = [&](){ set_parent_square(i, i + dy); };
    auto pair_square_down  = [&](){ set_parent_square(i, i - dy); };
    auto pair_square_left  = [&](){ set_parent_square(i, i - 1); };
    auto pair_square_right = [&](){ set_parent_square(i, i + 1); };

    // Mark the corners as critical, it will be overwritten if not
    i = 0; f = input(i);
    mark_vertex_critical(v_up_right());
    i = size_x; f = input(i);
    mark_vertex_critical(v_up_left());
    i = dy * size_y; f = input(i);
    mark_vertex_critical(v_down_right());
    i = size_x + dy * size_y; f = input(i);
    mark_vertex_critical(v_down_left());

    // Boundary nodes, 1st row
    for(Index x = 1; x < size_x; ++x) {
      i = x;
      f = input(x);
      if (has_larger_input(i + dy, i, f)) {
        auto up_left  = [&](){ return has_larger_input(i - 1, i, f) && has_larger_input(i + dy - 1, i, f); };
        auto up_right = [&](){ return has_larger_input(i + 1, i, f) && has_larger_input(i + dy + 1, i, f); };
        if (up_left()) {
          set_parent_vertex(v_up_left(), v_up_right());
          if (up_right()) mark_vertex_critical(v_up_right());
        } else if (up_right()) {
          set_parent_vertex(v_up_right(), v_up_left());
        } else {
          mark_edge_critical(v_up_left(), v_up_right());
        }
      }
    }
    // Internal rows
    for(Index y = 1; y < size_y; ++y) {
      // First column
      {
        i = y * dy;
        f = input(i);
        if (has_larger_input(i + 1, i, f)) {
          auto down_right = [&](){ return has_larger_input(i - dy, i, f) && has_larger_input(i + 1 - dy, i, f); };
          auto up_right   = [&](){ return has_larger_input(i + dy, i, f) && has_larger_input(i + 1 + dy, i, f); };
          if (down_right()) {
            set_parent_vertex(v_down_right(), v_up_right());
            if (up_right()) mark_vertex_critical(v_up_right());
          } else if (up_right()) {
            set_parent_vertex(v_up_right(), v_down_right());
          } else {
            mark_edge_critical(v_down_right(), v_up_right());
          }
        }
      }
      // Internal squares
      for(Index x = 1; x < size_x; ++x) {
        i = x + dy * y;
        f = input(i);
        // See what part of the boundary shares f
        auto left  = [&]() { return has_larger_input(i - 1, i, f); };
        auto right = [&]() { return has_larger_input(i + 1, i, f); };
        auto down  = [&]() { return has_larger_input(i - dy, i, f); };
        auto up    = [&]() { return has_larger_input(i + dy, i, f); };
        auto down_left  = [&]() { return has_larger_input(i - dy - 1, i, f); };
        auto up_left    = [&]() { return has_larger_input(i + dy - 1, i, f); };
        auto down_right = [&]() { return has_larger_input(i - dy + 1, i, f); };
        auto up_right   = [&]() { return has_larger_input(i + dy + 1, i, f); };
        if (up()) { // u
          if (left()) { // u l
            if (up_left()) { // u l ul
              set_parent_vertex(v_up_left(), v_up_right());
              if (down()) { // U l UL d
                if (down_left()) { // U l UL d dl
                  set_parent_vertex(v_down_left(), v_up_left());
                  if (right()) { // U L UL d DL r
                    if (down_right()) { // U L UL d DL r dr
                      set_parent_vertex(v_down_right(), v_down_left());
                      pair_square_right();
                      if (up_right()) { // U L UL D DL R DR ur - cr
                        mark_vertex_critical(v_up_right());
                      }
                    } else { // U L UL d DL r !dr
                      pair_square_down();
                      if (up_right()) { // U L UL D DL r !dr ur - cd
                        set_parent_vertex(v_up_right(), v_down_right());
                      } else { // U L UL D DL r !dr !ur - cd
                        mark_edge_critical(v_down_right(), v_up_right());
                      }
                    }
                  } else { // U L UL d DL !r
                    pair_square_down();
                  }
                } else { // U l UL d !dl
                  pair_square_left();
                  if (right()) { // U L UL d !dl r - cl
                    if (down_right()) { // U L UL d !dl r dr - cl
                      set_parent_vertex(v_down_right(), v_down_left());
                    } else { // U L UL d !dl r !dr - cl
                      mark_edge_critical(v_down_left(), v_down_right());
                    }
                    if (up_right()) { // U L UL D !dl r ur - cl
                      set_parent_vertex(v_up_right(), v_down_right());
                    } else { // U L UL D !dl r !ur - cl
                      mark_edge_critical(v_down_right(), v_up_right());
                    }
                  } else { // U L UL d !dl !r - cl
                    mark_edge_critical(v_down_left(), v_down_right());
                  }
                }
              } else { // U l UL !d
                pair_square_left();
                if (right()) { // U L UL !d r - cl
                  if (up_right()) { // U L UL !d r ur - cl
                    set_parent_vertex(v_up_right(), v_down_right());
                  } else { // U L UL !d r !ur - cl
                    mark_edge_critical(v_down_right(), v_up_right());
                  }
                } else {} // U L UL !d !r - cl
              }
            } else { // u l !ul
              pair_square_up();
              if (down()) { // U l !ul d - cu
                if (down_left()) { // U l !ul d dl - cu
                  set_parent_vertex(v_down_left(), v_up_left());
                } else { // U l !ul d !dl - cu
                  mark_edge_critical(v_down_left(), v_up_left());
                }
                if (right()) { // U L !ul d r - cu
                  if (down_right()) { // U L !ul d r dr - cu
                    set_parent_vertex(v_down_right(), v_down_left());
                  } else { // U L !ul d r !dr - cu
                    mark_edge_critical(v_down_left(), v_down_right());
                  }
                  if (up_right()) { // U L !ul D r ur - cu
                    set_parent_vertex(v_up_right(), v_down_right());
                  } else { // U L !ul D r !ur - cu
                    mark_edge_critical(v_down_right(), v_up_right());
                  }
                } else { // U L !ul d !r - cu
                  mark_edge_critical(v_down_left(), v_down_right());
                }
              } else { // U l !ul !d - cu
                mark_edge_critical(v_down_left(), v_up_left());
                if (right()) { // U L !ul !d r - cu
                  if (up_right()) { // U L !ul !d r ur - cu
                    set_parent_vertex(v_up_right(), v_down_right());
                  } else { // U L !ul !d r !ur - cu
                    mark_edge_critical(v_down_right(), v_up_right());
                  }
                } else {} // U L !ul !d !r - cu
              }
            }
          } else { // u !l
            pair_square_up();
            if (down()) { // U !l d - cu
              if (right()) { // U !l d r - cu
                if (down_right()) { // U !l d r dr - cu
                  set_parent_vertex(v_down_right(), v_down_left());
                } else { // U !l d r !dr - cu
                  mark_edge_critical(v_down_left(), v_down_right());
                }
                if (up_right()) { // U !l D r ur - cu
                  set_parent_vertex(v_up_right(), v_down_right());
                } else { // U !l D r !ur - cu
                  mark_edge_critical(v_down_right(), v_up_right());
                }
              } else { // U !l d !r - cu
                mark_edge_critical(v_down_left(), v_down_right());
              }
            } else { // U !l !d - cu
              if (right()) { // U !l !d r - cu
                if (up_right()) { // U !l !d r ur - cu
                  set_parent_vertex(v_up_right(), v_down_right());
                } else { // U !l !d r !ur - cu
                  mark_edge_critical(v_down_right(), v_up_right());
                }
              } else {} // U !l !d !r - cu
            }
          }
        } else { // !u
          if (left()) { // !u l
            if (down()) { // !u l d
              if (down_left()) { // !u l d dl
                set_parent_vertex(v_down_left(), v_up_left());
                if (right()) { // !u L d DL r
                  if (down_right()) { // !u L d DL r dr
                    set_parent_vertex(v_down_right(), v_down_left());
                  } else { // !u L d DL r !dr
                    mark_edge_critical(v_down_left(), v_down_right());
                  }
                  pair_square_right();
                } else { // !u L d DL !r
                  pair_square_down();
                }
              } else { // !u l d !dl
                pair_square_left();
                if (right()) { // !u L d !dl r - cl
                  if (down_right()) { // !u L d !dl r dr - cl
                    set_parent_vertex(v_down_right(), v_down_left());
                  } else { // !u L d !dl r !dr - cl
                    mark_edge_critical(v_down_left(), v_down_right());
                  }
                  mark_edge_critical(v_down_right(), v_up_right());
                } else { // !u L d !dl !r - cl
                  mark_edge_critical(v_down_left(), v_down_right());
                }
              }
            } else { // !u l !d
              pair_square_left();
              if (right()) { // !u L !d r - cl
                mark_edge_critical(v_down_right(), v_up_right());
              } else {} // !u L !d !r - cl
            }
          } else { // !u !l
            if (down()) { // !u !l d
              pair_square_down();
              if (right()) { // !u !l D r - cd
                if (down_right()) { // !u !l D r dr - cd
                  set_parent_vertex(v_down_right(), v_up_right());
                } else { // !u !l D r !dr - cd
                  mark_edge_critical(v_down_right(), v_up_right());
                }
              } else {} // !u !l D !r - cd
            } else { // !u !l !d
              if (right()) { // !u !l !d r
                pair_square_right();
              } else { // !u !l !d !r
                mark_square_critical();
              }
            }
          }
        }
      }
      // Last column
      {
        i = size_x + dy * y;
        f = input(i);
        if (has_larger_input(i - 1, i, f)) {
          auto down_left = [&](){ return has_larger_input(i - dy, i, f) && has_larger_input(i - 1 - dy, i, f); };
          auto up_left   = [&](){ return has_larger_input(i + dy, i, f) && has_larger_input(i - 1 + dy, i, f); };
          if (down_left()) {
            set_parent_vertex(v_down_left(), v_up_left());
            if (up_left()) mark_vertex_critical(v_up_left());
          } else if (up_left()) {
            set_parent_vertex(v_up_left(), v_down_left());
          } else {
            mark_edge_critical(v_down_left(), v_up_left());
          }
        }
      }
    }
    // Boundary nodes, last row
    for(Index x = 1; x < size_x; ++x) {
      i = size_y * dy + x;
      f = input(i);
      if (has_larger_input(i - dy, i, f)) {
        auto down_left  = [&](){ return has_larger_input(i - 1, i, f) && has_larger_input(i - dy - 1, i, f); };
        auto down_right = [&](){ return has_larger_input(i + 1, i, f) && has_larger_input(i - dy + 1, i, f); };
        if (down_left()) {
          set_parent_vertex(v_down_left(), v_down_right());
          if (down_right()) mark_vertex_critical(v_down_right());
        } else if (down_right()) {
          set_parent_vertex(v_down_right(), v_down_left());
        } else {
          mark_edge_critical(v_down_left(), v_down_right());
        }
      }
    }
  }

  void sort_edges(){
#ifdef GUDHI_USE_TBB
    // Parallelizing just this part is a joke. It would be possible to
    // parallelize the pairing (one edge list per thread) and run the dual in
    // parallel with the primal if we were motivated...
    tbb::parallel_sort(edges.begin(), edges.end());
#else
    std::sort(edges.begin(), edges.end());
#endif
#ifdef DEBUG_TRACES
    std::clog << "edges\n";
    for(auto&e : edges){ std::clog << e.v1 << '\t' << e.v2 << '\t' << e.filt() << '\n'; }
#endif
  }

  template<class Out>
  void primal(Out&&out){
    auto it = std::remove_if(edges.begin(), edges.end(), [&](Edge& e) {
        assert(e.v1 < e.v2);
        Index a = ds_find_set_vertex(e.v1);
        Index b = ds_find_set_vertex(e.v2);
        if (a == b) return false;
        if (data_vertex(b) < data_vertex(a)) std::swap(a, b);
        ds_parent_vertex(b) = a;
        out(data_vertex(b).out(), e.f.out());
        return true;
    });
    edges.erase(it, edges.end());
    global_min = data_vertex(ds_find_set_vertex(0)).out();
  }

  // In the dual, squares behave like vertices, and edges are rotated 90° around their middle.
  // To handle boundaries correctly, we imagine a single exterior cell with filtration +inf.
  template<class Out>
  void dual(Out&&out){
    for (auto e : boost::adaptors::reverse(edges)) {
      dualize_edge(e);
      Index a = ds_find_set_square(e.v1);
      Index b = ds_find_set_square(e.v2);
      GUDHI_CHECK(a != b, std::logic_error("Bug in Gudhi"));
      // This is more robust in case the input contains inf? I used to set the filtration of 0 to inf.
      if (b == 0 || (a != 0 && input(a) < input(b))) std::swap(a, b);
      ds_parent_square(b) = a;
      if constexpr (output_index)
        out(e.f.out(), b);
      else
        out(e.f.out(), input(b));
    }
  }
};
// Ideas for improvement:
// * for large hard (many intervals) inputs, primal/dual dominate the running time because of the random reads in
//   find_set and input(a/b). The input load would be cheaper if we stored it with parents, but then find_set would
//   be slower.
// * to increase memory locality, maybe pairing, which is currently arbitrary, could use some heuristic to favor
//   some pairs over others.
// * try to loosen tight dependency chains, load values several instructions before they are needed. Performing 2
//   find_set in lock step surprisingly doesn't help.
// * To handle very big instances, we could remove data_v_ and recompute it on demand as the min of the inputs of i,
//   i+1, i+dy and i+dy+1. On hard instances, it wouldn't save that much memory (and it is a bit slower). On easy
//   instances, if we also remove the calls to reserve(), the saving is less negligible, but we still have ds_parent_*_
//   that take about as much space as the input. We could, on a subarray, fill a dense ds_parent, then reduce it and
//   export only the critical vertices and boundary to some sparse datastructure, but it doesn't seem worth the trouble
//   for now.

/**
 * @private
 * Compute the persistence diagram of a function on a 2d cubical complex, defined as a lower-star filtration of the
 * values at the top-dimensional cells.
 *
 * @tparam output_index If false, each argument of the out functors is a filtration value. If true, it is instead the
 *   index of this filtration value in the input.
 * @tparam Filtration_value Must be comparable with `operator<`.
 * @tparam Index This is used to index the elements of `input`, so it must be large enough to represent the size
 *   of `input`.
 * @param[in] input Pointer to `n_rows*n_cols` filtration values for the square cells. Note that the values are assumed
 *   to be stored in C order, unlike `Gudhi::cubical_complex::Bitmap_cubical_complex` (you can exchange `n_rows` and
 *   `n_cols` for compatibility).
 * @param[in] n_rows number of rows of `input`.
 * @param[in] n_cols number of columns of `input`.
 * @param[out] out0 For each interval (b, d) in the persistence diagram of dimension 0, the function calls `out0(b, d)`.
 * @param[out] out1 Same as `out0` for persistence in dimension 1.
 * @returns The global minimum, which is not paired and is thus the birth of an infinite persistence interval of
 *   dimension 0.
 */
template <bool output_index = false, typename Filtration_value, typename Index, typename Out0, typename Out1>
auto persistence_on_rectangle_from_top_cells(Filtration_value const* input, Index n_rows, Index n_cols,
                                             Out0&&out0, Out1&&out1){
#ifdef GUDHI_DETAILED_TIMES
  Gudhi::Clock clock;
#endif
  GUDHI_CHECK(n_rows >= 2 && n_cols >= 2,
      std::domain_error("The complex must truly be 2d, i.e. at least 2 rows and 2 columns"));
  Persistence_on_rectangle<Filtration_value, Index, output_index> X;
  X.init(input, n_rows, n_cols);
#ifdef GUDHI_DETAILED_TIMES
    std::clog << "init: " << clock; clock.begin();
#endif
  X.fill_and_pair();
#ifdef GUDHI_DETAILED_TIMES
    std::clog << "fill and pair: " << clock; clock.begin();
#endif
  X.sort_edges();
#ifdef GUDHI_DETAILED_TIMES
    std::clog << "sort: " << clock; clock.begin();
#endif
  X.primal(out0);
#ifdef GUDHI_DETAILED_TIMES
    std::clog << "primal pass: " << clock; clock.begin();
#endif
  X.dual(out1);
#ifdef GUDHI_DETAILED_TIMES
    std::clog << "dual pass: " << clock;
#endif
  return X.global_min;
}
}  // namespace Gudhi::cubical_complex

#endif  // PERSISTENCE_ON_RECTANGLE_H
