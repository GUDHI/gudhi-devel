/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s): Marc Glisse (with Claude Sonnet 5)
 *
 *    Copyright (C) 2026 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

//  Adapted from Persistence_on_rectangle.h, which computes persistence for
//  the T-construction (values given at top-dimensional cells, extended to
//  lower-dimensional cells by the MINIMUM of incident top cells). This file
//  is its V-construction counterpart: values are given at vertices instead,
//  and extended to edges and squares by the MAXIMUM of incident vertices.
//  See Persistence_on_rectangle.h for details.

#ifndef GUDHI_PERSISTENCE_ON_RECTANGLE_V_H
#define GUDHI_PERSISTENCE_ON_RECTANGLE_V_H

#include <gudhi/Debug_utils.h>
#ifdef GUDHI_DETAILED_TIMES
 #include <gudhi/Clock.h>
 #include <iostream>
#endif

#include <boost/config.hpp>
#include <boost/range/adaptor/reversed.hpp>

#ifdef GUDHI_USE_TBB
 #include <tbb/parallel_sort.h>
#else
 #include <boost/sort/pdqsort/pdqsort.hpp>
#endif

#ifdef DEBUG_TRACES
 #include <iostream>
#endif
#include <vector>
#include <memory>
#include <algorithm>
#include <stdexcept>
#include <cstddef>
#include <cstdint>
#include <type_traits>

namespace Gudhi::cubical_complex {

template <class Filtration_value, class Index = std::size_t, bool output_index = false>
struct Persistence_on_rectangle_V {
  struct T_with_index {
    Filtration_value first; Index second;
    T_with_index() = default;
    T_with_index(Filtration_value f, Index i) : first(f), second(i) {}
    Index out() const { return second; }
  };
  struct T_no_index {
    Filtration_value first;
    T_no_index() = default;
    T_no_index(Filtration_value f, Index) : first(f) {}
    Filtration_value out() const { return first; }
  };
  typedef std::conditional_t<output_index, T_with_index, T_no_index> T;

  Filtration_value const* input_p;
  Filtration_value input(Index i) const { return input_p[i]; }

  // size_x/size_y count VERTICES in each direction (the full input grid,
  // nothing is dropped). dy is the row stride for the vertex grid.
  Index size_x, size_y, input_size;
  Index dy;

  // Squares are indexed by their own down-right corner, sharing the same
  // stride (dy) and index range as the vertices themselves: the four squares
  // around vertex i are then simply i, i+1, i+dy, i+dy+1.  This leaves the
  // entire top row and left column of this index range (including 0) unused by
  // any real square so they can safely represent the exterior cell.
  Index exterior = 0;

  // Birth value of each derived square = max of its 4 corners.
  // This is only needed for critical squares, which from fill_and_pair
  // always have the value of the same corner.
  T square_birth(Index s) const { return T(input(s - 1), s - 1); }

  // Union-find forests.
  std::unique_ptr<Index[]> ds_parent_vertex_;
  std::unique_ptr<Index[]> ds_parent_square_;
  Index& ds_parent_vertex(Index n) { return ds_parent_vertex_[n]; }
  Index& ds_parent_square(Index n) { return ds_parent_square_[n]; }

  // Birth of the single infinite H0 interval: the vertex's own index
  // when output_index, its raw value otherwise -- set in primal().
  std::conditional_t<output_index, Index, Filtration_value> global_min;

  template<class Parent>
  Index ds_find_set_(Index v, Parent&& ds_parent) {
    // Path halving
    Index parent = ds_parent(v);
    Index grandparent = ds_parent(parent);
    while (parent != grandparent) {
      ds_parent(v) = grandparent;
      v = grandparent;
      parent = ds_parent(v);
      grandparent = ds_parent(parent);
    }
    return parent;
  }
  Index ds_find_set_vertex(Index v) {
    return ds_find_set_(v, [this](Index i) -> Index& { return ds_parent_vertex(i); });
  }
  Index ds_find_set_square(Index v) {
    return ds_find_set_(v, [this](Index i) -> Index& { return ds_parent_square(i); });
  }

  void set_parent_vertex(Index child, Index parent) {
    GUDHI_CHECK(child != parent, std::logic_error("Bug: use a direct self-assignment to mark a vertex critical"));
    ds_parent_vertex(child) = parent;
  }
  void set_parent_square(Index child, Index parent) {
    ds_parent_square(child) = parent;
  }

  struct Edge {
    T f;  // derived (max-of-two-endpoints) value
    Index v1, v2;  // primal (real vertex) endpoints, v1 < v2
    Edge() = default;
    Edge(T f, Index v1, Index v2) : f(f), v1(v1), v2(v2) {}
    bool operator<(Edge const& other) const { return f.first < other.f.first; }
  };
  std::vector<Edge> edges;

  // This works thanks to the empty row and column of 0s from init.
  std::pair<Index, Index> dualize_edge(Index v1, Index v2) const {
    return { v2, v1 + dy + 1 };
  }

  void init(Filtration_value const* input_, Index n_rows, Index n_cols) {
    input_p = input_;
    size_x = n_cols; size_y = n_rows;
    dy = n_cols;
    input_size = n_rows * n_cols;

    ds_parent_vertex_.reset(new Index[input_size]);
    // Only `exterior` needs pre-initializing. In addition, we surround the
    // real squares with 0 so dualize_edge does not need to check boundary
    // cases.
    ds_parent_square_.reset(new Index[input_size+dy]);
    // ds_parent_square_[exterior] = exterior;
    for (Index i = 0; i < dy; ++i) ds_parent_square_[i] = exterior;
    for (Index i = 0; i < dy; ++i) ds_parent_square_[input_size + i] = exterior;
    for (Index i = 1; i < size_y; ++i) ds_parent_square_[dy * i] = exterior;

    edges.reserve(input_size / 2);  // same rough order-of-magnitude estimate as T's
  }

  bool beats(Index a, Filtration_value fa, Index b) const {
    // True if vertex a is functionally LARGER than vertex b.
    GUDHI_CHECK(a != b, std::logic_error("Bug: comparing a vertex to itself"));
    Filtration_value fb = input(b);
    if (fa > fb) return true;
    if (fa < fb) return false;
    return a > b;  // tie-break: larger index wins
  }
  // Same as beats when we already know the order of a and b.
  bool beats_before(Index a, Filtration_value fa, Index b) const {
    GUDHI_CHECK(a > b, std::logic_error("Bug: inconsistent order"));
    Filtration_value fb = input(b);
    return fa >= fb;
  }
  bool beats_after(Index a, Filtration_value fa, Index b) const {
    GUDHI_CHECK(a < b, std::logic_error("Bug: inconsistent order"));
    Filtration_value fb = input(b);
    return fa > fb;
  }

  BOOST_FORCEINLINE
  void fill_and_pair_interior(Index y, Index x) {
    Index i = y * dy + x;
    Filtration_value f = input(i);
    auto U = [&](){ return beats_before(i, f, i - dy); };
    auto D = [&](){ return beats_after (i, f, i + dy); };
    auto L = [&](){ return beats_before(i, f, i - 1); };
    auto R = [&](){ return beats_after (i, f, i + 1); };

    auto diagUL = [&](){ return beats_before(i, f, i - dy - 1); };
    auto diagUR = [&](){ return beats_before(i, f, i - dy + 1); };
    auto diagDL = [&](){ return beats_after (i, f, i + dy - 1); };
    auto diagDR = [&](){ return beats_after (i, f, i + dy + 1); };

    Index sUL = i, sUR = i + 1, sDL = i + dy, sDR = i + dy + 1;

    if (R()) {
      if (D()) {
        if (diagDR()) {
          set_parent_square(sDR, sDL);
          if (L()) {
            if (diagDL()) {
              set_parent_square(sDL, sUL);
              if (U()) {
                if (diagUL()) {
                  set_parent_square(sUL, sUR);
                  if (diagUR()) {
                    ds_parent_square(sUR) = sUR; // square_birth(sUR) = T(f, i);
                  }
                  // The following pair exists, but nothing will look at it
                  // set_parent_vertex(i, i + 1);
                } else {
                  if (diagUR()) {
                    set_parent_square(sUR, sDR);
                    // set_parent_vertex(i, i - dy);
                  } else {
                    set_parent_vertex(i, i + 1);
                    edges.emplace_back(T(f, i), i - dy, i);
                  }
                }
              } else {
                set_parent_vertex(i, i + 1);
              }
            } else {
              if (U()) {
                if (diagUL()) {
                  set_parent_square(sUL, sUR);
                } else {
                  edges.emplace_back(T(f, i), i - dy, i);
                }
                if (diagUR()) {
                  set_parent_square(sUR, sDR);
                  set_parent_vertex(i, i - 1);
                } else {
                  set_parent_vertex(i, i + 1);
                  edges.emplace_back(T(f, i), i - 1, i);
                }
              } else {
                set_parent_vertex(i, i + 1);
                edges.emplace_back(T(f, i), i - 1, i);
              }
            }
          } else {
            if (U()) {
              if (diagUR()) {
                set_parent_square(sUR, sDR);
                set_parent_vertex(i, i - dy);
              } else {
                set_parent_vertex(i, i + 1);
                edges.emplace_back(T(f, i), i - dy, i);
              }
            } else {
              set_parent_vertex(i, i + 1);
            }
          }
        } else {
          if (L()) {
            if (diagDL()) {
              set_parent_square(sDL, sUL);
            } else {
              edges.emplace_back(T(f, i), i - 1, i);
            }
            if (U()) {
              if (diagUL()) {
                set_parent_square(sUL, sUR);
              } else {
                edges.emplace_back(T(f, i), i - dy, i);
              }
              if (diagUR()) {
                set_parent_square(sUR, sDR);
                set_parent_vertex(i, i + dy);
              } else {
                set_parent_vertex(i, i + 1);
                edges.emplace_back(T(f, i), i, i + dy);
              }
            } else {
              set_parent_vertex(i, i + 1);
              edges.emplace_back(T(f, i), i, i + dy);
            }
          } else {
            if (U()) {
              if (diagUR()) {
                set_parent_square(sUR, sDR);
                set_parent_vertex(i, i + dy);
              } else {
                set_parent_vertex(i, i + 1);
                edges.emplace_back(T(f, i), i, i + dy);
              }
              edges.emplace_back(T(f, i), i - dy, i);
            } else {
              set_parent_vertex(i, i + 1);
              edges.emplace_back(T(f, i), i, i + dy);
            }
          }
        }
      } else {
        if (L()) {
          if (U()) {
            if (diagUL()) {
              set_parent_square(sUL, sUR);
            } else {
              edges.emplace_back(T(f, i), i - dy, i);
            }
            if (diagUR()) {
              set_parent_square(sUR, sDR);
              set_parent_vertex(i, i - 1);
            } else {
              set_parent_vertex(i, i + 1);
              edges.emplace_back(T(f, i), i - 1, i);
            }
          } else {
            set_parent_vertex(i, i + 1);
            edges.emplace_back(T(f, i), i - 1, i);
          }
        } else {
          if (U()) {
            if (diagUR()) {
              set_parent_square(sUR, sDR);
              set_parent_vertex(i, i - dy);
            } else {
              set_parent_vertex(i, i + 1);
              edges.emplace_back(T(f, i), i - dy, i);
            }
          } else {
            set_parent_vertex(i, i + 1);
          }
        }
      }
    } else {
      if (D()) {
        if (L()) {
          if (diagDL()) {
            set_parent_square(sDL, sUL);
          } else {
            edges.emplace_back(T(f, i), i - 1, i);
          }
          if (U()) {
            if (diagUL()) {
              set_parent_square(sUL, sUR);
            } else {
              edges.emplace_back(T(f, i), i - dy, i);
            }
          } else {
          }
        } else {
          if (U()) {
            edges.emplace_back(T(f, i), i - dy, i);
          } else {
          }
        }
        set_parent_vertex(i, i + dy);
      } else {
        if (L()) {
          if (U()) {
            if (diagUL()) {
              set_parent_square(sUL, sUR);
            } else {
              edges.emplace_back(T(f, i), i - dy, i);
            }
          } else {
          }
          set_parent_vertex(i, i - 1);
        } else {
          if (U()) {
            set_parent_vertex(i, i - dy);
          } else {
            ds_parent_vertex(i) = i;  // nothing matched at all
          }
        }
      }
    }
  }

  // Version for boundary vertices where some neighbors don't exist.
  void fill_and_pair_boundary(Index y, Index x) {
    Index i = y * dy + x;
    Filtration_value f = input(i);
    bool has_u = y > 0, has_d = y < size_y - 1, has_l = x > 0, has_r = x < size_x - 1;
    bool U = has_u && beats(i, f, i - dy);
    bool D = has_d && beats(i, f, i + dy);
    bool L = has_l && beats(i, f, i - 1);
    bool R = has_r && beats(i, f, i + 1);

    bool qUL = U && L && beats(i, f, i - dy - 1);
    bool qUR = U && R && beats(i, f, i - dy + 1);
    bool qDL = D && L && beats(i, f, i + dy - 1);
    bool qDR = D && R && beats(i, f, i + dy + 1);
    GUDHI_CHECK(!(qUL && qUR && qDL && qDR),
                std::logic_error("Bug in Persistence_on_rectangle_V: all4 on a boundary vertex"));

    uint8_t consumed = 0;
    if (qUL || qUR || qDL || qDR) {
      bool has_ul = has_u && has_l, has_ur = has_u && has_r;
      bool has_dl = has_d && has_l, has_dr = has_d && has_r;
      Index sUL = i, sUR = i + 1, sDL = i + dy, sDR = i + dy + 1;
      if (qDR) { set_parent_square(sDR, has_dl ? sDL : exterior); consumed |= 0x02; }
      if (qDL) { set_parent_square(sDL, has_ul ? sUL : exterior); consumed |= 0x04; }
      if (qUL) { set_parent_square(sUL, has_ur ? sUR : exterior); consumed |= 0x08; }
      if (qUR) { set_parent_square(sUR, has_dr ? sDR : exterior); consumed |= 0x01; }
    }

    uint8_t matched = (R ? 0x01 : 0) | (D ? 0x02 : 0) | (L ? 0x04 : 0) | (U ? 0x08 : 0);
    uint8_t leftover = matched & static_cast<uint8_t>(~consumed);

    if (leftover != 0) {
      uint8_t d = 0;
      while (!(leftover & (1u << d))) ++d;
      set_parent_vertex(i, (d == 0) ? i + 1 : (d == 1) ? i + dy : (d == 2) ? i - 1 : i - dy);
      leftover &= static_cast<uint8_t>(~(1u << d));
    } else {
      ds_parent_vertex(i) = i;
    }
    if (leftover & 0x01) edges.emplace_back(T(f, i), i, i + 1);
    if (leftover & 0x02) edges.emplace_back(T(f, i), i, i + dy);
    if (leftover & 0x04) edges.emplace_back(T(f, i), i - 1, i);
    if (leftover & 0x08) edges.emplace_back(T(f, i), i - dy, i);
  }

  void fill_and_pair() {
    // Boundary ring (top row, then per interior row the two side cells
    // around the interior fast path, then the bottom row), in the same
    // row-major order as a single unified loop would visit them.
    for (Index x = 0; x < size_x; ++x) fill_and_pair_boundary(0, x);
    for (Index y = 1; y < size_y - 1; ++y) {
      fill_and_pair_boundary(y, 0);
      for (Index x = 1; x < size_x - 1; ++x) fill_and_pair_interior(y, x);
      fill_and_pair_boundary(y, size_x - 1);
    }
    for (Index x = 0; x < size_x; ++x) fill_and_pair_boundary(size_y - 1, x);
  }

  void sort_edges() {
#ifdef GUDHI_USE_TBB
    tbb::parallel_sort(edges.begin(), edges.end());
#else
    // std::sort(edges.begin(), edges.end());
    boost::sort::pdqsort_branchless(edges.begin(), edges.end());
#endif
  }

#ifdef __GNUC__
#define GUDHI_NOP asm(""::);
#else
#define GUDHI_NOP
#endif

  template<class Out>
  void primal(Out&& out) {
    auto it = std::remove_if(edges.begin(), edges.end(), [&](Edge& e) {
      Index a = ds_find_set_vertex(e.v1);
      Index b = ds_find_set_vertex(e.v2);
      if (a == b) return false;  // not needed by primal; dual will need it -> keep
      if (input(a) > input(b)) { GUDHI_NOP std::swap(a, b); }
      ds_parent_vertex(b) = a;
      if constexpr (output_index) out(b, e.f.out());
      else out(input(b), e.f.out());
      return true;  // used by primal -> remove from the list
    });
    edges.erase(it, edges.end());
    if constexpr (output_index) global_min = ds_find_set_vertex(0);
    else global_min = input(ds_find_set_vertex(0));
  }

  template<class Out>
  void dual(Out&& out) {
    for (auto& e : boost::adaptors::reverse(edges)) {
      auto [dv1, dv2] = dualize_edge(e.v1, e.v2);
      Index a = ds_find_set_square(dv1), b = ds_find_set_square(dv2);
      // fill_and_pair's local matching is maximal so dual() never actually
      // sees an edge whose two sides are already connected.
      GUDHI_CHECK(a != b, std::logic_error("Bug in Persistence_on_rectangle_V"));
      // a should end up as the survivor: exterior, or the larger birth.
      if(a != exterior && (b == exterior || square_birth(a).first < square_birth(b).first))
      { GUDHI_NOP std::swap(a, b);}
      ds_parent_square(b) = a;
      out(e.f.out(), square_birth(b).out());  // b is never exterior here
    }
  }
};

#undef GUDHI_NOP

/**
 * @private
 * Compute the persistence diagram of a function on a 2d cubical complex,
 * defined as an upper-star (V-construction) filtration: values are given
 * at the vertices, and each higher-dimensional cell (edge, square) takes
 * the MAXIMUM of the values of its incident vertices.
 *
 * @tparam output_index If false, each argument of the out functors is a
 *   filtration value. If true, it is instead the index (into `input`) of
 *   the vertex that achieves that value -- useful e.g. for autodiff,
 *   where the caller wants to trace a birth/death back to the specific
 *   input entry it came from rather than just its (already-detached)
 *   value.
 * @param[in] input Pointer to n_rows*n_cols filtration values, one per
 *   vertex, stored in C order.
 * @param[in] n_rows, n_cols grid dimensions (number of vertices).
 * @param[out] out0 called as out0(birth, death) for each finite interval
 *   of dimension 0.
 * @param[out] out1 called as out1(birth, death) for each interval of
 *   dimension 1
 * @returns The global minimum, i.e. the birth of the single infinite
 *   interval of dimension 0 (an index if output_index, else a value).
 */
template <bool output_index = false, typename Filtration_value, typename Index, typename Out0, typename Out1>
auto persistence_on_rectangle_from_vertices(
    Filtration_value const* input, Index n_rows, Index n_cols, Out0&& out0, Out1&& out1) {
#ifdef GUDHI_DETAILED_TIMES
  Gudhi::Clock clock;
#endif
  GUDHI_CHECK(n_rows >= 2 && n_cols >= 2,
      std::domain_error("The complex must truly be 2d, i.e. at least 2 rows and 2 columns"));
  Persistence_on_rectangle_V<Filtration_value, Index, output_index> X;
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

// Ideas for improvement:
// * see ideas in the T file.
// * currently, in fill_and_pair, we test edges clockwise, and pair a square with the next edge clockwise. Using opposing directions for those 2 things would have the advantage that we could sometimes simplify the parent tree directly, for instance if we have UR->UL then DR->UR, we could write directly DR->UL (we already do it a bit in the T construction).
// * specialize fill_and_pair_boundary (as in the T construction), probably only relevant for very small (or at least thin) input.
// * Reduce tests in dualize_edge: if we fill the unused boundary with 0 (=exterior), returning one of those boundary squares may be safe, although it may not be faster if it just delays the work to the next find.

#endif  // GUDHI_PERSISTENCE_ON_RECTANGLE_V_H
