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
//  See Persistence_on_rectangle.h's own header comment for the general
//  T-construction algorithm this one mirrors (square-centered pairing,
//  dualized here to be vertex-centered); comments below note where and why
//  this file's approach differs.

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

/**
 * @private
 * V-construction analog of Persistence_on_rectangle: persistence of a
 * function given at the vertices of a 2d grid, extended to edges and
 * squares by the MAXIMUM of incident vertices (as opposed to the
 * T-construction's minimum-of-top-cells).
 *
 * All n_rows*n_cols input values are genuine vertices; the derived
 * squares are indexed on an (n_rows-1) x (n_cols-1) grid. Unlike the
 * T-construction, no layer of cells can be dropped: every input value is
 * a genuine vertex of the complex and must be represented. The derived
 * squares number (n_rows-1)*(n_cols-1), exactly the mirror image of the
 * T-construction's square/vertex counts.
 *
 * fill_and_pair() computes a maximal local pairing around each vertex,
 * mirroring the T-construction's square-centered recipe with vertex and
 * square swapped (duality). For vertex i: an incident EDGE "matches" i
 * when i is the argmax of its two endpoints (the edge's value equals
 * input(i)); an incident quadrant SQUARE "matches" i when i is the
 * argmax of all four of its corners. Because `beats` is a strict total
 * order, a given edge only ever matches from the ONE endpoint that beats
 * the other, so there is no ambiguity to resolve and no need to keep two
 * independently-pruned edge lists; each matching edge is handled exactly
 * once, by its one relevant vertex.
 *
 * Walking the four squares around i in cyclic order [DR, DL, UL, UR], a
 * square that matches AND whose two bordering edges both match gets
 * paired with the edge that follows it in the cycle: concretely, that
 * square is merged with its neighbor across that edge (or the exterior,
 * at the grid boundary) in ds_parent_square_ (the H1/dual union-find).
 * If all four qualify, pairing all four would close a cycle, so the
 * last one (UR's, "R") is left alone instead: the square is marked
 * critical directly (self-rooted in ds_parent_square_, with its birth
 * recorded), and R -- the one direction with nowhere else to go -- still
 * pairs the vertex into its R-neighbor in ds_parent_vertex_, exactly as
 * it would have if the last square hadn't qualified at all.
 * Otherwise, among i's matching edges not already consumed by a square-pairing,
 * one is used to merge i into that neighbor in ds_parent_vertex_ (H0);
 * any further such edges are added directly to the shared `edges` list,
 * the same way Persistence_on_rectangle.h's own fill_and_pair adds a
 * leftover edge as critical.
 *
 * That single shared `edges` list plays the same role as in the
 * T-construction: primal() is handed the whole list and physically
 * removes (via remove_if) whatever it actually uses for a real H0
 * merge; dual() then only ever looks at the true remainder, relying on
 * the standard planar-graph fact that a spanning tree's edges are
 * exactly the complement of a spanning tree of the dual graph.
 *
 * Squares are indexed by their own down-right corner, sharing the same
 * index range as vertices (see the class-level comment by `exterior`),
 * which safely allows index 0 to double as the point at infinity for
 * the dual (H1) union-find, exactly as the T-construction does.
 */
template <class Filtration_value, class Index = std::size_t, bool output_index = false>
struct Persistence_on_rectangle_V {
  // As in the T-construction: when output_index is requested, a
  // derived (edge/square) filtration value is paired with the index of
  // the vertex that achieved it, so a caller doing autodiff can trace a
  // birth/death back to the specific input entry it came from, instead
  // of only getting the (already-differentiated-away) numeric value.
  // Vertex births need no such wrapper at all: unlike the T-construction
  // (where a vertex's value is itself derived, as the min of incident
  // squares, so which square achieved it must be tracked), here a
  // vertex's birth is simply input(v) directly, so v's own index already
  // *is* the vertex's "index" for free -- see primal()/dual() below,
  // which branch on output_index directly for those, with no T needed.
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
  // stride (dy) and index range as the vertices themselves, exactly as
  // the T-construction indexes its derived vertices relative to squares:
  // the four squares around vertex i are then simply i, i+1, i+dy,
  // i+dy+1, with no separate square stride or multiplication needed.
  // This leaves the entire top row and left column of this index range
  // unused by any real square (there is no square whose down-right
  // corner is on the grid's own top row or left column) -- in
  // particular, index 0 is never a real square, so it can safely double
  // as the point at infinity for the dual (H1) union-find, exactly as
  // the T-construction does.
  Index exterior = 0;

  // Birth value of each derived square = max of its 4 corners.
  // This is only needed for critical squares, which from fill_and_pair
  // always have the value of the same corner.
  T square_birth(Index s) const { return T(input(s - 1), s - 1); }

  // Union-find forests. ds_parent_vertex_ needs no special "point at
  // infinity": H0 never needs one, only H1 does (Alexander duality).
  std::unique_ptr<Index[]> ds_parent_vertex_;
  std::unique_ptr<Index[]> ds_parent_square_;
  Index& ds_parent_vertex(Index n) { return ds_parent_vertex_[n]; }
  Index& ds_parent_square(Index n) { return ds_parent_square_[n]; }

  // Birth of the single infinite H0 interval: the vertex's own index
  // when output_index, its raw value otherwise -- set in primal().
  std::conditional_t<output_index, Index, Filtration_value> global_min;

  template<class Parent>
  Index ds_find_set_(Index v, Parent&& ds_parent) {
    // Path halving, as in Persistence_on_rectangle's own ds_find_set_.
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
  // Used for a genuine H0 merge (never for marking a vertex critical,
  // which instead sets ds_parent_vertex(i) = i directly): despite taking
  // two vertices, this is not a symmetric union -- child and parent play
  // fixed, distinct roles (child's forest entry is what gets written).
  // Named set_parent_vertex, not e.g. "merge_vertices", for the same
  // reason set_parent_square below isn't "merge_squares".
  void set_parent_vertex(Index child, Index parent) {
    GUDHI_CHECK(child != parent, std::logic_error("Bug: use a direct self-assignment to mark a vertex critical"));
    ds_parent_vertex(child) = parent;
  }
  // Every real square is the unique argmax of its own 4 corners for
  // exactly one of them (a strict total order has a unique max), so
  // this is called with `child` equal to that square from exactly that
  // one vertex's own processing, and never again -- ds_parent_square_
  // is write-only here, same as the T-construction's own fill_and_pair,
  // with no find() needed (path-compression during later reads in
  // primal()/dual() is what actually resolves the chains this builds).
  void set_parent_square(Index child, Index parent) {
    ds_parent_square(child) = parent;
  }

  struct Edge {
    T f;  // derived (max-of-two-endpoints) value, paired with the achieving vertex's index
    Index v1, v2;  // primal (real vertex) endpoints, v1 < v2
    Edge() = default;
    Edge(T f, Index v1, Index v2) : f(f), v1(v1), v2(v2) {}
    bool operator<(Edge const& other) const { return f.first < other.f.first; }
  };
  // A single shared list, mirroring the T-construction's own
  // fill_and_pair: since `beats` is a strict total order, a given edge
  // only ever "matches" (in the fill_and_pair sense) from the one
  // endpoint that beats the other, so there is no genuine ambiguity to
  // resolve and no need for two separately-pruned edge lists. Each
  // matching edge, from its one relevant endpoint's local processing,
  // ends up in exactly one of three states: consumed by a
  // square-pairing (dual), chosen as that vertex's one primal leftover,
  // or -- if matched but neither of those -- added here. primal() then
  // physically removes whatever it actually uses (exactly as the
  // T-construction's own primal() does), leaving only the true
  // remainder for dual() to look at.
  std::vector<Edge> edges;

  // The two squares flanking a real grid edge (v1,v2) are fully
  // determined by the edge's endpoints and direction, so they need not
  // be stored in Edge at all -- just recomputed here, once per edge,
  // only when dual() actually needs them. With the down-right-corner square
  // indexing, the in-bounds case is a plain offset from v1, exactly as for
  // sUL/sUR/sDL/sDR in fill_and_pair; only the boundary is tricky, but that's
  // handled by filling a ring of 0 around the squares in init.
  std::pair<Index, Index> dualize_edge(Index v1, Index v2) const {
    return { v2, v1 + dy + 1 };
  }

  void init(Filtration_value const* input_, Index n_rows, Index n_cols) {
    input_p = input_;
    size_x = n_cols; size_y = n_rows;
    dy = n_cols;
    input_size = n_rows * n_cols;

    // Both arrays are sized like the vertex grid itself (input_size),
    // not the smaller true square count: down-right-corner indexing
    // leaves the top row and left column unused, in exchange for the
    // simpler offset arithmetic used throughout (see the class-level
    // comment by `exterior`).
    ds_parent_vertex_.reset(new Index[input_size]);
    // Every real square's ds_parent_square_ entry gets written exactly
    // once, by its own unique argmax corner, during fill_and_pair (see
    // set_parent_square) -- so, like the T-construction's own
    // fill_and_pair, there is no need to eagerly fill it here.
    // `exterior` itself is never a "child" in any merge, so it
    // needs pre-initializing. In addition, we surround the real squares
    // with 0 so dualize_edge does not need to check boundary cases.
    ds_parent_square_.reset(new Index[input_size+dy]);
    // ds_parent_square_[exterior] = exterior;
    for (Index i = 0; i < dy; ++i) ds_parent_square_[i] = exterior;
    for (Index i = 0; i < dy; ++i) ds_parent_square_[input_size + i] = exterior;
    for (Index i = 1; i < size_y; ++i) ds_parent_square_[dy * i] = exterior;

    edges.reserve(input_size / 2);  // same rough order-of-magnitude estimate as T's
  }

  bool beats(Index a, Filtration_value fa, Index b) const {
    // True if vertex a is functionally LARGER than vertex b. Takes fa
    // (input(a)) rather than looking it up itself: every call site is of
    // the form beats(i, f, neighbor) where f = input(i) is already sitting
    // in a local the caller computed once at function entry, so passing it
    // in avoids re-reading input(a) on every single call; input(b) is the
    // one genuinely-fresh read, so it stays looked up in here.
    GUDHI_CHECK(a != b, std::logic_error("Bug: comparing a vertex to itself"));
    Filtration_value fb = input(b);
    if (fa > fb) return true;
    if (fa < fb) return false;
    return a > b;  // tie-break: larger index wins ("beats")
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

  // See the class-level comment above for the general fill_and_pair
  // algorithm this implements, including the all4 cycle-avoidance case.
  //
  // Interior fast path: has_u/has_d/has_l/has_r (and all four has_ul/ur/
  // dl/dr) are unconditionally true away from the boundary, so every
  // existence check and every exterior fallback in the general version
  // below is dead weight here -- skipped entirely.
  //
  // R/D/L/U and the four diagUL/UR/DL/DR quadrant checks are lambdas,
  // not eagerly-computed bools, matching the style Persistence_on_
  // rectangle.h itself uses (its left()/right()/up_left() etc.): the
  // decision tree below tests each one at most once per cell, in a fixed
  // order (R, D, DR, L, DL, U, UL, UR), so a lambda costs nothing beyond
  // an eager bool on the happy path, but never computes one that a given
  // cell's branch never needs. Because the tree always confirms U and L
  // (resp. U/R, D/L, D/R) true before ever testing diagUL (resp. diagUR,
  // diagDL, diagDR) -- checked exhaustively, for all 256 combinations of
  // the 8 booleans this tree branches on -- the diag lambdas need no
  // redundant "U() && L() &&" prefix of their own.
  //
  // Within any one leaf, the emitted actions (a set_parent_square call,
  // the single ds_parent_vertex(i) assignment, and any
  // edges.emplace_back calls) touch disjoint state -- each of the four
  // quadrant squares gets set_parent_square'd at most once,
  // ds_parent_vertex(i) is set exactly once, and edge insertion order is
  // unobserved (edges are only ever later compared by value, never by
  // position) -- so they have no ordering dependency on one another, and
  // neither does testing a diagXX condition (a pure read). That's what
  // makes it sound to physically hoist whatever's common between two
  // sibling branches out of their if/else, wherever it happens to sit
  // in each
  // branch: it can only remove a duplicate, never change what runs, on
  // any given execution path. Where two branches differ only in one
  // direction's leftover fate (e.g. whether U ends up consumed by a
  // square-pairing or left as its own critical edge), this collapses
  // what would otherwise be two full duplicated copies of the *other*
  // three directions' handling down to a shared copy plus a two-line
  // if/else for just that one direction -- e.g. the qUL/qUR interplay
  // below is shared between the qUL-true and qUL-false cases wherever
  // the two aren't entangled by the (rare) all4 interaction.
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

  // General version, for boundary vertices where some neighbors may not
  // exist. all4 (see fill_and_pair_interior) is impossible here -- a
  // boundary vertex is always missing at least one of has_u/has_d/has_l/
  // has_r -- so, unlike the interior version, there is no all4 branch at
  // all here (see the GUDHI_CHECK inside).
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
    // Unlike fill_and_pair_interior, all4 (all four quadrants qualifying
    // at once) cannot happen here: it needs all four of has_u/has_d/
    // has_l/has_r true, which is false by construction for any boundary
    // vertex (at least one side is missing). So qUR's own square-pairing
    // always fires unconditionally when qUR holds -- there's no cycle to
    // avoid, and so no all4 special case to write here at all. Kept as a
    // GUDHI_CHECK (free: it only reuses already-computed booleans) so the
    // invariant stays documented and would fail loudly if this function
    // were ever called on a non-boundary vertex by mistake.
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

  // H0: standard elder-rule Kruskal directly on the real vertices/edges.
  // Vertex birth is simply input(v) -- no derived value to reconstruct,
  // unlike the T-construction's data_vertex.
  //
  // As in the T-construction's own primal(), edges that ARE used here (a
  // genuine merge) are physically removed from `edges` via remove_if
  // before returning; dual() then only has to look at whatever primal
  // did NOT need -- valid by the standard planar-graph fact that a
  // spanning tree's edges are exactly the complement of a spanning tree
  // of the dual graph.
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

  // H1: Alexander duality on the square graph + point at infinity
  // (exterior), processing edges from LARGEST to smallest, exactly as
  // the T-construction's own dual() does -- only the source of the
  // square values and edges differs.
  template<class Out>
  void dual(Out&& out) {
    for (auto& e : boost::adaptors::reverse(edges)) {
      auto [dv1, dv2] = dualize_edge(e.v1, e.v2);
      Index a = ds_find_set_square(dv1), b = ds_find_set_square(dv2);
      // fill_and_pair's local matching is maximal so dual()
      // never actually sees an edge whose two sides are already
      // connected -- so, exactly as in the
      // T-construction's own dual(), a==b here indicates real
      // corruption, not an expected case to skip.
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
 *   dimension 1 (all necessarily finite: H1 needs no point at infinity
 *   beyond the one already absorbed into the exterior).
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
