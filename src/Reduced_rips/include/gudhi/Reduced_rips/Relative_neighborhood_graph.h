/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Thomas Burnett, Musashi Koyama
 *
 *    Copyright (C) 2026 Thomas Burnett, Musashi Koyama
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file Relative_neighborhood_graph.h
 * @author Thomas Burnett, Musashi Koyama
 * @brief The relative-neighborhood-graph cycle rank, used as the reduction's early-stop target (the number of
 * finite degree-1 bars), via a Delaunay triangulation in 2D/3D and a direct construction otherwise.
 */

#ifndef REDUCED_RIPS_RELATIVE_NEIGHBORHOOD_GRAPH_H_
#define REDUCED_RIPS_RELATIVE_NEIGHBORHOOD_GRAPH_H_

#include <algorithm>
#include <cstddef>
#include <type_traits>
#include <utility>
#include <vector>

#include <gudhi/Reduced_rips/Helpers.h>

namespace Gudhi {

namespace reduced_rips {

// ---- Relative neighborhood graph (RNG) cycle rank -----------------------------------------------------
// total_death = the RNG cycle rank |E| - |V| + 1 (the RNG is connected): the number of finite H1 bars, which
// drives the reduction's early-stop. Each function below returns this rank for one geometry. |V| is n, except
// the 2D/3D Delaunay can merge coincident input points into one vertex (see rng_cycle_rank_delaunay).
//
// |E| is the *true* RNG = the open lune: edge (a,b) is in the RNG iff no point lies *strictly* inside its lune
// (strict `<` on both endpoint distances); a point exactly on the boundary does NOT remove the edge. This is
// deliberately different from the lune-occupancy test inside the reduction (in_lune, in Lune_builder), whose
// lexical tie-break assigns boundary points to 2-simplices for the homology algorithm and must NOT decide RNG
// membership: using in_lune here would over-eliminate boundary edges and undercount. Do not "unify" them.
//
// The RNG is computed with the fastest strategy for the ambient dimension (Delaunay-based in 2D/3D, direct
// otherwise).

// RNG cycle rank from the Delaunay edges (an Urquhart superset of the RNG, see Delaunay_edges.h): discard every edge whose open
// lune contains a point, then return |E| - |V| + 1. |V| is the number of *participating* vertices, counted
// from the surviving edges rather than taken as n, because CGAL's 2D/3D Delaunay merges coincident input
// points into one vertex; the n - |V| merged duplicates contribute only zero-persistence cycles, which the
// reduction drops, so counting them out keeps the rank exact. One ball query around an endpoint suffices,
// since each candidate is re-tested against the other endpoint.
template <class T, class KdTree, typename DelaunayEdges>
std::size_t rng_cycle_rank_delaunay(const detail::Cloud& pm, const KdTree& kd_tree, DelaunayEdges delaunay_edges) {
  std::vector<std::pair<std::size_t, std::size_t>> possible_edges = delaunay_edges(pm);
  std::vector<char> seen(pm.n, 0);
  std::size_t kept = 0, vertices = 0;
  for (const auto& edge : possible_edges) {
    std::size_t a = edge.first, b = edge.second;
    T r = detail::l2_dist_2<T>(pm[a], pm[b], pm.dim);
    // widen_radius widens only the ball-query radius (so the strict test below sees every candidate); occupancy
    // is the open lune: a point strictly inside both endpoint balls. Boundary points do not remove the edge.
    auto ball = kd_tree.points_in_squared_ball(pm[a], detail::widen_radius(r));
    bool lune_occupied = std::any_of(ball.begin(), ball.end(), [&](const auto& pr) {
      std::size_t k = pr.first;
      if (k == a || k == b) return false;
      // The ball query is centered at a, so pr.second already is d(k,a)^2 in T; testing it first spares
      // the b-side distance for the candidates in the widened shell.
      if (!(pr.second < r)) return false;
      return detail::l2_dist_2<T>(pm[b], pm[k], pm.dim) < r;
    });
    if (!lune_occupied) {
      ++kept;
      for (std::size_t v : {a, b})
        if (seen[v] == 0) {
          seen[v] = 1;
          ++vertices;
        }
    }
  }
  if (kept == 0) return 0;     // all points coincident: one merged vertex, no edges, rank 0
  return kept - vertices + 1;  // |E| - |V| + 1 for the connected RNG; >= 0
}

// A candidate edge as (lo, hi, length in the geometry's scale), ordered by length then index for the set operations. Shared
// by the direct (general-dimension) and distance-matrix RNG supergraph builds. The length is in the scalar
// type T supplied by the distance callable.
template <class T>
struct Rng_edge {
  std::size_t i, j;
  T length;
  Rng_edge(std::size_t a, std::size_t b, T len) : i(std::min(a, b)), j(std::max(a, b)), length(len) {}
  bool operator<(const Rng_edge& o) const {
    if (length != o.length) return length < o.length;
    if (i != o.i) return i < o.i;
    return j < o.j;
  }
  bool operator==(const Rng_edge& o) const { return i == o.i && j == o.j; }
};

// Phase 1 of the direct RNG construction: an O(n^2) supergraph of the RNG, returned sorted and deduped.
// `dist` is a distance-like callable (i, j) to the geometry's ordering scale; both the coordinate and matrix paths
// share this body and differ only in how they supply distances. e_i is a sorted vector, front() is the
// current minimum, and the survivors of each pruning pass are compacted in place, preserving sorted order
// with no per-iteration allocation. e_all is accumulated then deduped once at the end. An edge (i,j) can
// be emitted from both its endpoint passes, which the old std::set merged.
//
// An edge is eliminated only when a point lies *strictly* inside its open lune, mirroring the strict
// occupancy test of the pruning phase below: eliminating on a boundary tie (equal distances broken by
// index, as the Rng_edge ordering would) drops edges the strict open-lune RNG keeps, undercounting the
// cycle rank and disagreeing with the 2D/3D Delaunay path on tied inputs.
template <class Dist>
auto rng_supergraph(std::size_t n, Dist dist) {
  using T = std::invoke_result_t<Dist, std::size_t, std::size_t>;
  std::vector<Rng_edge<T>> e_all;
  for (std::size_t i = 0; i < n; ++i) {
    std::vector<Rng_edge<T>> e_i;
    e_i.reserve(n - 1);
    for (std::size_t j = 0; j < n; ++j)
      if (i != j) e_i.emplace_back(i, j, dist(i, j));
    std::sort(e_i.begin(), e_i.end());
    while (!e_i.empty()) {
      Rng_edge<T> min_edge = e_i.front();
      e_all.push_back(min_edge);
      // Every entry of e_i has i as an endpoint, so `other` (min_edge's non-i endpoint) is loop-invariant.
      const std::size_t other = min_edge.i == i ? min_edge.j : min_edge.i;
      std::size_t w = 0;
      for (std::size_t t = 1; t < e_i.size(); ++t) {
        const Rng_edge<T>& edge = e_i[t];
        const std::size_t far = edge.i == i ? edge.j : edge.i;
        // `other` lies strictly inside the open lune of (i, far) iff both d(i, other) (== min_edge.length)
        // and d(other, far) fall strictly below d(i, far) (== edge.length, already cached).
        const bool strictly_dominated = min_edge.length < edge.length && dist(other, far) < edge.length;
        if (!strictly_dominated) e_i[w++] = edge;
      }
      e_i.erase(e_i.begin() + static_cast<std::ptrdiff_t>(w), e_i.end());
    }
  }
  std::sort(e_all.begin(), e_all.end());
  e_all.erase(std::unique(e_all.begin(), e_all.end()), e_all.end());
  return e_all;
}

// Direct RNG construction for general dimension. Phase 1 builds an O(n^2) RNG superset, phase 2 prunes
// edges whose lune is non-empty; returns the cycle rank |E| - n + 1. Edge lengths are squared distances
// throughout. No vertex merging here (unlike the Delaunay path), so |V| = n.
template <class T, class KdTree>
std::size_t rng_cycle_rank_general(const detail::Cloud& pm, const KdTree& kd_tree) {
  auto e_all = rng_supergraph(pm.n, [&pm](std::size_t i, std::size_t j) { return detail::l2_dist_2<T>(pm[i], pm[j], pm.dim); });

  // Phase 2: eliminate edges whose open lune contains a point (strictly inside both endpoint balls). widen_radius
  // widens only the ball-query radius; boundary points do not remove the edge.
  std::size_t count = 0;
  for (const auto& edge : e_all) {
    std::size_t a = edge.i, b = edge.j;
    auto ball = kd_tree.points_in_squared_ball(pm[a], detail::widen_radius(edge.length));
    bool lune_occupied = std::any_of(ball.begin(), ball.end(), [&](const auto& pr) {
      std::size_t k = pr.first;
      if (k == a || k == b) return false;
      // pr.second is d(k,a)^2 (the query is centered at a); test it before computing the b side.
      if (!(pr.second < edge.length)) return false;
      return detail::l2_dist_2<T>(pm[b], pm[k], pm.dim) < edge.length;
    });
    if (!lune_occupied) ++count;
  }
  return count - pm.n + 1;
}

// Exact early-stop target for a bare distance matrix: the RNG cycle rank |E| - n + 1, i.e. the number of
// finite H1 bars. Mirrors rng_cycle_rank_general (O(n^2) supergraph build, then O(|E| n) open-lune pruning,
// reading every distance from the geometry). No vertex merging here, so |V| = n.
template <class Geom>
std::size_t rng_cycle_rank_matrix(const Geom& g) {
  std::size_t n = g.size();
  auto e_all = rng_supergraph(n, [&g](std::size_t i, std::size_t j) { return g.dist(i, j); });

  // Occupancy is the open lune: a point strictly inside both endpoint balls.
  std::size_t kept = 0;
  for (const auto& edge : e_all) {
    std::size_t a = edge.i, b = edge.j;
    bool occupied = false;
    for (std::size_t k = 0; k < n && !occupied; ++k) {
      if (k == a || k == b) continue;
      if (g.dist(a, k) < edge.length && g.dist(b, k) < edge.length) occupied = true;
    }
    if (!occupied) ++kept;
  }
  return kept - n + 1;
}

}  // namespace reduced_rips

}  // namespace Gudhi

#endif  // REDUCED_RIPS_RELATIVE_NEIGHBORHOOD_GRAPH_H_
