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
 * @file Lune_builder.h
 * @author Thomas Burnett, Musashi Koyama
 * @brief Per-edge lune computation: turns a candidate edge into the 2-simplices it contributes (one apparent
 * simplex, or one boundary column per lune connected component).
 */

#ifndef REDUCED_RIPS_LUNE_BUILDER_H_
#define REDUCED_RIPS_LUNE_BUILDER_H_

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <random>
#include <utility>
#include <vector>

#include <boost/pending/disjoint_sets.hpp>

#include <gudhi/Reduced_rips/Helpers.h>

namespace Gudhi {

namespace reduced_rips {

namespace detail {

// 1-simplices are keyed by a single integer lo*n + hi (lo < hi). The key stays 64-bit: with up to 2^32 points
// the product lo*n overflows a 32-bit index, so the endpoints widen to std::size_t here.
inline std::size_t pack_edge(std::size_t lo, std::size_t hi, std::size_t n) { return (lo * n) + hi; }

}  // namespace detail

// ---- Per-edge "lune" computation --------------------------------------------------------------------------
// For a candidate edge (a,b) of length r, this determines the 2-simplices it contributes: either one
// apparent 2-simplex (whose boundary column is returned ready to file under its own pivot), or one boundary
// column per connected component of the lune (returned to be reduced).
// The same struct serves as the engine's heap entry and as the lune input. While an edge sits in the heap,
// `id` transiently holds the neighbor-list frontier position `t`. When the edge is popped into a batch that
// slot is overwritten with the assigned 1-simplex id (see Persistence_engine::pop_batch), so by the time the
// lune sees it `id` is the 1-simplex id. `Index` is the module's stored-index type (point indices and edge
// ids alike; std::uint32_t by default), narrow to keep the heap and the neighbor lists compact.
template <class T, class Index = std::uint32_t>
struct Batch_edge {
  Index a, b;  // endpoints, a < b
  T r;         // length of (a,b), in the geometry's scale
  Index id;    // 1-simplex id (heap frontier position t while queued); also the eye-sampling RNG seed
};

// The 2-simplices a candidate edge contributes, as parallel boundary columns and their death values (in the geometry's
// scale). The number of columns tells the engine how to file them, so no separate flags are needed:
//   - 0 columns: the lune contributes no 2-simplex (empty);
//   - 1 column:  a single apparent 2-simplex, filed directly under its own pivot (no reduction, no bar);
//   - >1 column: one boundary column per lune component, each reduced.
template <class T, class Index = std::uint32_t>
struct Lune_result {
  std::vector<std::vector<Index>> cols;  // boundary column(s), each an ascending list of edge ids
  std::vector<T> deaths;                 // 2-simplex diameter (geometry scale) for each column in `cols`

  Lune_result() = default;  // empty: no 2-simplex
  // A single 2-simplex from its ascending 3-edge boundary column and diameter (in the geometry's scale).
  Lune_result(std::vector<Index> column, T death) {
    cols.push_back(std::move(column));
    deaths.push_back(death);
  }

  // The candidate edge contributed no 2-simplex.
  [[nodiscard]] bool empty() const { return cols.empty(); }
  // Exactly one boundary column: a single apparent 2-simplex, filed under its own pivot with no recorded bar.
  [[nodiscard]] bool is_apparent() const { return cols.size() == 1; }
};

// Turns the lune points of one candidate edge into a Lune_result. The geometry-specific front-ends
// (build_euclidean / build_matrix) gather the lune points and, for the Euclidean geometry, may apply the
// lens/eye certificates of the reference paper. The work that does not depend on coordinates (the connected
// components of the lune points thresholded at r, and the boundary columns they yield) is shared by both.
template <class T, class Index = std::uint32_t>
class Lune_builder {
  using Cloud = detail::Cloud;
  template <class K, class V>
  using Edge_map = detail::Edge_map<K, V>;

 public:
  Lune_builder(const Batch_edge<T, Index>& e, const Edge_map<std::size_t, Index>& one_simp_to_idx, std::size_t n)
      : a_(e.a), b_(e.b), r_(e.r), id_(e.id), one_simp_to_idx_(&one_simp_to_idx), n_(n) {}

  // Euclidean front-end: gather the lune points from a midpoint ball query, apply the lens-ball fast path and
  // the wide-angle ("eye") single-component certificate, then defer the component analysis to the shared code.
  // The kd-tree is a template parameter (in practice Euclidean_kd_tree) so this header, which the CGAL-free
  // distance-matrix path also uses, carries no CGAL includes.
  template <class KdTree>
  [[nodiscard]] Lune_result<T, Index> build_euclidean(const Cloud& pm, const KdTree& kd_tree) const {
    const std::size_t a = a_, b = b_, dim = pm.dim;
    const T r = r_;

    // A zero-length edge (coincident endpoints) bounds no 2-simplex of positive diameter, so its lune is empty.
    if (r == T(0)) return Lune_result<T, Index>{};

    // Midpoint of (a,b): center of the candidate ball query.
    thread_local std::vector<double> mid_point;
    mid_point.resize(dim);
    for (std::size_t i = 0; i < dim; ++i) mid_point[i] = (pm[a][i] + pm[b][i]) / 2.0;
    // 0.75 * r is the (squared) radius necessary for a ball centered at the midpoint to cover the lune
    std::vector<std::pair<Index, T>> ball_mid =
        kd_tree.points_in_squared_ball(mid_point.data(), detail::widen_radius(T(0.75) * r));

    // Fast path: a candidate inside the (inscribed) lens ball certifies the lune is a single component.
    auto lens = std::find_if(ball_mid.begin(), ball_mid.end(), [&](const std::pair<Index, T>& pr) {
      return pr.second <= detail::lens_ball_factor<T> * r && both_edges_id(pr.first);
    });
    if (lens != ball_mid.end()) return single(lens->first);

    // Slower path: prune the midpoint-ball candidates to the lune, then sort the survivors.
    thread_local std::vector<Index> r_ab;  // per-worker scratch, reused across lune queries
    r_ab.clear();
    r_ab.reserve(ball_mid.size());
    for (const auto& pr : ball_mid) {
      Index k = pr.first;
      if (k == a) continue;
      T dist_ka_sq = detail::l2_dist_2<T>(pm[a], pm[k], dim);
      T dist_kb_sq = detail::l2_dist_2<T>(pm[b], pm[k], dim);
      if (in_lune(dist_ka_sq, dist_kb_sq, r, a, b, k) && both_edges_id(k)) r_ab.push_back(k);
    }
    std::sort(r_ab.begin(), r_ab.end());
    if (r_ab.empty()) return Lune_result<T, Index>{};

    // Heuristic (only worth sampling with more than two lune points): a wide-angle ("eye") point certifies a
    // single component without union-find. The seed is the edge id (id_), so the sampling is deterministic and
    // each parallel worker carries its own; a wrong guess only forgoes the shortcut (the exact union-find still
    // runs below) while emitting the same simplex, so the barcode does not depend on the random number generator.
    bool single_component_hint = false;
    if (r_ab.size() > 2) {
      std::mt19937 rng(id_);
      auto n_check = static_cast<std::size_t>(std::sqrt(double(r_ab.size())));
      for (std::size_t j = 0; j < n_check && !single_component_hint; ++j) {
        std::size_t temp_idx = r_ab[rng() % r_ab.size()];
        const double *pa = pm[a], *pb = pm[b], *pt = pm[temp_idx];
        double dot = 0.0, uu = 0.0, vv = 0.0;
        for (std::size_t k = 0; k < dim; ++k) {
          double u = pa[k] - pt[k], v = pb[k] - pt[k];
          dot += u * v;
          uu += u * u;
          vv += v * v;
        }
        // A single wide-angle sample (angle > 5*pi/6) certifies one component; stop sampling once one is found.
        // Note: cos^2(5*pi/6) = 3/4
        single_component_hint = dot < 0.0 && 4.0 * dot * dot > 3.0 * uu * vv;
      }
    }

    return from_lune_points(
        r_ab, [&pm, dim](std::size_t i, std::size_t j) { return detail::l2_dist_2<T>(pm[i], pm[j], dim); },
        single_component_hint);
  }

  // Matrix front-end: with no coordinates there is no midpoint, hence none of the Euclidean shortcuts.
  // Candidates are gathered from the row of endpoint a (the lune is contained in the closed ball of radius r
  // about a), filtered exactly by in_lune, and the components are always found by the exact union-find.
  template <class Geom>
  [[nodiscard]] Lune_result<T, Index> build_matrix(const Geom& g) const {
    const std::size_t a = a_, b = b_, np = g.size();
    const T r = r_;

    // A zero-length edge (coincident endpoints) bounds no 2-simplex of positive diameter, so its lune is empty.
    if (r == T(0)) return Lune_result<T, Index>{};

    // Lune points of (a,b): scan a's row (lune is a subset of the closed ball of radius r about a), keeping
    // those that pass the exact in_lune test against b. k increases, so r_ab is already ascending in index.
    std::vector<Index> r_ab;
    for (std::size_t k = 0; k < np; ++k) {
      if (k == a || k == b) continue;
      T dist_ka = g.dist(a, k);
      if (dist_ka > r) continue;
      T dist_kb = g.dist(b, k);
      if (in_lune(dist_ka, dist_kb, r, a, b, k)) r_ab.push_back(static_cast<Index>(k));
    }
    if (r_ab.empty()) return Lune_result<T, Index>{};

    return from_lune_points(
        r_ab, [&g](std::size_t i, std::size_t j) { return g.dist(i, j); }, /*single_component_hint=*/false);
  }

 private:
  // Lexicographic order of the sorted index pairs {x,y} < {p,q}. Used for the lune-boundary tie-break.
  [[nodiscard]] static bool sorted_pair_less(std::size_t x, std::size_t y, std::size_t p, std::size_t q) {
    std::size_t xlo = std::min(x, y), xhi = std::max(x, y);
    std::size_t plo = std::min(p, q), phi = std::max(p, q);
    return xlo != plo ? xlo < plo : xhi < phi;
  }

  // True if a point k lies in the lune of edge (a,b) at threshold `thresh`: within `thresh` of both endpoints,
  // where a point sitting exactly on a boundary is admitted only when the index tie-break assigns it to this
  // edge. dist_ka and dist_kb are the distances from k to endpoints a and b (the paper's d(x,y) and d(x,z)) and
  // thresh is the edge length (the paper's r).
  [[nodiscard]] static bool in_lune(T dist_ka, T dist_kb, T thresh, std::size_t a, std::size_t b, std::size_t k) {
    if (dist_ka < thresh && dist_kb < thresh) return true;
    if (dist_ka == thresh && dist_kb < thresh) return sorted_pair_less(a, k, a, b);
    if (dist_ka < thresh && dist_kb == thresh) return sorted_pair_less(b, k, a, b);
    if (dist_ka == thresh && dist_kb == thresh) return sorted_pair_less(a, k, a, b) && sorted_pair_less(b, k, a, b);
    return false;
  }

  // True when both edges (a,c) and (b,c) already have a 1-simplex id, i.e. were popped no later than (a,b).
  // Subtlety: ids for the WHOLE current batch are assigned serially before the lunes run, so this lookup can
  // also see edges popped after (a,b) in the same batch. That never puts a later id into an emitted column,
  // for two reasons that must be preserved together. First, lens fast-path candidates are strictly interior,
  // so both their edges are strictly shorter than r and hence already assigned. Second, the slow path tests
  // in_lune BEFORE this guard, and in_lune admits a boundary tie only when the tied edge lexicographically
  // precedes (a,b), which is exactly the heap's tie order. Reordering the `in_lune(...) && both_edges_id(k)`
  // conjunction, or loosening the lens threshold, would let a same-batch later edge become the column pivot
  // and silently corrupt the apparent-pivot filing in Phase C.
  [[nodiscard]] bool both_edges_id(std::size_t c) const {
    return one_simp_to_idx_->contains(detail::pack_edge(std::min<std::size_t>(a_, c), std::max<std::size_t>(a_, c), n_)) &&
           one_simp_to_idx_->contains(detail::pack_edge(std::min<std::size_t>(b_, c), std::max<std::size_t>(b_, c), n_));
  }

  // A single 2-simplex (a, b, c): its boundary column paired with the diameter r (in the geometry's scale; c lies in
  // the closed lune, so (a,b) is the longest edge and r is the death).
  [[nodiscard]] Lune_result<T, Index> single(std::size_t c) const { return {column_of(c), r_}; }

  // Reduce the lune points r_ab (ascending, non-empty) to the final Lune_result. `dist` is a distance-like
  // callable (i, j) to the geometry's ordering scale; `single_component_hint` lets the Euclidean eye certificate skip
  // the union-find when it has already certified that the lune points form a single component.
  template <class Dist>
  [[nodiscard]] Lune_result<T, Index> from_lune_points(const std::vector<Index>& r_ab, Dist dist,
                                                       bool single_component_hint) const {
    std::size_t n_rab = r_ab.size();
    if (n_rab == 1 || single_component_hint) return single(r_ab[0]);

    // One representative (a global lune-point index) per connected component of the lune points thresholded
    // at r. Two points is the cheap special case; otherwise defer to the all-pairs union-find.
    std::vector<Index> reps;
    if (n_rab == 2) {
      if (dist(r_ab[0], r_ab[1]) < r_) return single(r_ab[0]);  // one component
      reps = {r_ab[0], r_ab[1]};
    } else {
      reps = component_representatives(r_ab, dist);
      if (reps.size() == 1) return single(r_ab[0]);
    }

    // Multiple components: one boundary column per component, each dying at the diameter r (each
    // representative lies in the closed lune of (a,b)).
    Lune_result<T, Index> res;
    res.cols.reserve(reps.size());
    res.deaths.reserve(reps.size());
    for (Index rep : reps) {
      res.cols.push_back(column_of(rep));
      res.deaths.push_back(r_);
    }
    return res;
  }
  // One representative lune-point index per connected component of r_ab of distance less than r. An all-pairs
  // union-find over the local positions 0..|r_ab| merges two positions whose distance is below r. The
  // roots are then deduped to one per component and mapped back to global indices.
  template <class Dist>
  [[nodiscard]] std::vector<Index> component_representatives(const std::vector<Index>& r_ab, Dist dist) const {
    const std::size_t n_rab = r_ab.size();
    thread_local std::vector<std::size_t> rank, parent;  // per-worker union-find scratch, reused across lunes
    rank.assign(n_rab, 0);
    parent.assign(n_rab, 0);
    boost::disjoint_sets<std::size_t*, std::size_t*> ds(rank.data(), parent.data());
    for (std::size_t q = 0; q < n_rab; ++q) ds.make_set(q);
    // Count only real merges so we can stop once all points coalesce into one component (after n_rab-1 merges).
    std::size_t remaining = n_rab;
    for (std::size_t k = 0; k < n_rab && remaining > 1; ++k)
      for (std::size_t l = k + 1; l < n_rab && remaining > 1; ++l)
        if (dist(r_ab[k], r_ab[l]) < r_) {
          std::size_t rk = ds.find_set(k), rl = ds.find_set(l);
          if (rk != rl) {
            ds.link(rk, rl);
            --remaining;
          }
        }
    // One local root per component, deduped, then mapped back to global lune-point indices.
    std::vector<std::size_t> roots(n_rab);
    for (std::size_t z = 0; z < n_rab; ++z) roots[z] = ds.find_set(z);
    std::sort(roots.begin(), roots.end());
    roots.erase(std::unique(roots.begin(), roots.end()), roots.end());
    std::vector<Index> reps(roots.size());
    for (std::size_t z = 0; z < roots.size(); ++z) reps[z] = r_ab[roots[z]];
    return reps;
  }

  // Boundary column of the 2-simplex (a, b, c): the ids of its three edges, sorted ascending so the pivot is
  // the last element. a_ < b_ by construction, so the three edges are (a_,b_) and the two (min,max) pairs with
  // c, formed without sorting the vertices. The three ids are distinct, so a fixed min/max network orders them
  // branch-free (no std::sort), without assuming which edge is longest, so length ties are harmless.
  [[nodiscard]] std::vector<Index> column_of(std::size_t c) const {
    const Index id_ab = one_simp_to_idx_->at(detail::pack_edge(a_, b_, n_));
    const Index id_ac = one_simp_to_idx_->at(detail::pack_edge(std::min<std::size_t>(a_, c), std::max<std::size_t>(a_, c), n_));
    const Index id_bc = one_simp_to_idx_->at(detail::pack_edge(std::min<std::size_t>(b_, c), std::max<std::size_t>(b_, c), n_));
    const Index lo = std::min(id_ab, id_ac), hi = std::max(id_ab, id_ac);
    const Index top = std::max(hi, id_bc), mid_hi = std::min(hi, id_bc);
    return {std::min(lo, mid_hi), std::max(lo, mid_hi), top};
  }

  Index a_, b_;
  T r_;
  Index id_;                                                    // 1-simplex id of (a,b); seeds the eye-sampling RNG
  const Edge_map<std::size_t, Index>* one_simp_to_idx_;         // non-owning, never null; outlives this
  std::size_t n_;
};

}  // namespace reduced_rips

}  // namespace Gudhi

#endif  // REDUCED_RIPS_LUNE_BUILDER_H_
