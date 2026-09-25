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
 * @file Euclidean_geometry.h
 * @author Thomas Burnett, Musashi Koyama
 * @brief The Euclidean geometry policy: coordinates plus a kd-tree, working in squared distances, with the
 * paper's lens-ball / wide-angle accelerations.
 */

#ifndef REDUCED_RIPS_EUCLIDEAN_GEOMETRY_H_
#define REDUCED_RIPS_EUCLIDEAN_GEOMETRY_H_

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

#include <gudhi/Reduced_rips/Delaunay_edges.h>
#include <gudhi/Reduced_rips/Euclidean_kd_tree.h>
#include <gudhi/Reduced_rips/Helpers.h>
#include <gudhi/Reduced_rips/Lune_builder.h>
#include <gudhi/Reduced_rips/Relative_neighborhood_graph.h>

namespace Gudhi {

namespace reduced_rips {

// Geometry policy (a model of the geometry concept; see concept/Geometry.h) over a Euclidean point cloud. It
// works in squared distances: the algorithm relies only on the *ordering* of dist, so dist returns d^2 and
// to_distance is sqrt. Coordinates stay double, but every candidate is re-checked with l2_dist_2<T>, so the barcode
// keeps T's precision.
template <class T, class Index_ = std::uint32_t>
class Euclidean_geometry {
  using Cloud = detail::Cloud;

 public:
  using Filtration_value = T;
  using Index = Index_;  // stored-index type (point indices and edge ids), std::uint32_t by default
  Euclidean_geometry(const Cloud& pm, const Euclidean_kd_tree<Index>& kd) : pm_(pm), kd_(&kd) {}
  [[nodiscard]] std::size_t size() const { return pm_.n; }
  // Ordering scale: the squared Euclidean distance.
  [[nodiscard]] T dist(std::size_t i, std::size_t j) const { return detail::l2_dist_2<T>(pm_[i], pm_[j], pm_.dim); }
  // Maps the squared ordering scale back to a true distance for the output barcode.
  [[nodiscard]] static T to_distance(T squared) { return std::sqrt(squared); }
  [[nodiscard]] std::vector<Index> nearest(std::size_t i, std::size_t budget) const {
    return kd_->template nearest_neighbors<T>(pm_[i], budget);
  }
  // The nearest points to i restricted to index > i, ascending by distance (ties by index; a boundary tie group
  // is dropped rather than split, so the list is always a prefix of the neighbors_above ordering). The count
  // shrinks as i rises. Brute-force mode searches the above-i tail directly under a scaled budget. The kd-tree
  // can't bound a spatial query by index, so it searches globally and trims.
  [[nodiscard]] std::vector<Index> nearest_neighbors_above(std::size_t i, std::size_t budget) const {
    if (kd_->brute_force()) return brute_nearest_above(i, detail::above_budget(budget, i, pm_.n));
    std::vector<Index> result = nearest(i, budget + 1);
    detail::keep_above(i, result);
    return result;
  }
  // The k nearest points to i with index > i, ascending by squared distance (ties by index): a prefix of the
  // full above-i ordering, computed by a direct brute scan of the above-i tail (not the kd-tree). The engine
  // grows this on demand with a doubling k when the heap frontier outruns the k-nearest prefetch; a brute
  // partial-sort of the tail is far cheaper here than repeated large kd k-nearest queries, and storing only k
  // (rather than the whole tail) keeps peak memory bounded. When k >= the above-i count it is the full list.
  [[nodiscard]] std::vector<Index> neighbors_above(std::size_t i, std::size_t k) const {
    const std::size_t off = i + 1, count = pm_.n - off;
    std::vector<T> dist(count);
    for (std::size_t j = off; j < pm_.n; ++j) dist[j - off] = detail::l2_dist_2<T>(pm_[i], pm_[j], pm_.dim);
    return detail::smallest_indices_by<Index>(off, count, k, [&dist, off](std::size_t j) { return dist[j - off]; });
  }
  // Early-stop target: the RNG cycle rank (the number of finite H1 bars), from the per-dimension routine.
  [[nodiscard]] std::size_t rng_early_stop_target() const {
    if (pm_.dim == 2) return rng_cycle_rank_delaunay<T>(pm_, *kd_, delaunay_edges_2d);
    if (pm_.dim == 3) return rng_cycle_rank_delaunay<T>(pm_, *kd_, delaunay_edges_3d);
    return rng_cycle_rank_general<T>(pm_, *kd_);
  }
  [[nodiscard]] Lune_result<T, Index> lune(const Batch_edge<T, Index>& e,
                                           const detail::Edge_map<std::size_t, Index>& one_simp_to_idx,
                                           std::size_t n) const {
    return Lune_builder<T, Index>(e, one_simp_to_idx, n).build_euclidean(pm_, *kd_);
  }

 private:
  // Brute-force: the `budget` nearest points above i, distances computed only for the above-i tail. Used when the
  // kd-tree is disengaged. The budget is already scaled by the caller (detail::above_budget).
  [[nodiscard]] std::vector<Index> brute_nearest_above(std::size_t i, std::size_t budget) const {
    const std::size_t off = i + 1, count = pm_.n - off;
    std::vector<T> dist(count);
    for (std::size_t j = off; j < pm_.n; ++j) dist[j - off] = detail::l2_dist_2<T>(pm_[i], pm_[j], pm_.dim);
    return detail::nearest_within_budget<Index>(off, count, budget, [&dist, off](std::size_t j) { return dist[j - off]; });
  }

  Cloud pm_;                            // non-owning view; pointee outlives this
  const Euclidean_kd_tree<Index>* kd_;  // non-owning, never null; kd-tree is move-only so stored by pointer
};

}  // namespace reduced_rips

}  // namespace Gudhi

#endif  // REDUCED_RIPS_EUCLIDEAN_GEOMETRY_H_
