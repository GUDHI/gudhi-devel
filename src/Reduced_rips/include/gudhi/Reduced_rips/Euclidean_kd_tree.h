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
 * @file Euclidean_kd_tree.h
 * @author Thomas Burnett, Musashi Koyama
 * @brief Spatial neighbor search over a Euclidean point cloud: a CGAL kd-tree in low ambient dimension, a
 * flat brute-force scan otherwise.
 */

#ifndef REDUCED_RIPS_EUCLIDEAN_KD_TREE_H_
#define REDUCED_RIPS_EUCLIDEAN_KD_TREE_H_

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <optional>
#include <utility>
#include <vector>

#include <CGAL/Epick_d.h>

#include <gudhi/Kd_tree_search.h>
#include <gudhi/Reduced_rips/Helpers.h>

namespace Gudhi {

namespace reduced_rips {

// Spatial search over the cloud. In dimension <= 3 it uses GUDHI's Kd_tree_search. In dimension >= 4,
// it falls back to a flat brute-force scan over the Cloud. Every distance comparison uses l2_dist_2.
// `Index` is the module's stored-index type (std::uint32_t by default): the returned point indices are
// narrowed to it, since every index is < n.
template <class Index = std::uint32_t>
class Euclidean_kd_tree {
  using Cloud = detail::Cloud;

 public:
  Euclidean_kd_tree(const Cloud& pm, bool brute_force) : pm_(pm) {
    if (brute_force) return;  // leave tree_ disengaged; queries then fall back to a brute-force scan
    kd_points_.reserve(pm.n);
    for (std::size_t i = 0; i < pm.n; ++i) kd_points_.emplace_back(pm[i], pm[i] + pm.dim);
    // Kd_tree_search keeps a reference to kd_points_, which outlives the tree (both are members below).
    tree_.emplace(kd_points_);
  }

  Euclidean_kd_tree(const Euclidean_kd_tree&) = delete;
  Euclidean_kd_tree& operator=(const Euclidean_kd_tree&) = delete;
  Euclidean_kd_tree(Euclidean_kd_tree&&) = delete;
  Euclidean_kd_tree& operator=(Euclidean_kd_tree&&) = delete;
  ~Euclidean_kd_tree() = default;

  // True when queries fall back to a flat scan (no kd-tree). The tree can't restrict a spatial query by point
  // index, so index-bounded searches take a different path in brute-force mode.
  [[nodiscard]] bool brute_force() const { return !tree_; }

  // The nearest points to query within an approximate `budget`, ascending by squared distance with ties broken
  // by index. This is the same (distance, index) order the exhaustive per-point scans use, and the returned list
  // is an initial segment of it: the engine's heap frontier resumes positionally inside a refreshed full list, so
  // this list must be an exact prefix of it. Exactly `budget` in the common tie-free case, fewer when a tie
  // straddles the boundary (the whole boundary tie group is dropped). The engine backstops a short list with the
  // full neighbor list. A query that is itself a Cloud point is returned as the nearest (distance 0). Callers
  // filter that out. The ordering keys are computed in T, the scalar the heap orders by, so this prefix agrees
  // with the exhaustive T-keyed scans. The CGAL candidate search itself runs in double and only selects which points
  // reach the tie-trim below.
  template <class T>
  std::vector<Index> nearest_neighbors(const double* query, std::size_t budget) const {
    if (!tree_) {
      std::vector<T> dist(pm_.n);
      for (std::size_t i = 0; i < pm_.n; ++i) dist[i] = detail::l2_dist_2<T>(pm_[i], query, pm_.dim);
      return detail::nearest_within_budget<Index>(0, pm_.n, budget, [&dist](Index x) { return dist[x]; });
    }
    // CGAL's k-nearest search orders equal distances arbitrarily, so its raw top-k is not the canonical
    // (distance, index) prefix. Ask for one extra point and drop the whole tie group at the farthest distance:
    // every point strictly nearer than the (budget+1)-th is unambiguously among the nearest whatever order CGAL
    // used, so what remains is a deterministic prefix of the l2_dist_2 ordering (exactly budget unless ties shorten
    // it). This mirrors the nearest_within_budget rule used by the brute-force and matrix paths.
    Kd_point center(query, query + pm_.dim);
    std::vector<std::pair<Index, T>> near;
    T d_max = T(0);
    for (auto nb : tree_->k_nearest_neighbors(center, static_cast<unsigned int>(budget + 1), true)) {
      auto d = detail::l2_dist_2<T>(pm_[static_cast<std::size_t>(nb.first)], query, pm_.dim);
      near.emplace_back(static_cast<Index>(nb.first), d);
      d_max = std::max(d_max, d);
    }
    near.erase(
        std::remove_if(near.begin(), near.end(), [d_max](const std::pair<Index, T>& pr) { return pr.second >= d_max; }),
        near.end());
    std::sort(near.begin(), near.end(), [](const std::pair<Index, T>& x, const std::pair<Index, T>& y) {
      return x.second != y.second ? x.second < y.second : x.first < y.first;
    });
    std::vector<Index> result;
    result.reserve(near.size());
    for (const auto& pr : near) result.push_back(pr.first);
    return result;
  }

  // All points within the given squared radius of query, as (index, squared distance) pairs. The radius and the
  // returned squared distances are in the scalar type T; the CGAL tree searches in double, and each survivor's
  // distance is recomputed in T afterwards.
  template <class T>
  std::vector<std::pair<Index, T>> points_in_squared_ball(const double* query, T squared_radius) const {
    std::vector<std::pair<Index, T>> result;
    if (!tree_) {
      for (std::size_t i = 0; i < pm_.n; ++i) {
        T d = detail::l2_dist_2<T>(pm_[i], query, pm_.dim);
        if (d <= squared_radius) result.emplace_back(static_cast<Index>(i), d);
      }
      return result;
    }
    Kd_point center(query, query + pm_.dim);
    thread_local std::vector<std::size_t> found;  // per-worker scratch, reused across lune queries
    found.clear();
    tree_->all_near_neighbors2(center, squared_radius, squared_radius, std::back_inserter(found));
    result.reserve(found.size());
    for (std::size_t idx : found)
      result.emplace_back(static_cast<Index>(idx), detail::l2_dist_2<T>(pm_[idx], query, pm_.dim));
    return result;
  }

 private:
  using Kd_kernel = CGAL::Epick_d<CGAL::Dynamic_dimension_tag>;
  using Kd_point = Kd_kernel::Point_d;
  using Kd_tree = Gudhi::spatial_searching::Kd_tree_search<Kd_kernel, std::vector<Kd_point>>;

  Cloud pm_;  // non-owning view (pointer + sizes); pointee outlives this
  std::vector<Kd_point> kd_points_;
  std::optional<Kd_tree> tree_;  // disengaged in brute-force mode; engaged holds the kd-tree over kd_points_
};

}  // namespace reduced_rips

}  // namespace Gudhi

#endif  // REDUCED_RIPS_EUCLIDEAN_KD_TREE_H_
