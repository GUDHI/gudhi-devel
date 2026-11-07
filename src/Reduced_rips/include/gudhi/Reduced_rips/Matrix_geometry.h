/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Thomas Burnett
 *
 *    Copyright (C) 2026 Thomas Burnett, Musashi Koyama
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file Matrix_geometry.h
 * @author Thomas Burnett
 * @brief The distance-matrix geometry policy: a bare symmetric matrix of dissimilarities, carried unchanged.
 */

#ifndef REDUCED_RIPS_MATRIX_GEOMETRY_H_
#define REDUCED_RIPS_MATRIX_GEOMETRY_H_

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include <gudhi/Reduced_rips/Helpers.h>
#include <gudhi/Reduced_rips/Lune_builder.h>
#include <gudhi/Reduced_rips/Relative_neighborhood_graph.h>

namespace Gudhi {

namespace reduced_rips {

// Geometry policy (a model of the geometry concept; see concept/Geometry.h) over a bare symmetric distance
// matrix (n by n, row-major). With no coordinates it keeps the supplied distances as-is (dist returns d,
// to_distance is the identity), gathers lune candidates by scanning a matrix row, and always runs union-find.
template <class T, class Index_ = std::uint32_t>
class Matrix_geometry {
 public:
  using Filtration_value = T;
  using Index = Index_;  // stored-index type (point indices and edge ids), std::uint32_t by default
  Matrix_geometry(std::vector<T> distances, std::size_t n) : dist_(std::move(distances)), n_(n) {}
  [[nodiscard]] std::size_t size() const { return n_; }
  // Ordering scale: the raw distance; comparisons and the barcode use it directly.
  [[nodiscard]] T dist(std::size_t i, std::size_t j) const { return dist_[(i * n_) + j]; }
  // The value already is the distance, so mapping it back to an output distance is the identity.
  [[nodiscard]] static T to_distance(T distance) { return distance; }

  // The nearest points to i restricted to index > i, ascending by distance (ties by index). Searches the above-i
  // tail of the matrix row directly under a budget scaled so the count is roughly what it would have been if we
  // had computed all nearest `budget' points, and then trimmed just to those above `i', (see detail::above_budget and
  // detail::nearest_within_budget).
  [[nodiscard]] std::vector<Index> nearest_neighbors_above(std::size_t i, std::size_t budget) const {
    const T* row = &dist_[i * n_];
    const std::size_t off = i + 1, count = n_ - off;
    return detail::nearest_within_budget<Index>(off, count, detail::above_budget(budget, i, n_),
                                                [row](std::size_t x) { return row[x]; });
  }

  // The k nearest indices > i, ascending by distance from i (ties by index): a prefix of the full above-i
  // ordering, grown on demand by the engine with a doubling k. When k >= the above-i count it is the full list.
  [[nodiscard]] std::vector<Index> neighbors_above(std::size_t i, std::size_t k) const {
    const T* row = &dist_[i * n_];
    const std::size_t count = n_ - i - 1;
    return detail::smallest_indices_by<Index>(i + 1, count, k, [row](std::size_t x) { return row[x]; });
  }

  // Early-stop target: the exact RNG cycle rank.
  [[nodiscard]] std::size_t rng_early_stop_target() const { return rng_cycle_rank_matrix(*this); }
  [[nodiscard]] Lune_result<T, Index> lune(const Batch_edge<T, Index>& e,
                                           const detail::Edge_map<std::size_t, Index>& one_simp_to_idx,
                                           std::size_t n) const {
    return Lune_builder<T, Index>(e, one_simp_to_idx, n).build_matrix(*this);
  }

 private:
  std::vector<T> dist_;
  std::size_t n_;
};

}  // namespace reduced_rips

}  // namespace Gudhi

#endif  // REDUCED_RIPS_MATRIX_GEOMETRY_H_
