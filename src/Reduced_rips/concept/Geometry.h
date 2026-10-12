/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Thomas Burnett
 *
 *    Copyright (C) 2026 Thomas Burnett, Musashi Koyama
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef CONCEPT_REDUCED_RIPS_GEOMETRY_H_
#define CONCEPT_REDUCED_RIPS_GEOMETRY_H_

namespace Gudhi {

namespace reduced_rips {

/** \brief The concept Geometry describes the metric interface the Reduced_rips core
 * (@ref Gudhi::reduced_rips::Persistence_engine) runs against. Two models are included with the module:
 * @ref Gudhi::reduced_rips::Euclidean_geometry (coordinates and a kd-tree, working in squared distances) and
 * @ref Gudhi::reduced_rips::Matrix_geometry (a bare symmetric distance matrix).
 *
 * The algorithm relies only on the *ordering* of the values returned by `dist`, so a model is free to work in
 * whatever scale it computes most cheaply and exactly. `to_distance` maps that scale back to a true distance
 * for the output barcode.
 */
struct Geometry {
  /** \brief Arithmetic type of the filtration / barcode values. Must be comparable with < and ==. */
  typedef unspecified Filtration_value;

  /** \brief Unsigned integer type used to store point indices and 1-simplex ids (e.g. `std::uint32_t`). It must
   * be wide enough to number the processed 1-simplices, which can far exceed the point count (the engine throws
   * when they no longer fit). The neighbor lists and boundary columns are stored in it. */
  typedef unspecified Index;

  /** \brief Returns the number of points. */
  std::size_t size();

  /** \brief Returns the ordering scale between points i and j: any value monotone in the true distance (for
   * instance the squared distance). Only its ordering is used by the reduction. */
  Filtration_value dist(std::size_t i, std::size_t j);

  /** \brief Maps a value on the ordering scale (as returned by `dist`) back to a true distance, for the output
   * barcode. May be a static member. */
  Filtration_value to_distance(Filtration_value ordering_value);

  /** \brief Returns approximately the k nearest points to i whose index is > i, ascending by distance (ties
   * broken by ascending index). May return fewer than k points (or none), or more when distances tie; the
   * list must be a prefix of the `neighbors_above` ordering, as the reduction resumes positionally in that
   * list after a refresh. */
  std::vector<Index> nearest_neighbors_above(std::size_t i, std::size_t k);

  /** \brief Returns the k nearest points with index > i, ascending by distance from i (ties by ascending
   * index): a prefix of the full above-i ordering. When k reaches the above-i count it is the whole tail. The
   * engine grows this on demand (with a doubling k) when the heap frontier outruns the `nearest_neighbors_above`
   * prefetch, so the returned prefixes for increasing k must be consistent (each a prefix of the next). */
  std::vector<Index> neighbors_above(std::size_t i, std::size_t k);

  /** \brief Returns the relative-neighborhood-graph cycle rank: the number of finite degree-1 bars. This serves
   * as the early-stop target. */
  std::size_t rng_early_stop_target();

  /** \brief Returns the 2-simplices the candidate edge `e` contributes, as a `Lune_result`: either a single
   * apparent 2-simplex, or one boundary column per connected component of the edge's lune. `one_simp_to_idx`
   * maps a packed edge (a 64-bit key) to its 1-simplex id, and `n` is the point count. */
  Lune_result<Filtration_value, Index> lune(const Batch_edge<Filtration_value, Index>& e,
                                            const Edge_map<std::size_t, Index>& one_simp_to_idx, std::size_t n);
};

}  // namespace reduced_rips

}  // namespace Gudhi

#endif  // CONCEPT_REDUCED_RIPS_GEOMETRY_H_
