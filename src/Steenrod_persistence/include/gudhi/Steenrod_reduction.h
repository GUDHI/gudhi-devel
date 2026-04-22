/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Umberto Lupo, Anibal M. Medina-Mardones, Guillaume Tauzin
 *                     (port to Gudhi: Maximiliano Alvarez)
 *
 *    Copyright (C) 2026 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef STEENROD_REDUCTION_H_
#define STEENROD_REDUCTION_H_

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <functional>
#include <unordered_map>
#include <utility>
#include <vector>

#include <gudhi/Steenrod_types.h>
#include <gudhi/Steenrod_utils.h>

namespace Gudhi {

namespace steenrod_persistence {

/** \brief A matrix: a list of columns, each a sorted ascending ``Column``.
 * \ingroup steenrod_persistence
 */
using Matrix = std::vector<Column>;

/** \brief Hash functor for a sorted vertex list.
 * \ingroup steenrod_persistence
 */
struct Simplex_hash {
  std::size_t operator()(const Simplex& s) const noexcept {
    // Boost-style combining hash.
    std::size_t seed = s.size();
    for (Index x : s) {
      seed ^= std::hash<Index>{}(x) + 0x9e3779b9u + (seed << 6) + (seed >> 2);
    }
    return seed;
  }
};

/** \brief Map from a simplex (sorted vertex list) to its relative index within a dimension.
 * \ingroup steenrod_persistence
 */
using Spx_to_idx = std::unordered_map<Simplex, Index, Simplex_hash>;

/** \brief All data for one dimension of the filtration. \ingroup steenrod_persistence */
struct Dim_filtration {
  std::vector<Index> idxs;     /**< Absolute filtration positions of simplices. */
  std::vector<Simplex> tups;   /**< Sorted vertex lists, one per simplex. */
};

/** \brief Full filtration grouped by dimension: ``filtration_by_dim[d]`` is the dim-``d`` data.
 * \ingroup steenrod_persistence
 */
using Filtration_by_dim = std::vector<Dim_filtration>;

/** \brief Output of \ref compute_reduced_triangular. \ingroup steenrod_persistence */
struct Reduced_triangular {
  std::vector<Spx_to_idx> spx2idx;            /**< Simplex-to-relative-index maps, per dim. */
  std::vector<std::vector<Index>> idxs;       /**< Absolute indices, per dim. */
  std::vector<Matrix> reduced;                /**< R-matrix columns, per dim. */
  std::vector<Matrix> triangular;             /**< V-matrix columns, per dim. */
};

/** \brief Return the codimension-1 faces of ``spx`` (each sorted).
 * \ingroup steenrod_persistence
 */
inline std::vector<Simplex> drop_elements(const Simplex& spx) {
  const std::size_t n = spx.size();
  std::vector<Simplex> faces;
  faces.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    Simplex face;
    face.reserve(n - 1);
    for (std::size_t j = 0; j < n; ++j) {
      if (j != i) face.push_back(spx[j]);
    }
    // Face is sorted because spx is sorted and we skipped exactly one element.
    faces.push_back(std::move(face));
  }
  return faces;
}

/** \brief Persistent cohomology matrix reduction with the clearing optimization.
 *
 * \ingroup steenrod_persistence
 *
 * Columns are stored in *ascending* order; the pivot of a column is its first
 * (lowest) element. Both ``coboundary`` and ``triangular`` are modified in
 * place, preserving the invariant ``R = D V``.
 *
 * @param[in,out] coboundary     Matrix ``D`` on entry, reduced matrix ``R`` on exit.
 * @param[in,out] triangular     Triangular matrix ``V``, updated in place.
 * @param[in,out] pivots_lookup  ``pivots_lookup[row] = col`` — column owning
 *                               this row as its pivot, or ``-1``. Indexed by
 *                               relative row in the *next* dimension.
 * @return Relative indices (in the next dimension) that were claimed as
 *         pivots — these will be cleared in the next-dimension coboundary.
 */
inline std::vector<Index> twist_reduction(Matrix& coboundary,
                                          Matrix& triangular,
                                          std::vector<Index>& pivots_lookup) {
  const int n = static_cast<int>(coboundary.size());
  std::vector<Index> rel_idxs_to_clear;

  for (int j = n - 1; j >= 0; --j) {
    Index lowest = pivot(coboundary[j]);
    Index pivot_col = (lowest >= 0) ? pivots_lookup[lowest] : -1;

    while (lowest >= 0 && pivot_col >= 0) {
      // GF(2) column addition; both columns share ``lowest`` as first
      // element, so it cancels and ascending order is preserved.
      symm_diff_inplace(coboundary[j], coboundary[pivot_col]);
      symm_diff_inplace(triangular[j], triangular[pivot_col]);
      lowest = pivot(coboundary[j]);
      pivot_col = (lowest >= 0) ? pivots_lookup[lowest] : -1;
    }

    if (lowest >= 0) {
      pivots_lookup[lowest] = j;
      rel_idxs_to_clear.push_back(lowest);
    }
  }

  return rel_idxs_to_clear;
}

/** \brief Bundles every output of a single-dimension reduction. \ingroup steenrod_persistence */
struct Single_dim_reduction {
  Spx_to_idx spx2idx;
  Matrix reduced;
  Matrix triangular;
  std::vector<Index> rel_idxs_to_clear;   /**< Pivot rows — to clear in next dim. */
  std::vector<Index> pivots_lookup;       /**< Needed by \ref fix_triangular_after_clearing. */
};

/** \brief Build and reduce the coboundary matrix for one dimension.
 *
 * \ingroup steenrod_persistence
 *
 * @param[in] idxs_dim              Absolute filtration positions for this dimension.
 * @param[in] tups_dim              Vertex lists for this dimension.
 * @param[in] rel_idxs_to_clear_in  Relative column indices to pre-zero before reduction.
 * @param[in] idxs_next_dim         Absolute filtration positions for ``dim+1`` (optional).
 * @param[in] tups_next_dim         Vertex lists for ``dim+1`` (optional).
 * @return Reduction output for this dimension.
 */
inline Single_dim_reduction reduce_single_dim(const std::vector<Index>& idxs_dim,
                                              const std::vector<Simplex>& tups_dim,
                                              const std::vector<Index>& rel_idxs_to_clear_in,
                                              const std::vector<Index>* idxs_next_dim = nullptr,
                                              const std::vector<Simplex>* tups_next_dim = nullptr) {
  const int n = static_cast<int>(idxs_dim.size());

  // 1. Simplex-to-relative-index lookup for this dimension.
  Spx_to_idx spx2idx;
  spx2idx.reserve(n);
  for (int i = 0; i < n; ++i) {
    spx2idx[tups_dim[i]] = static_cast<Index>(i);
  }

  // 2. Initialise reduced (coboundary) and triangular (V = identity start).
  Matrix reduced(n);
  Matrix triangular(n);
  for (int i = 0; i < n; ++i) {
    triangular[i] = {static_cast<Index>(i)};
  }

  std::vector<Index> pivots_lookup;
  std::vector<Index> rel_idxs_to_clear_out;

  if (idxs_next_dim != nullptr && tups_next_dim != nullptr) {
    const int n_next = static_cast<int>(idxs_next_dim->size());

    // 3. Populate coboundary: for each (d+1)-simplex j, append j to each
    //    d-face's column. Iterating j in ascending order keeps each column
    //    sorted without an extra pass.
    for (int j = 0; j < n_next; ++j) {
      for (const Simplex& face : drop_elements((*tups_next_dim)[j])) {
        auto it = spx2idx.find(face);
        if (it != spx2idx.end()) {
          reduced[it->second].push_back(static_cast<Index>(j));
        }
      }
    }

    // 4. Apply clearing.
    for (Index rel_idx : rel_idxs_to_clear_in) {
      reduced[rel_idx].clear();
    }

    // 5. Reduce.
    pivots_lookup.assign(n_next, -1);
    rel_idxs_to_clear_out = twist_reduction(reduced, triangular, pivots_lookup);
  }
  // If no next dimension: coboundary is all-zero; no reduction needed.

  return Single_dim_reduction{std::move(spx2idx),
                              std::move(reduced),
                              std::move(triangular),
                              std::move(rel_idxs_to_clear_out),
                              std::move(pivots_lookup)};
}

/** \brief Restore ``R = D V`` after the clearing optimisation.
 *
 * \ingroup steenrod_persistence
 *
 * See \cite Chen11persistenthomology §3.2. For each column cleared in the current
 * dimension, sets
 * ``triangular[rel_idx] = reduced_prev[pivots_lookup_prev[rel_idx]]``.
 */
inline void fix_triangular_after_clearing(Matrix& triangular,
                                          const Matrix& reduced_prev,
                                          const std::vector<Index>& rel_idxs_to_clear,
                                          const std::vector<Index>& pivots_lookup_prev) {
  for (Index rel_idx : rel_idxs_to_clear) {
    const Index prev_col = pivots_lookup_prev[rel_idx];
    if (prev_col >= 0) {
      triangular[rel_idx] = reduced_prev[prev_col];
    }
  }
}

/** \brief Full ``R = D V`` decomposition over all dimensions.
 * \ingroup steenrod_persistence
 */
inline Reduced_triangular compute_reduced_triangular(const Filtration_by_dim& filtration_by_dim) {
  const int maxdim = static_cast<int>(filtration_by_dim.size()) - 1;
  assert(maxdim >= 0 && "filtration_by_dim must have at least one dimension");

  Reduced_triangular result;
  result.spx2idx.resize(maxdim + 1);
  result.idxs.resize(maxdim + 1);
  result.reduced.resize(maxdim + 1);
  result.triangular.resize(maxdim + 1);

  std::vector<Index> rel_idxs_to_clear;   // empty at dim 0
  std::vector<Index> pivots_lookup_prev;  // from previous dim's reduction

  for (int dim = 0; dim <= maxdim; ++dim) {
    const auto& df = filtration_by_dim[dim];

    const std::vector<Index>* idxs_next = nullptr;
    const std::vector<Simplex>* tups_next = nullptr;
    if (dim < maxdim) {
      idxs_next = &filtration_by_dim[dim + 1].idxs;
      tups_next = &filtration_by_dim[dim + 1].tups;
    }

    auto r = reduce_single_dim(df.idxs, df.tups, rel_idxs_to_clear, idxs_next, tups_next);

    // At dim == 0 there is nothing to clear, so this is a no-op.
    if (dim > 0) {
      fix_triangular_after_clearing(r.triangular, result.reduced[dim - 1], rel_idxs_to_clear, pivots_lookup_prev);
    }

    result.spx2idx[dim] = std::move(r.spx2idx);
    result.idxs[dim] = df.idxs;
    result.reduced[dim] = std::move(r.reduced);
    result.triangular[dim] = std::move(r.triangular);

    rel_idxs_to_clear = std::move(r.rel_idxs_to_clear);
    pivots_lookup_prev = std::move(r.pivots_lookup);
  }

  return result;
}

}  // namespace steenrod_persistence

}  // namespace Gudhi

#endif  // STEENROD_REDUCTION_H_
