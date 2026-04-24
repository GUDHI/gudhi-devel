/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Anibal M. Medina-Mardones
 *
 *    Copyright (C) 2026 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef STEENROD_SQ_H_
#define STEENROD_SQ_H_

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <gudhi/Steenrod_types.h>
#include <gudhi/Steenrod_cocycles.h>
#include <gudhi/Steenrod_reduction.h>

namespace Gudhi {

namespace steenrod_persistence {

/** \brief Sq^k matrix by dimension.
 * \ingroup steenrod_persistence
 *
 * ``Steenrod_matrix[d][j]`` is a sorted column of relative (dim-``d``) indices
 * representing Sq^k applied to the ``j``-th cohomology representative in
 * degree ``d-k``. Indices ``0..k-1`` are empty padding (Sq^k has no source
 * in those degrees).
 */
using Steenrod_matrix = std::vector<Matrix>;

/** \brief Decide whether the simplex pair ``(a, b)`` contributes to Sq^k.
 * \ingroup steenrod_persistence
 *
 * Preconditions (guaranteed by the caller):
 *   - ``a``, ``b``, ``u`` are sorted and ``u`` equals ``sorted(a ∪ b)``.
 *   - ``|u| == dim + k + 1`` (caller already checked).
 *   - ``u`` exists in the target dimension's simplex-index map (caller checked).
 */
inline bool stsq_parity(const Simplex& a, const Simplex& b, const Simplex& u) {
  // Asymmetric differences.
  Simplex a_bar;
  Simplex b_bar;
  a_bar.reserve(a.size());
  b_bar.reserve(b.size());
  std::set_difference(a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(a_bar));
  std::set_difference(b.begin(), b.end(), a.begin(), a.end(), std::back_inserter(b_bar));

  // u_bar = merge of two disjoint sorted ranges.
  Simplex u_bar;
  u_bar.reserve(a_bar.size() + b_bar.size());
  std::merge(a_bar.begin(), a_bar.end(), b_bar.begin(), b_bar.end(), std::back_inserter(u_bar));

  // All vertices in a_bar must share the same (pos_u + pos_bar) parity.
  int a_idx = -1;
  for (Index v : a_bar) {
    const int pos_u = static_cast<int>(std::lower_bound(u.begin(), u.end(), v) - u.begin());
    const int pos_bar = static_cast<int>(std::lower_bound(u_bar.begin(), u_bar.end(), v) - u_bar.begin());
    const int idx_v = (pos_u + pos_bar) & 1;
    if (a_idx < 0) {
      a_idx = idx_v;
    } else if (idx_v != a_idx) {
      return false;
    }
  }

  // All vertices in b_bar must have the complementary parity.
  const int b_expected = 1 - a_idx;
  for (Index v : b_bar) {
    const int pos_u = static_cast<int>(std::lower_bound(u.begin(), u.end(), v) - u.begin());
    const int pos_bar = static_cast<int>(std::lower_bound(u_bar.begin(), u_bar.end(), v) - u_bar.begin());
    const int idx_v = (pos_u + pos_bar) & 1;
    if (idx_v != b_expected) return false;
  }

  return true;
}

/** \brief Apply Sq^k to every cohomology representative in one dimension.
 * \ingroup steenrod_persistence
 *
 * @param[in] coho_reps_dim Cohomology representatives in this dimension.
 * @param[in] tups_dim      Vertex lists for every dim-simplex.
 * @param[in] spx2idx_dpk   Map from (dim+k)-simplex to its relative index.
 * @param[in] length        ``dim + k + 1``.
 * @param[in] n_jobs        Number of OpenMP threads; ``0`` or negative = all.
 * @return One column per representative (sorted relative indices in dim+k).
 */
inline Matrix populate_steenrod_matrix_single_dim(const std::vector<CohoRep>& coho_reps_dim,
                                                  const std::vector<Simplex>& tups_dim,
                                                  const Spx_to_idx& spx2idx_dpk,
                                                  int length,
                                                  int n_jobs = -1) {
  const int n_reps = static_cast<int>(coho_reps_dim.size());
#ifdef _OPENMP
  const int n_threads = (n_jobs <= 0) ? omp_get_max_threads() : n_jobs;
#else
  (void)n_jobs;
#endif

  Matrix result(n_reps);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(n_threads)
#endif
  for (int idx = 0; idx < n_reps; ++idx) {
    const CohoRep& rep = coho_reps_dim[idx];

    // Gather actual simplex vertex lists for this representative.
    std::vector<Simplex> cocycle;
    cocycle.reserve(rep.size());
    for (Index rel : rep) cocycle.push_back(tups_dim[rel]);

    // Accumulate (dim+k)-relative indices in a flat vector; each occurrence
    // toggles GF(2) presence, reduced by run-length-mod-2 after sorting.
    std::vector<Index> flat_indices;
    Simplex u;
    u.reserve(length);

    const int nc = static_cast<int>(cocycle.size());
    for (int i = 0; i < nc; ++i) {
      for (int j = i + 1; j < nc; ++j) {
        u.clear();
        std::set_union(cocycle[i].begin(), cocycle[i].end(),
                       cocycle[j].begin(), cocycle[j].end(),
                       std::back_inserter(u));

        if (static_cast<int>(u.size()) != length) continue;

        auto it = spx2idx_dpk.find(u);
        if (it == spx2idx_dpk.end()) continue;
        if (!stsq_parity(cocycle[i], cocycle[j], u)) continue;

        flat_indices.push_back(it->second);
      }
    }

    // Reduce mod 2: keep indices appearing an odd number of times.
    std::sort(flat_indices.begin(), flat_indices.end());
    Column col;
    col.reserve(flat_indices.size());
    const std::size_t n_flat = flat_indices.size();
    for (std::size_t p = 0; p < n_flat;) {
      std::size_t q = p + 1;
      while (q < n_flat && flat_indices[q] == flat_indices[p]) ++q;
      if ((q - p) & 1) col.push_back(flat_indices[p]);
      p = q;
    }
    result[idx] = std::move(col);
  }

  return result;
}

/** \brief Compute the Sq^k Steenrod matrices across all dimensions.
 * \ingroup steenrod_persistence
 */
inline Steenrod_matrix compute_steenrod_matrix(int k,
                                               const std::vector<std::vector<CohoRep>>& coho_reps,
                                               const Filtration_by_dim& filtration_by_dim,
                                               const std::vector<Spx_to_idx>& spx2idx,
                                               int n_jobs = -1) {
  const int n_dims = static_cast<int>(coho_reps.size());
  Steenrod_matrix st_matrix(k);  // k empty padding matrices

  for (int dim = 0; dim < n_dims - k; ++dim) {
    const int dim_plus_k = dim + k;
    const int length = dim_plus_k + 1;
    st_matrix.push_back(populate_steenrod_matrix_single_dim(coho_reps[dim],
                                                            filtration_by_dim[dim].tups,
                                                            spx2idx[dim_plus_k],
                                                            length,
                                                            n_jobs));
  }

  return st_matrix;
}

}  // namespace steenrod_persistence

}  // namespace Gudhi

#endif  // STEENROD_SQ_H_
