/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Anibal M. Medina-Mardones
 *
 *    Copyright (C) 2026 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef STEENROD_BARCODE_IMPL_H_
#define STEENROD_BARCODE_IMPL_H_

#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <gudhi/Steenrod_types.h>
#include <gudhi/Steenrod_utils.h>
#include <gudhi/Steenrod_reduction.h>
#include <gudhi/Steenrod_cocycles.h>
#include <gudhi/Steenrod_sq.h>

namespace Gudhi {

namespace steenrod_persistence {

/** \brief Compute the Sq^k barcode for one cohomological degree.
 * \ingroup steenrod_persistence
 *
 * Augmented column reduction. The augmented matrix is the concatenation of
 * ``reduced_prev_dim`` and ``steenrod_matrix_dim``. Process (dim-1)-simplices
 * in reverse filtration order, registering their R-column pivots permanently.
 * After each (dim-1)-simplex, reduce any newly unlocked Steenrod column (those
 * with birth ≤ the current absolute index) against the permanent pivots.
 * A column that reduces to zero records a finite bar; surviving columns yield
 * essential bars.
 *
 * @param[in] steenrod_matrix_dim   Steenrod columns in dimension ``dim``.
 * @param[in] n_idxs_dim            Number of dim-simplices.
 * @param[in] idxs_prev_dim         Absolute indices of (dim-1)-simplices, sorted.
 * @param[in] reduced_prev_dim      R-matrix columns from dimension ``dim-1``.
 * @param[in] births_dim_minus_k    Birth indices of (dim-k)-classes, decreasing.
 * @return Barcode of finite + essential bars. ``death == -1`` marks essentials.
 */
inline Barcode steenrod_barcode_single_dim(const Matrix& steenrod_matrix_dim,
                                           int n_idxs_dim,
                                           const std::vector<Index>& idxs_prev_dim,
                                           const Matrix& reduced_prev_dim,
                                           const std::vector<Index>& births_dim_minus_k) {
  const int n_cols_red = static_cast<int>(reduced_prev_dim.size());
  const int n_cols_st = static_cast<int>(steenrod_matrix_dim.size());
  const int n_births = static_cast<int>(births_dim_minus_k.size());

  // Augmented = [reduced_prev_dim | steenrod_matrix_dim]
  std::vector<Column> augmented;
  augmented.reserve(n_cols_red + n_cols_st);
  for (const Column& c : reduced_prev_dim) augmented.push_back(c);
  for (const Column& c : steenrod_matrix_dim) augmented.push_back(c);

  std::vector<int> pivots_lookup(n_idxs_dim, -1);
  std::vector<bool> alive(n_cols_st, true);
  Barcode result;

  if (n_births != 0 && n_cols_red != 0) {
    int n_cols_st_curr = 0;
    {
      const Index max_prev = idxs_prev_dim.back();
      for (int ii = 0; ii < n_births; ++ii) {
        n_cols_st_curr = ii;
        if (max_prev >= births_dim_minus_k[ii]) break;
      }
      // If the loop ran to completion, n_cols_st_curr stays at n_births - 1
      // (mirrors Python ``enumerate`` semantics).
    }

    for (int i = n_cols_red - 1; i >= 0; --i) {
      const Index idx = idxs_prev_dim[i];

      // Unlock Steenrod columns whose birth lies between idxs_prev_dim[i-1]
      // (exclusive) and idx (inclusive).
      if (n_cols_st_curr < n_cols_st) {
        const Index next_birth = births_dim_minus_k[n_cols_st_curr];
        if (idx >= next_birth) {
          if (i > 0) {
            const Index next_idx = idxs_prev_dim[i - 1];
            for (int jj = n_cols_st_curr; jj < n_cols_st; ++jj) {
              if (next_idx >= births_dim_minus_k[jj]) break;
              ++n_cols_st_curr;
            }
          } else {
            n_cols_st_curr = n_cols_st;
          }
        }
      }

      // Register the R-column pivot of this (dim-1)-simplex permanently.
      if (!augmented[i].empty()) {
        pivots_lookup[augmented[i].front()] = i;
      }

      // Reduce each unlocked Steenrod column against the permanent pivots.
      std::vector<Index> st_pivots_claimed;
      for (int ii = n_cols_red; ii < n_cols_red + n_cols_st_curr; ++ii) {
        Index highest_one = augmented[ii].empty() ? -1 : augmented[ii].front();
        int pc = (highest_one >= 0) ? pivots_lookup[highest_one] : -1;

        while (highest_one >= 0 && pc >= 0) {
          symm_diff_inplace(augmented[ii], augmented[pc]);
          highest_one = augmented[ii].empty() ? -1 : augmented[ii].front();
          pc = (highest_one >= 0) ? pivots_lookup[highest_one] : -1;
        }

        if (highest_one >= 0) {
          // Column did not reduce to zero: claim pivot temporarily.
          pivots_lookup[highest_one] = ii;
          st_pivots_claimed.push_back(highest_one);
        } else {
          // Column reduced to zero: the class is killed here.
          const int st_idx = ii - n_cols_red;
          if (alive[st_idx]) {
            alive[st_idx] = false;
            const Index birth = births_dim_minus_k[st_idx];
            if (idx < birth) result.push_back({idx, birth});
          }
        }
      }

      // Reset pivots claimed temporarily by Steenrod columns.
      for (Index p : st_pivots_claimed) pivots_lookup[p] = -1;
    }
  }

  // Surviving alive Steenrod classes yield essential bars.
  for (int i = 0; i < n_cols_st; ++i) {
    if (alive[i]) result.push_back({Index{-1}, births_dim_minus_k[i]});
  }

  return result;
}

/** \brief Compute the Sq^k persistence barcode from the Steenrod matrix.
 * \ingroup steenrod_persistence
 *
 * Each per-dimension slot is independent: the outer loop is parallelised
 * with OpenMP.
 *
 * @param[in] k                  Steenrod square exponent.
 * @param[in] steenrod_matrix    Output of \ref compute_steenrod_matrix.
 * @param[in] rt                 Output of \ref compute_reduced_triangular.
 * @param[in] br                 Output of \ref compute_barcode_and_coho_reps.
 * @param[in] filtration_values  Optional filtration values; when provided,
 *                               zero-length finite bars are discarded.
 * @param[in] n_jobs             OpenMP threads; ``0`` or negative = all.
 * @return ``BarcodeByDim`` where indices ``0..k-1`` are empty.
 */
inline BarcodeByDim compute_steenrod_barcode(int k,
                                             const Steenrod_matrix& steenrod_matrix,
                                             const Reduced_triangular& rt,
                                             const Barcode_result& br,
                                             const std::vector<double>* filtration_values = nullptr,
                                             int n_jobs = -1) {
  const int n_st_dims = static_cast<int>(steenrod_matrix.size());
  BarcodeByDim sb(n_st_dims);

#ifdef _OPENMP
  const int n_threads = (n_jobs <= 0) ? omp_get_max_threads() : n_jobs;
#else
  (void)n_jobs;
#endif

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(n_threads)
#endif
  for (int dim = k; dim < n_st_dims; ++dim) {
    const Barcode& barcode_src = br.barcode[dim - k];
    std::vector<Index> births;
    births.reserve(barcode_src.size());
    for (const Bar& bar : barcode_src) births.push_back(bar[1]);

    const bool has_prev = (dim - 1 < static_cast<int>(rt.idxs.size()));
    const bool has_dim = (dim < static_cast<int>(rt.idxs.size()));

    Barcode sbd;
    if (has_prev && has_dim && !births.empty()) {
      sbd = steenrod_barcode_single_dim(steenrod_matrix[dim],
                                        static_cast<int>(rt.idxs[dim].size()),
                                        rt.idxs[dim - 1],
                                        rt.reduced[dim - 1],
                                        births);
    } else {
      for (const Index b : births) sbd.push_back({Index{-1}, b});
    }

    if (filtration_values != nullptr) {
      Barcode filtered;
      filtered.reserve(sbd.size());
      for (const Bar& bar : sbd) {
        if (bar[0] == -1) {
          filtered.push_back(bar);
        } else if ((*filtration_values)[bar[0]] != (*filtration_values)[bar[1]]) {
          filtered.push_back(bar);
        }
      }
      sbd = std::move(filtered);
    }

    sb[dim] = std::move(sbd);
  }

  return sb;
}

}  // namespace steenrod_persistence

}  // namespace Gudhi

#endif  // STEENROD_BARCODE_IMPL_H_
