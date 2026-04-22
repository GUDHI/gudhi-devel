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

#ifndef STEENROD_COCYCLES_H_
#define STEENROD_COCYCLES_H_

#include <algorithm>
#include <numeric>
#include <unordered_set>
#include <utility>
#include <vector>

#include <gudhi/Steenrod_types.h>
#include <gudhi/Steenrod_reduction.h>

namespace Gudhi {

namespace steenrod_persistence {

/** \brief Cohomology barcode together with cocycle representatives.
 * \ingroup steenrod_persistence
 */
struct Barcode_result {
  BarcodeByDim barcode;                           /**< ``barcode[d]`` = bars in dim ``d``. */
  std::vector<std::vector<CohoRep>> coho_reps;    /**< ``coho_reps[d][j]`` = rep for bar ``j``. */
};

/** \brief Return indices that sort bars by *decreasing* birth.
 * \ingroup steenrod_persistence
 */
inline std::vector<int> lexsort_barcode(const std::vector<Bar>& bars) {
  std::vector<int> idx(bars.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::stable_sort(idx.begin(), idx.end(), [&](int a, int b) { return bars[a][1] > bars[b][1]; });
  return idx;
}

/** \brief Apply a permutation to a vector.
 * \ingroup steenrod_persistence
 */
template <typename T>
inline std::vector<T> apply_permutation(const std::vector<T>& v, const std::vector<int>& order) {
  std::vector<T> out(order.size());
  for (std::size_t i = 0; i < order.size(); ++i) out[i] = v[order[i]];
  return out;
}

/** \brief Extract the relative cohomology barcode and cocycle representatives from ``R = D V``.
 *
 * \ingroup steenrod_persistence
 *
 * Bars use the convention ``{death_index, birth_index}``. In relative
 * cohomology the birth index is greater than the death index. A ``death``
 * of ``-1`` marks an essential (infinite) bar.
 *
 * @param[in] rt                 Output of \ref compute_reduced_triangular.
 * @param[in] filtration_values  Optional per-simplex filtration values. When
 *                               provided, finite bars whose birth and death
 *                               simplices share a filtration value (zero
 *                               length) are discarded.
 * @return Barcode together with cocycle representatives.
 */
inline Barcode_result compute_barcode_and_coho_reps(const Reduced_triangular& rt,
                                                    const std::vector<double>* filtration_values = nullptr) {
  const int n_dims = static_cast<int>(rt.idxs.size());

  Barcode_result result;
  result.barcode.resize(n_dims);
  result.coho_reps.resize(n_dims);

  // Dim 0: every bar is essential (no smaller simplices).
  {
    std::vector<Bar> pairs;
    std::vector<CohoRep> reps;

    for (int i = 0; i < static_cast<int>(rt.idxs[0].size()); ++i) {
      if (rt.reduced[0][i].empty()) {
        pairs.push_back({Index{-1}, rt.idxs[0][i]});
        reps.push_back(rt.triangular[0][i]);
      }
    }
    const auto order = lexsort_barcode(pairs);
    result.barcode[0] = apply_permutation(pairs, order);
    result.coho_reps[0] = apply_permutation(reps, order);
  }

  // Dims 1..maxdim.
  for (int dim = 1; dim < n_dims; ++dim) {
    std::vector<Bar> pairs;
    std::vector<CohoRep> reps;
    std::unordered_set<Index> all_birth_idxs;

    // Finite bars: each non-empty column of reduced[dim-1] pairs a (dim-1)-simplex
    // (death) with a dim-simplex (birth).
    for (int i = 0; i < static_cast<int>(rt.idxs[dim - 1].size()); ++i) {
      const Column& col = rt.reduced[dim - 1][i];
      if (!col.empty()) {
        const Index birth_rel = col.front();
        const Index b = rt.idxs[dim][birth_rel];
        const Index d = rt.idxs[dim - 1][i];

        // Mark birth as claimed even if the bar is filtered out.
        all_birth_idxs.insert(b);

        bool include = true;
        if (filtration_values != nullptr) {
          include = ((*filtration_values)[b] != (*filtration_values)[d]);
        }

        if (include) {
          pairs.push_back({d, b});
          reps.push_back(col);
        }
      }
    }

    // Essential bars: unpaired dim-simplices with an empty reduced column.
    for (int i = 0; i < static_cast<int>(rt.idxs[dim].size()); ++i) {
      const Index abs_idx = rt.idxs[dim][i];
      if (all_birth_idxs.count(abs_idx) == 0 && rt.reduced[dim][i].empty()) {
        pairs.push_back({Index{-1}, abs_idx});
        reps.push_back(rt.triangular[dim][i]);
      }
    }

    const auto order = lexsort_barcode(pairs);
    result.barcode[dim] = apply_permutation(pairs, order);
    result.coho_reps[dim] = apply_permutation(reps, order);
  }

  return result;
}

}  // namespace steenrod_persistence

}  // namespace Gudhi

#endif  // STEENROD_COCYCLES_H_
