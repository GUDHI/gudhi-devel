/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Anibal M. Medina-Mardones
 *
 *    Copyright (C) 2026 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef STEENROD_BARCODE_H_
#define STEENROD_BARCODE_H_

#include <stdexcept>
#include <utility>
#include <vector>

#include <gudhi/Steenrod_types.h>
#include <gudhi/Steenrod_reduction.h>
#include <gudhi/Steenrod_cocycles.h>
#include <gudhi/Steenrod_sq.h>
#include <gudhi/Steenrod_barcode_impl.h>

namespace Gudhi {

namespace steenrod_persistence {

/** \brief Ordinary persistence barcode + Sq^k Steenrod barcode over GF(2).
 * \ingroup steenrod_persistence
 */
struct Barcodes_result {
  BarcodeByDim ordinary;   /**< Ordinary persistence barcode, per degree. */
  BarcodeByDim steenrod;   /**< Sq^k Steenrod barcode,        per degree. */
};

/** \brief Run the full pipeline and return both barcodes.
 *
 * \ingroup steenrod_persistence
 *
 * Coefficients are in GF(2) — Steenrod squares exist only over \f$\mathbb{F}_2\f$.
 *
 * @param[in] k                  Steenrod squaring operation Sq^k (non-negative).
 * @param[in] fbd                Filtration grouped by dimension.
 * @param[in] filtration_values  Optional filtration values. When provided,
 *                               zero-length bars (birth and death at the
 *                               same value) are discarded from both barcodes.
 * @param[in] n_jobs             OpenMP thread count for the parallelised
 *                               stages (compute_steenrod_matrix and
 *                               compute_steenrod_barcode). ``-1`` = all cores.
 * @param[in] maxdim             Highest cohomological degree at which to
 *                               compute the Steenrod barcode. Dimensions
 *                               above ``maxdim`` are dropped from ``fbd``
 *                               before the pipeline runs. ``-1`` keeps the
 *                               full filtration. Note that
 *                               ``ordinary[maxdim]`` is incomplete when
 *                               ``maxdim`` is used: some classes that should
 *                               die at ``maxdim`` remain essential because
 *                               their coboundary simplices were dropped.
 *                               The Steenrod barcode through ``maxdim`` is
 *                               correct.
 * @return ``Barcodes_result`` with the ordinary and Steenrod barcodes.
 */
inline Barcodes_result barcodes(int k,
                                const Filtration_by_dim& fbd,
                                const std::vector<double>* filtration_values = nullptr,
                                int n_jobs = -1,
                                int maxdim = -1) {
  if (k < 0) {
    throw std::invalid_argument("Sq^k exponent k must be non-negative");
  }
  const int n_dims_in = static_cast<int>(fbd.size());
  const bool truncate = (maxdim >= 0) && (maxdim + 1 < n_dims_in);

  Filtration_by_dim fbd_capped;
  const Filtration_by_dim* fbd_used = &fbd;
  if (truncate) {
    fbd_capped.reserve(maxdim + 1);
    for (int d = 0; d <= maxdim; ++d) fbd_capped.push_back(fbd[d]);
    fbd_used = &fbd_capped;
  }

  auto rt = compute_reduced_triangular(*fbd_used);
  auto br = compute_barcode_and_coho_reps(rt, filtration_values);

  // Sq^0 is the identity, so the Steenrod barcode is the ordinary barcode.
  // Short-circuit before compute_steenrod_matrix / compute_steenrod_barcode:
  // those routines loop ``for (int dim = k; ...)`` and read ``rt.idxs[dim - 1]``,
  // which underflows the bounds check when ``k == 0``.
  if (k == 0) {
    auto barcode_copy = br.barcode;
    return Barcodes_result{std::move(br.barcode), std::move(barcode_copy)};
  }

  auto sm = compute_steenrod_matrix(k, br.coho_reps, *fbd_used, rt.spx2idx, n_jobs);
  auto sb = compute_steenrod_barcode(k, sm, rt, br, filtration_values, n_jobs);

  return Barcodes_result{std::move(br.barcode), std::move(sb)};
}

}  // namespace steenrod_persistence

}  // namespace Gudhi

#endif  // STEENROD_BARCODE_H_
