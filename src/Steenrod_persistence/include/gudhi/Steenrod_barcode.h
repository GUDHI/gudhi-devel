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

#include <gudhi/Debug_utils.h>
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
  // Debug-mode guard against negative k.  Python users already get a clean
  // ValueError from the nanobind interface in any build mode; this guard is
  // the gudhi house-style way to flag direct C++ misuse in Debug builds.
  GUDHI_CHECK(k >= 0, std::invalid_argument("Sq^k exponent k must be non-negative"));
  // Pass ``maxdim`` down to the reduction as the loop bound; no copy of the
  // filtration is needed (compute_reduced_triangular clamps negative or
  // out-of-range values to fbd.size() - 1, so passing -1 keeps the
  // "process all dimensions" default).  All downstream functions derive
  // their own ``n_dims`` from ``rt.idxs.size()`` (or its descendants
  // ``br.coho_reps.size()`` and ``sm.size()``), so the truncation
  // propagates naturally without further plumbing.
  auto rt = compute_reduced_triangular(fbd, maxdim);
  auto br = compute_barcode_and_coho_reps(rt, filtration_values);

  // Sq^0 is the identity, so the Steenrod barcode is the ordinary barcode.
  // Short-circuit before compute_steenrod_matrix / compute_steenrod_barcode:
  // those routines loop ``for (int dim = k; ...)`` and read ``rt.idxs[dim - 1]``,
  // which underflows the bounds check when ``k == 0``.
  if (k == 0) {
    auto barcode_copy = br.barcode;
    return Barcodes_result{std::move(br.barcode), std::move(barcode_copy)};
  }

  auto sm = compute_steenrod_matrix(k, br.coho_reps, fbd, rt.spx2idx, n_jobs);
  auto sb = compute_steenrod_barcode(k, sm, rt, br, filtration_values, n_jobs);

  return Barcodes_result{std::move(br.barcode), std::move(sb)};
}

}  // namespace steenrod_persistence

}  // namespace Gudhi

#endif  // STEENROD_BARCODE_H_
