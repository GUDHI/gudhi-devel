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

#include <cstddef>
#include <ostream>
#include <stdexcept>
#include <utility>
#include <vector>

#ifdef GUDHI_USE_TBB
#include <tbb/task_arena.h>
#endif

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

  /** \brief Print bars-per-dimension summary for both barcodes. */
  friend std::ostream& operator<<(std::ostream& os, const Barcodes_result& r) {
    os << "ordinary barcode:\n";
    for (std::size_t d = 0; d < r.ordinary.size(); ++d) {
      os << "  dim " << d << " : " << r.ordinary[d].size() << " bar(s)\n";
    }
    os << "steenrod barcode:\n";
    for (std::size_t d = 0; d < r.steenrod.size(); ++d) {
      os << "  dim " << d << " : " << r.steenrod[d].size() << " bar(s)\n";
    }
    return os;
  }
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
 * @param[in] n_jobs             Maximum TBB worker threads for the
 *                               parallelised stages.  ``-1`` (default) uses
 *                               the global TBB scheduler default.  Honoured
 *                               via ``tbb::task_arena`` only when the build
 *                               has ``GUDHI_USE_TBB`` defined.
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

  // Empty filtration: every downstream stage assumes at least one
  // dimension exists, and would dereference rt.idxs[0] in
  // compute_barcode_and_coho_reps.  Return empty barcodes instead.
  if (fbd.empty()) return Barcodes_result{};

  // Body of the reduction.  When TBB is enabled we run it inside a
  // ``tbb::task_arena`` so that ``n_jobs`` limits the concurrency for this
  // call.  ``tbb::task_arena::automatic`` (= -1) means "use the global
  // default" — the same behaviour as a plain ``tbb::parallel_for`` outside
  // any explicit arena.  Inner routines
  // (``populate_steenrod_matrix_single_dim``, ``compute_steenrod_barcode``)
  // use ``tbb::parallel_for`` and inherit this arena.  When TBB is not
  // enabled, ``n_jobs`` is ignored and the inner loops run serially.
  auto body = [&]() -> Barcodes_result {
    // Pass ``maxdim`` down to the reduction as the loop bound; no copy of
    // the filtration is needed (compute_reduced_triangular clamps
    // negative or out-of-range values to fbd.size() - 1).
    auto rt = compute_reduced_triangular(fbd, maxdim);
    auto br = compute_barcode_and_coho_reps(rt, filtration_values);

    // Sq^0 is the identity, so the Steenrod barcode is the ordinary one.
    // Short-circuit before compute_steenrod_matrix / compute_steenrod_barcode:
    // those routines loop ``for (int dim = k; ...)`` and read
    // ``rt.idxs[dim - 1]``, which underflows the bounds check when
    // ``k == 0``.
    if (k == 0) {
      auto barcode_copy = br.barcode;
      return Barcodes_result{std::move(br.barcode), std::move(barcode_copy)};
    }

    auto sm = compute_steenrod_matrix(k, br.coho_reps, fbd, rt.spx2idx, n_jobs);
    auto sb = compute_steenrod_barcode(k, sm, rt, br, filtration_values, n_jobs);
    return Barcodes_result{std::move(br.barcode), std::move(sb)};
  };

#ifdef GUDHI_USE_TBB
  Barcodes_result result;
  tbb::task_arena arena(n_jobs > 0 ? n_jobs : tbb::task_arena::automatic);
  arena.execute([&] { result = body(); });
  return result;
#else
  return body();
#endif
}

}  // namespace steenrod_persistence

}  // namespace Gudhi

#endif  // STEENROD_BARCODE_H_
