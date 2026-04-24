/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Anibal M. Medina-Mardones
 *
 *    Copyright (C) 2026 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef INCLUDE_STEENROD_BARCODE_INTERFACE_H_
#define INCLUDE_STEENROD_BARCODE_INTERFACE_H_

#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include <nanobind/nanobind.h>

#include <gudhi/Steenrod_barcode.h>

namespace Gudhi {

/** \brief Python-interface wrapper that exposes Steenrod persistence barcodes
 *         to the nanobind SimplexTree bindings.
 *
 * Usage mirrors ``Persistent_cohomology_interface``:
 *   1. Construct with a pointer to a ``Simplex_tree_interface`` and the Sq^k exponent.
 *   2. Call ``compute()``; receive a pair ``(ordinary, steenrod)`` of
 *      ``vector<vector<pair<double, double>>>`` with one inner vector per
 *      homological dimension.
 *
 * ``death = +infinity`` marks essential bars. Coefficients are GF(2); no field
 * parameter is exposed.
 */
template <class Simplex_tree_t>
class Steenrod_barcode_interface {
 public:
  using Filtration_value = typename Simplex_tree_t::Filtration_value;
  using Birth_death = std::pair<double, double>;
  using Bars_one_dim = std::vector<Birth_death>;
  using Bars_by_dim = std::vector<Bars_one_dim>;

  Steenrod_barcode_interface(Simplex_tree_t* stptr, int k) : stptr_(stptr), k_(k) {
    if (k < 0) throw std::invalid_argument("Sq^k exponent k must be non-negative");
  }

  /** Returns ``(ordinary, steenrod)`` barcodes.
   *
   * @param[in] n_jobs  Number of OpenMP threads for the parallelised stages.
   *                    ``-1`` (default) uses all available cores.
   */
  std::pair<Bars_by_dim, Bars_by_dim> compute(int n_jobs = -1) {
    using steenrod_persistence::Index;
    using steenrod_persistence::Filtration_by_dim;
    using steenrod_persistence::Dim_filtration;
    using steenrod_persistence::Simplex;

    // 1. Walk the filtration and bucket each simplex by dimension.
    Filtration_by_dim fbd;
    std::vector<double> filt_values;

    Index abs_idx = 0;
    for (auto sh : stptr_->filtration_simplex_range()) {
      const int dim = static_cast<int>(stptr_->dimension(sh));
      while (static_cast<int>(fbd.size()) <= dim) fbd.emplace_back();

      Simplex verts;
      for (auto v : stptr_->simplex_vertex_range(sh)) verts.push_back(static_cast<Index>(v));
      // Gudhi returns vertices in decreasing order; we need them sorted ascending.
      std::sort(verts.begin(), verts.end());

      fbd[dim].idxs.push_back(abs_idx);
      fbd[dim].tups.push_back(std::move(verts));
      filt_values.push_back(static_cast<double>(stptr_->filtration(sh)));
      ++abs_idx;
    }

    // 2. Run the steenroder pipeline (GF(2) implicit).
    auto result = steenrod_persistence::barcodes(k_, fbd, &filt_values, n_jobs, /*maxdim=*/-1);

    // 3. Convert cohomology {death_idx, birth_idx} bars to gudhi
    //    (birth_filt_value, death_filt_value) pairs per dimension.
    Bars_by_dim ordinary = to_filt_pairs(result.ordinary, filt_values);
    Bars_by_dim steenrod = to_filt_pairs(result.steenrod, filt_values);
    return {std::move(ordinary), std::move(steenrod)};
  }

 private:
  static Bars_by_dim to_filt_pairs(const steenrod_persistence::BarcodeByDim& bars_by_dim,
                                   const std::vector<double>& fv) {
    const double inf = std::numeric_limits<double>::infinity();
    Bars_by_dim out(bars_by_dim.size());
    for (std::size_t d = 0; d < bars_by_dim.size(); ++d) {
      out[d].reserve(bars_by_dim[d].size());
      for (const auto& bar : bars_by_dim[d]) {
        const auto death_idx = bar[0];
        const auto birth_idx = bar[1];
        const double birth_f = fv[birth_idx];
        const double death_f = (death_idx < 0) ? inf : fv[death_idx];
        out[d].emplace_back(birth_f, death_f);
      }
    }
    return out;
  }

  Simplex_tree_t* stptr_;
  int k_;
};

}  // namespace Gudhi

#endif  // INCLUDE_STEENROD_BARCODE_INTERFACE_H_
