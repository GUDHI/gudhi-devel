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

#include <algorithm>
#include <array>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

#include <nanobind/nanobind.h>

#include <gudhi/Steenrod_barcode.h>

#include <python_interfaces/numpy_utils.h>

namespace Gudhi {

/** \brief Python-interface wrapper that exposes Steenrod persistence barcodes
 *         to the nanobind SimplexTree bindings.
 *
 * Usage mirrors ``Persistent_cohomology_interface``:
 *   1. Construct with a pointer to a ``Simplex_tree_interface`` and the Sq^k exponent.
 *   2. Call ``_compute()``; receive a nanobind tuple
 *      ``(ordinary, steenrod)`` where each element is itself a 2-tuple
 *      ``(finites_per_dim, infinites_per_dim)``:
 *
 *        * ``finites_per_dim``  — ``list`` of numpy arrays of shape ``(n_bars, 2)``,
 *          rows are ``(death_value, birth_value)`` with ``death < birth``,
 *          encoding the relative-cohomology bar ``[a_p, a_{q+1})`` of
 *          Lupo, Medina-Mardones, Tauzin (2022) §2.4.
 *        * ``infinites_per_dim`` — ``list`` of numpy arrays of shape ``(n_bars,)``,
 *          birth values of essential bars; the implicit lower endpoint is
 *          ``-inf`` (a_0).
 *
 * Coefficients are GF(2); no field parameter is exposed.
 */
template <class Simplex_tree_t>
class Steenrod_barcode_interface {
 public:
  using Filtration_value = typename Simplex_tree_t::Filtration_value;

  Steenrod_barcode_interface(Simplex_tree_t* stptr, int k) : stptr_(stptr), k_(k) {
    if (k < 0) throw std::invalid_argument("Sq^k exponent k must be non-negative");
  }

  /** Returns ``(ordinary, steenrod)`` in the (finite, infinite) per-dim
   * format in the **relative-cohomology convention**.
   *
   * @param[in] n_jobs  Number of OpenMP threads for the parallelised stages.
   *                    ``-1`` (default) uses all available cores.
   */
  nanobind::tuple compute(int n_jobs = -1) {
    std::vector<double> filt_values;
    auto result = run_pipeline(n_jobs, filt_values);

    auto ord_split = split_to_vectors(result.ordinary, filt_values);
    auto st_split  = split_to_vectors(result.steenrod, filt_values);

    nanobind::gil_scoped_acquire acquire;
    return nanobind::make_tuple(wrap_split(std::move(ord_split)),
                                wrap_split(std::move(st_split)));
  }

  /** Returns ``(ordinary, steenrod, problematic_dims)`` in the
   * **absolute-cohomology convention**.
   *
   * Conversion from the relative output is purely structural: a finite bar
   * at relative dim ``d`` represents an absolute class in degree ``d - 1``,
   * so finite bars shift down one dimension; essential bars represent
   * absolute classes in the same degree and stay where they are.  The
   * stored numerical values do not change.
   *
   * The conversion is well-defined only when the relative ordinary barcode
   * has no essential bars at degrees ``[1, max_dim]`` (the duality
   * condition; see LMT 2022 §3.3).  ``problematic_dims`` lists the
   * degrees that violate the condition; the Python wrapper turns it into
   * a ``UserWarning``.
   *
   * @param[in] n_jobs   OpenMP threads.  ``-1`` (default) uses all cores.
   * @param[in] max_dim  Largest absolute dimension to return.  Negative or
   *                     out-of-range values keep every dimension of the
   *                     complex.  Higher-dimensional simplices remain in
   *                     the reduction so classes at ``H^max_dim`` can be
   *                     killed correctly; only the returned arrays are
   *                     truncated.
   */
  nanobind::tuple compute_absolute(int n_jobs = -1, int max_dim = -1) {
    std::vector<double> filt_values;
    auto result = run_pipeline(n_jobs, filt_values);

    auto ord_split = split_to_vectors(result.ordinary, filt_values);
    auto st_split  = split_to_vectors(result.steenrod, filt_values);

    const int n_dims = static_cast<int>(ord_split.infinite.size());
    const int eff_max = (max_dim < 0 || max_dim >= n_dims) ? n_dims - 1 : max_dim;

    // Duality-condition check: any essential ordinary bar at relative
    // dim d in [1, eff_max] is a violation.  d = 0 is excluded — H^0
    // essentials (connected-component bars) cannot affect any absolute
    // Sq^k bar with k >= 1.
    std::vector<int> problematic_dims;
    for (int d = 1; d <= eff_max; ++d) {
      if (!ord_split.infinite[d].empty()) problematic_dims.push_back(d);
    }

    auto ord_abs = shift_to_absolute(std::move(ord_split), eff_max);
    auto st_abs  = shift_to_absolute(std::move(st_split),  eff_max);

    nanobind::gil_scoped_acquire acquire;
    nanobind::list problematic_list;
    for (int d : problematic_dims) problematic_list.append(d);
    return nanobind::make_tuple(wrap_split(std::move(ord_abs)),
                                wrap_split(std::move(st_abs)),
                                std::move(problematic_list));
  }

 private:
  /** Walk the simplex tree, build a ``Filtration_by_dim``, and run the
   * full steenroder reduction.  Pure GIL-free work: this function does
   * not touch any Python state.
   */
  steenrod_persistence::Barcodes_result run_pipeline(int n_jobs,
                                                     std::vector<double>& filt_values) {
    using steenrod_persistence::Index;
    using steenrod_persistence::Filtration_by_dim;
    using steenrod_persistence::Simplex;

    Filtration_by_dim fbd;
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

    return steenrod_persistence::barcodes(k_, fbd, &filt_values, n_jobs, /*maxdim=*/-1);
  }


  /** Per-dim partitioning of the raw {death_idx, birth_idx} bars into
   * finite ``(a_p, a_{q+1})`` 2-tuples and infinite single ``a_{q+1}``
   * values.  No Python state touched.
   */
  struct Split {
    std::vector<std::vector<std::array<double, 2>>> finite;    // per dim
    std::vector<std::vector<double>>                infinite;  // per dim
  };

  static Split split_to_vectors(const steenrod_persistence::BarcodeByDim& bars_by_dim,
                                const std::vector<double>& fv) {
    Split out;
    const std::size_t n_dims = bars_by_dim.size();
    out.finite.resize(n_dims);
    out.infinite.resize(n_dims);
    for (std::size_t d = 0; d < n_dims; ++d) {
      out.finite[d].reserve(bars_by_dim[d].size());
      for (const auto& bar : bars_by_dim[d]) {
        const auto death_idx = bar[0];
        const auto birth_idx = bar[1];
        if (death_idx < 0) {
          // Essential bar: half-open interval [-inf, a_{q+1}) in relative
          // cohomology.  Expose only the upper endpoint a_{q+1}.
          out.infinite[d].push_back(fv[birth_idx]);
        } else {
          // Finite bar: ``(death_value, birth_value)``, encoding the
          // relative-cohomology interval [a_p, a_{q+1}); the reduction
          // guarantees fv[death_idx] < fv[birth_idx].
          out.finite[d].push_back({fv[death_idx], fv[birth_idx]});
        }
      }
    }
    return out;
  }

  /** Apply the relative→absolute dim shift to a ``Split`` and truncate to
   * absolute dimensions ``[0, eff_max]``.
   *
   * Mapping:
   *   * finite_abs[D] = finite_rel[D + 1] for D in [0, eff_max]  (shift)
   *   * infinite_abs[D] = infinite_rel[D] for D in [0, eff_max]  (no shift)
   *
   * Numerical values are preserved.  Pure GIL-free work.
   */
  static Split shift_to_absolute(Split&& rel, int eff_max) {
    const int n_dims = static_cast<int>(rel.infinite.size());
    Split abs;
    const int out_size = eff_max + 1;
    abs.finite.resize(out_size);
    abs.infinite.resize(out_size);
    for (int D = 0; D < out_size; ++D) {
      abs.infinite[D] = std::move(rel.infinite[D]);
      if (D + 1 < n_dims) {
        abs.finite[D] = std::move(rel.finite[D + 1]);
      }
    }
    return abs;
  }

  /** Wrap each per-dim buffer as a zero-copy numpy array via
   * ``_wrap_as_numpy_array``.  Must be called with the GIL held.
   */
  static nanobind::tuple wrap_split(Split&& s) {
    nanobind::list finite_per_dim;
    nanobind::list infinite_per_dim;
    for (std::size_t d = 0; d < s.finite.size(); ++d) {
      // Capture sizes before std::move so the shape argument cannot be
      // sequenced after the move.
      const auto inf_size = s.infinite[d].size();
      finite_per_dim.append(_wrap_as_numpy_array(std::move(s.finite[d])));
      infinite_per_dim.append(_wrap_as_numpy_array(std::move(s.infinite[d]), inf_size));
    }
    return nanobind::make_tuple(finite_per_dim, infinite_per_dim);
  }

  Simplex_tree_t* stptr_;
  int k_;
};

}  // namespace Gudhi

#endif  // INCLUDE_STEENROD_BARCODE_INTERFACE_H_
