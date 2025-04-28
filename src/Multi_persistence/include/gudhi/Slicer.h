/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       David Loiseaux
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - 2025/04 Hannah Schreiber: Reorganization.
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file Slicer.h
 * @author David Loiseaux
 * @brief Contains the @ref Gudhi::multi_persistence::Slicer class.
 */

#ifndef MP_SLICER_H_INCLUDED
#define MP_SLICER_H_INCLUDED

#include <initializer_list>
#include <numeric>  //std::iota
#include <utility>  //std::move
#include <cstdint>  //std::int32_t
#include <vector>
#include <ostream>
#include <sstream>  //std::stringstream, to remove when to_str gets removed

#ifdef GUDHI_USE_TBB
#include <oneapi/tbb/enumerable_thread_specific.h>
#include <oneapi/tbb/parallel_for.h>
#endif

#include <gudhi/Debug_utils.h>
#include <gudhi/Multi_filtration/multi_filtration_utils.h>
#include <gudhi/Multi_parameter_filtered_complex.h>
#include <gudhi/Multi_persistence/Line.h>
#include <gudhi/Thread_safe_slicer.h>

namespace Gudhi {
namespace multi_persistence {

template <class Multi_filtration_value, class Persistence_algorithm>
class Slicer
{
 public:
  using Persistence = Persistence_algorithm;
  using Filtration_value = Multi_filtration_value;
  using T = typename Filtration_value::value_type;
  using Complex = Multi_parameter_filtered_complex<Filtration_value>;
  using Index = typename Complex::Index;
  using Dimension = typename Complex::Dimension;
  // TODO: again another barcode type...
  // Homogenize with other barcode types and use the flat barcode for Python instead?
  template <typename Value = T>
  using Barcode = std::vector<std::vector<std::pair<Value, Value>>>;
  template <typename Value = T>
  using Flat_barcode = std::vector<Value>;
  using Boundary = std::vector<Index>;
  using Cycle = std::vector<Boundary>;
  using Thread_safe = Thread_safe_slicer<Slicer>;

  // CONSTRUCTORS

  Slicer() {}

  Slicer(const Complex& complex)
      : complex_(complex),
        slice_(complex.get_number_of_cycle_generators()),
        generatorOrder_(complex.get_number_of_cycle_generators()),
        persistence_()
  {}

  Slicer(Complex&& complex)
      : complex_(std::move(complex)),
        slice_(complex.get_number_of_cycle_generators()),
        generatorOrder_(complex.get_number_of_cycle_generators()),
        persistence_()
  {}

  // ACCESS

  Thread_safe weak_copy() const { return Thread_safe(*this); }

  Index get_number_of_cycle_generators() const { return complex_.get_number_of_cycle_generators(); }

  Index get_number_of_parameters() const { return complex_.get_number_of_parameters(); }

  // only used for scc io TODO: see if still necessary
  const Complex& get_chain_complex() const { return complex_; }

  const std::vector<Index>& get_current_order() const { return generatorOrder_; }

  const std::vector<T>& get_slice() const { return slice_; }

  const Persistence& get_persistence_algorithm() const { return persistence_; }

  std::pair<Filtration_value, Filtration_value> get_bounding_box() const
  {
    Filtration_value a = Filtration_value::inf(get_number_of_parameters());
    Filtration_value b = -a;
    for (const Filtration_value& filtration_value : complex_.get_filtration_values()) {
      if (filtration_value.num_generators() > 1) {
        a.pull_to_greatest_common_lower_bound(factorize_below(filtration_value));
        b.push_to_least_common_upper_bound(factorize_above(filtration_value));
      } else {
        a.pull_to_greatest_common_lower_bound(filtration_value);
        b.push_to_least_common_upper_bound(filtration_value);
      }
    }
    return std::make_pair(std::move(a), std::move(b));
  }

  const typename Complex::Filtration_value_container& get_filtration_values() const
  {
    return complex_.get_filtration_values();
  }

  typename Complex::Filtration_value_container& get_filtration_values() { return complex_.get_filtration_values(); }

  const std::vector<Dimension>& get_dimensions() const { return complex_.get_dimensions(); }

  // std::vector<Dimension>& get_dimensions() { return complex_.get_dimensions(); }

  const typename Complex::Boundary_container& get_boundaries() const { return complex_.get_boundaries(); }

  // Complex::Boundary_container& get_boundaries() { return complex_.get_boundaries(); }

  Dimension get_dimension(Index i) const { return complex_.get_dimension(i); }

  // TODO: only used to print info in python, so put in some interface instead
  std::string to_str() const
  {
    std::stringstream stream;
    stream << *this;
    return stream.str();
  }

  // MODIFIERS

  template <class Array = std::initializer_list<T>>
  void set_slice(const Array& slice)
  {
    // TODO: resize with inf values instead of gudhi_check ?
    GUDHI_CHECK(slice.size() == get_number_of_cycle_generators(), "The given slice does not have the right size.");
    slice_ = slice;
  }

  template <class Line>
  void push_to(const Line& line)
  {
    _push_to(this, line);
  }

  void prune_above_dimension(int maxDim)
  {
    int idx = complex_.prune_above_dimension(maxDim);
    generatorOrder_.resize(idx);
    slice_.resize(idx);
  }

  void coarsen_on_grid(const std::vector<std::vector<T>>& grid, bool coordinate = true)
  {
    complex_.coarsen_on_grid(grid, coordinate);
  }

  // PERSISTENCE

  bool persistence_computation_is_initialized() const { return persistence_.is_initialized(); }

  void initialize_persistence_computation(const bool ignoreInf = true)
  {
    _initialize_persistence_computation(this, ignoreInf);
  }

  void vineyard_update()
  {
    for (std::size_t i = 0; i < get_number_of_cycle_generators(); i++) {
      auto j = i;
      while (j > 0 && persistence_.get_dimension(j) == persistence_.get_dimension(j - 1) &&
             slice_[generatorOrder_[j]] < slice_[generatorOrder_[j - 1]]) {
        persistence_.vine_swap(j - 1);
        std::swap(generatorOrder_[j - 1], generatorOrder_[j]);
        j--;
      }
    }
  }

  template <typename Value = T>
  Barcode<Value> get_barcode(int maxDim = -1)
  {
    if (maxDim < 0) maxDim = complex_.get_max_dimension();
    auto barcode_indices = persistence_.get_barcode();
    // TODO : This doesn't allow for negative dimensions
    // Hannah: not sure what this comment means ?
    Barcode<Value> out(maxDim + 1);
    const Value inf = Gudhi::multi_filtration::MF_T_inf<Value>;
    const auto null_death = Persistence::nullDeath;
    for (const auto& bar : barcode_indices) {
      Value birth_filtration = slice_[bar.birth];
      Value death_filtration = inf;
      if (bar.death != null_death) death_filtration = slice_[bar.death];
      if (birth_filtration <= death_filtration)
        out[bar.dim].emplace_back({birth_filtration, death_filtration});
      else {
        out[bar.dim].emplace_back({inf, inf});
      }
    }
    return out;
  }

  template <bool withDim = false, typename Value = T>
  Flat_barcode<Value> get_flat_barcode()
  {
    auto barcode_indices = persistence_.get_barcode();
    Flat_barcode<Value> out(barcode_indices.size() * (withDim ? 3 : 2));
    const Value inf = Gudhi::multi_filtration::MF_T_inf<Value>;
    const auto null_death = Persistence::nullDeath;
    Index i = 0;
    for (const auto& bar : barcode_indices) {
      Value birth_filtration = slice_[bar.birth];
      Value death_filtration = inf;
      if (bar.death != null_death) death_filtration = slice_[bar.death];
      if (birth_filtration > death_filtration) {
        birth_filtration = inf;
        death_filtration = inf;
      }
      if constexpr (withDim) {
        out[i] = bar.dim;
        out[i + 1] = birth_filtration;
        out[i + 2] = death_filtration;
        i += 3;
      } else {
        out[i] = birth_filtration;
        out[i + 1] = death_filtration;
        i += 2;
      }
    }
    return out;
  }

  std::vector<Barcode<T>> persistence_on_lines(const std::vector<std::vector<T>>& basePoints, bool ignoreInf)
  {
    return _batch_persistence_on_lines([](const std::vector<T>& bp) { return Line<T>(bp); }, basePoints, ignoreInf);
  }

  std::vector<Barcode<T>> persistence_on_lines(
      const std::vector<std::pair<std::vector<T>, std::vector<T>>>& basePointsWithDirections,
      bool ignoreInf)
  {
    return _batch_persistence_on_lines(
        [](const std::pair<std::vector<T>, std::vector<T>>& bpwd) { return Line<T>(bpwd.first, bpwd.second); },
        basePointsWithDirections,
        ignoreInf);
  }

  std::vector<std::vector<Cycle>> get_representative_cycles(bool update = true, bool detailed = false)
  {
    return _get_representative_cycles(this, update, detailed);
  }

  // FRIENDS

  friend Slicer build_permuted_slicer(const Slicer& slicer, const std::vector<Index>& permutation)
  {
    GUDHI_CHECK(permutation.size() > slicer.get_number_of_cycle_generators(),
                "Too many elements in permutation vector.");
    return Slicer(build_permuted_complex(slicer.complex_, permutation));
  }

  friend std::pair<Slicer, std::vector<Index>> build_permuted_slicer(const Slicer& slicer)
  {
    auto [complex, permutation] = build_permuted_complex(slicer.complex_);
    return std::make_pair(Slicer(std::move(complex)), std::move(permutation));
  }

  friend auto build_slicer_coarsen_on_grid(const Slicer& slicer, const std::vector<std::vector<T>> grid)
  {
    using return_filtration_value = decltype(std::declval<Filtration_value>().template as_type<std::int32_t>());
    return Slicer<return_filtration_value, Persistence>(build_complex_coarsen_on_grid(slicer.complex_, grid));
  }

  friend std::ostream& operator<<(std::ostream& stream, const Slicer& slicer)
  {
    stream << "-------------------- Slicer \n";

    stream << "--- Filtered complex \n";
    stream << slicer.complex_;

    stream << "--- Order \n";
    stream << "{";
    for (const auto& idx : slicer.generatorOrder_) stream << idx << ", ";
    stream << "}" << std::endl;

    stream << "--- Current slice filtration\n";
    stream << "{";
    for (const auto& val : slicer.slice_) stream << val << ", ";
    stream << "\b" << "\b";
    stream << "}" << std::endl;

    stream << "--- PersBackend \n";
    stream << slicer.persistence_;

    return stream;
  }

 protected:
  // For ThreadSafe version
  Slicer(const std::vector<T>& slice, const std::vector<Index>& generatorOrder, const Persistence& persistence)
      : complex_(), slice_(slice), generatorOrder_(generatorOrder), persistence_(persistence)
  {
    persistence_.update_permutation_pointer(&generatorOrder_);
  }

  template <class Line>
  void _push_to(const Slicer* slicer, const Line& line)
  {
    const auto& filtrationValues = slicer->complex_.get_filtration_values();
    for (std::size_t i = 0u; i < filtrationValues.size(); i++) {
      slice_[i] = line.template compute_forward_intersection<T>(filtrationValues[i]);
    }
  }

  void _initialize_persistence_computation(const Slicer* slicer, [[maybe_unused]] const bool ignoreInf = true)
  {
    _initialize_order();
    // if is_vine is true, ignoring the inf values causes an indexation problem when translating back the barcode values
    // can be corrected by storing additional data, but does not seem to be worth it
    if constexpr (!Persistence::is_vine) {
      if (ignoreInf) {
        for (Index& i : generatorOrder_) {
          if (slice_[i] == Filtration_value::T_inf) {
            i = std::numeric_limits<Index>::max();
          }
        }
      }
    }
    persistence_ = Persistence(slicer->complex_, generatorOrder_);
  }

  std::vector<std::vector<Cycle>> _get_representative_cycles(const Slicer* slicer,
                                                             bool update = true,
                                                             bool detailed = false)
  {
    static_assert(Persistence::has_rep_cycles, "Representative cycles not enabled.");

    // iterable iterable simplex key
    auto cycles_key = persistence_.get_representative_cycles(update, detailed);
    auto num_cycles = cycles_key.size();
    std::vector<std::vector<Cycle>> out(slicer->complex_.get_max_dimension() + 1);
    for (auto& cycles_of_dim : out) cycles_of_dim.reserve(num_cycles);
    for (const auto& cycle : cycles_key) {
      int cycle_dim = 0;        // for more generality, should be minimal dimension instead
      if (!cycle[0].empty()) {  // if empty, cycle has no border -> assumes dimension 0 even if it could be min dim
        cycle_dim = slicer->complex_.dimension(cycle[0][0]) + 1;  // all faces have the same dim
      }
      out[cycle_dim].push_back(cycle);
    }
    return out;
  }

 private:
  Complex complex_;
  std::vector<T> slice_;  // filtration of the current slice
  std::vector<Index> generatorOrder_;
  Persistence persistence_;

  void _initialize_order()
  {
    std::iota(generatorOrder_.begin(), generatorOrder_.end(), 0);
    std::sort(generatorOrder_.begin(), generatorOrder_.end(), [&](Index i, Index j) {
      if (complex_.get_dimension(i) > complex_.get_dimension(j)) return false;
      if (complex_.get_dimension(i) < complex_.get_dimension(j)) return true;
      return slice_[i] < slice_[j];
    });
  }

  template <typename F, typename F_arg>
  std::vector<Barcode<T>> _batch_persistence_on_lines(F&& get_line,
                                                      const std::vector<F_arg>& basePoints,
                                                      [[maybe_unused]] const bool ignoreInf = true)
  {
    if (basePoints.size() == 0) return {};

    std::vector<Barcode<T>> out(basePoints.size());

    if constexpr (Persistence::is_vine) {
      push_to(get_line(basePoints[0]));
      initialize_persistence_computation();
      out[0] = get_barcode();
      for (auto i = 1u; i < basePoints.size(); ++i) {
        push_to(get_line(basePoints[i]));
        vineyard_update();
        out[i] = get_barcode();
      }
    } else {
#ifdef GUDHI_USE_TBB
      Thread_safe local_template = weak_copy();
      tbb::enumerable_thread_specific<Thread_safe> thread_locals(local_template);
      tbb::parallel_for(static_cast<std::size_t>(0), basePoints.size(), [&](const std::size_t& i) {
        Thread_safe& s = thread_locals.local();
        s.push_to(get_line(basePoints[i]));
        s.compute_persistence(ignoreInf);
        out[i] = s.get_barcode();
      });
#else
      for (auto i = 0u; i < basePoints.size(); ++i) {
        push_to(get_line(basePoints[i]));
        initialize_persistence_computation(ignoreInf);
        out[i] = get_barcode();
      }
#endif
    }
    return out;
  }
};

}  // namespace multi_persistence
}  // namespace Gudhi

#endif  // MP_SLICER_H_INCLUDED
