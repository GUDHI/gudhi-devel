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

#include <cstddef>
#include <initializer_list>
#include <type_traits>
#include <numeric>  //std::iota
#include <utility>  //std::move
#include <cstdint>  //std::int32_t
#include <vector>
#include <ostream>
// #include <sstream>  //std::stringstream, to remove when to_str gets removed

#ifdef GUDHI_USE_TBB
#include <oneapi/tbb/enumerable_thread_specific.h>
#include <oneapi/tbb/parallel_for.h>
#endif

#include <gudhi/Debug_utils.h>
#include <gudhi/Multi_filtration/multi_filtration_utils.h>
#include <gudhi/Multi_parameter_filtered_complex.h>
#include <gudhi/Multi_persistence/Line.h>
#include <gudhi/Thread_safe_slicer.h>
#include <gudhi/Projective_cover_kernel.h>
#include <gudhi/persistence_interval.h>
#include <gudhi/slicer_helpers.h>

namespace Gudhi {
namespace multi_persistence {

template <class MultiFiltrationValue, class PersistenceAlgorithm>
class Slicer
{
 public:
  using Persistence = PersistenceAlgorithm;
  using Filtration_value = MultiFiltrationValue;
  using T = typename Filtration_value::value_type;
  using Complex = Multi_parameter_filtered_complex<Filtration_value>;
  using Index = typename Complex::Index;
  using Dimension = typename Complex::Dimension;
  template <typename Value = T>
  using Bar = Gudhi::persistence_matrix::Persistence_interval<Dimension, Value>;
  template <typename Value = T>
  using Barcode = std::vector<Bar<Value>>;
  // TODO: replace by std::vector<std::array<Value,2> > to avoid double push_back for multi dim version?
  template <typename Value = T>
  using Flat_barcode = std::vector<Value>;
  template <typename Value = T>
  using Multi_dimensional_barcode = std::vector<Barcode<Value>>;
  template <typename Value = T>
  using Multi_dimensional_flat_barcode = std::vector<Flat_barcode<Value>>;
  using Cycle = std::vector<Index>;
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
        slice_(complex_.get_number_of_cycle_generators()),
        generatorOrder_(complex_.get_number_of_cycle_generators()),
        persistence_()
  {}

  Slicer(const Slicer& other)
      : complex_(other.complex_), slice_(other.slice_), generatorOrder_(other.generatorOrder_), persistence_()
  {}

  Slicer(Slicer&& other)
      : complex_(std::move(other.complex_)),
        slice_(std::move(other.slice_)),
        generatorOrder_(std::move(other.generatorOrder_)),
        persistence_()
  {}

  // TODO: swap + assign operator?

  // ACCESS

  Thread_safe weak_copy() const { return Thread_safe(*this); }

  Index get_number_of_cycle_generators() const { return complex_.get_number_of_cycle_generators(); }

  Index get_number_of_parameters() const { return complex_.get_number_of_parameters(); }

  // // only used for scc io for now
  // const Complex& get_chain_complex() const { return complex_; }

  // initialized with `initialize_persistence_computation`
  // if ignoreInf was true, indices at inf are not contained in the vector
  const std::vector<Index>& get_current_order() const { return generatorOrder_; }

  const std::vector<T>& get_slice() const { return slice_; }

  const Persistence& get_persistence_algorithm() const { return persistence_; }

  std::pair<Filtration_value, Filtration_value> get_bounding_box() const
  {
    Filtration_value a = Filtration_value::inf(get_number_of_parameters());
    Filtration_value b = -a;
    for (const Filtration_value& fil : complex_.get_filtration_values()) {
      if (fil.num_generators() > 1) {
        a.pull_to_greatest_common_lower_bound(factorize_below(fil));
        b.push_to_least_common_upper_bound(factorize_above(fil));
      } else {
        a.pull_to_greatest_common_lower_bound(fil);
        b.push_to_least_common_upper_bound(fil);
      }
    }
    return std::make_pair(std::move(a), std::move(b));
  }

  const typename Complex::Filtration_value_container& get_filtration_values() const
  {
    return complex_.get_filtration_values();
  }

  typename Complex::Filtration_value_container& get_filtration_values() { return complex_.get_filtration_values(); }

  const Filtration_value& get_filtration_value(Index i) const { return complex_.get_filtration_values()[i]; }

  Filtration_value& get_filtration_value(Index i) { return complex_.get_filtration_values()[i]; }

  const std::vector<Dimension>& get_dimensions() const { return complex_.get_dimensions(); }

  Dimension get_dimension(Index i) const { return complex_.get_dimensions()[i]; }

  Dimension get_max_dimension() const { return complex_.get_max_dimension(); }

  const typename Complex::Boundary_container& get_boundaries() const { return complex_.get_boundaries(); }

  const typename Complex::Boundary& get_boundary(Index i) const { return complex_.get_boundaries()[i]; }

  // // TODO: only used to print info in python, so put in some interface instead
  // std::string to_str() const
  // {
  //   std::stringstream stream;
  //   stream << *this;
  //   return stream.str();
  // }

  // MODIFIERS

  template <class Array = std::initializer_list<T>>
  void set_slice(const Array& slice)
  {
    slice_ = slice;
  }

  template <class Line>
  void push_to(const Line& line)
  {
    _push_to(complex_, line);
  }

  // Warning: initialize_persistence_computation needs to be recalled if barcode was and is still needed
  // warning: shuffles order if not ordered by dimension
  void prune_above_dimension(int maxDim)
  {
    int idx = complex_.prune_above_dimension(maxDim);
    generatorOrder_.resize(idx);
    generatorOrder_.shrink_to_fit();
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
    _initialize_persistence_computation(complex_, ignoreInf);
  }

  // if ignoreInf was used when initializing the persistence computation, any update of slice has to keep at inf
  // the boundaries which were before, otherwise the behaviour is undefined (will throw with high probability)
  // preferable with complex ordered by dim
  void vineyard_update()
  {
    const bool is_ordered_by_dim = complex_.is_ordered_by_dimension();
    // speed up when ordered by dim, to avoid unnecessary swaps
    auto dim_condition = [&](int curr) {
      if (is_ordered_by_dim) {
        return persistence_.get_dimension(curr) == persistence_.get_dimension(curr - 1);
      }
      return true;
    };
    for (Index i = 1; i < generatorOrder_.size(); i++) {
      int curr = i;
      while (curr > 0 && dim_condition(curr) && slice_[generatorOrder_[curr]] < slice_[generatorOrder_[curr - 1]]) {
        persistence_.vine_swap(curr - 1);
        std::swap(generatorOrder_[curr - 1], generatorOrder_[curr]);
        --curr;
      }
    }
  }

  template <bool byDim = true, typename Value = T>
  std::conditional_t<byDim, Multi_dimensional_barcode<Value>, Barcode<Value>> get_barcode(int maxDim = -1)
  {
    if (maxDim < 0) maxDim = get_max_dimension();
    if constexpr (byDim) {
      return _get_barcode_by_dim<Value>(maxDim);
    } else {
      return _get_barcode<Value>(maxDim);
    }
  }

  template <bool byDim = false, typename Value = T>
  std::conditional_t<byDim, Multi_dimensional_flat_barcode<Value>, Flat_barcode<Value>> get_flat_barcode(
      int maxDim = -1)
  {
    if (maxDim < 0) maxDim = get_max_dimension();
    if constexpr (byDim) {
      return _get_flat_barcode_by_dim<Value>(maxDim);
    } else {
      return _get_flat_barcode<Value>(maxDim);
    }
  }

  std::vector<Multi_dimensional_barcode<T>> persistence_on_lines(const std::vector<std::vector<T>>& basePoints,
                                                                 [[maybe_unused]] bool ignoreInf = true)
  {
    // TODO: Thread_safe has to use his own version of weak_copy(), so I had to decompose everything to factorize
    // but it is quite ugly. Does someone has a more elegant solution?
    if constexpr (Persistence::is_vine) {
      return _batch_persistence_on_lines_with_vine(
          complex_, [](const std::vector<T>& bp) { return Line<T>(bp); }, basePoints);
    } else {
#ifdef GUDHI_USE_TBB
      return _batch_persistence_on_lines(
          weak_copy(), [](const std::vector<T>& bp) { return Line<T>(bp); }, basePoints, ignoreInf);
#else
      return _batch_persistence_on_lines(
          complex_, [](const std::vector<T>& bp) { return Line<T>(bp); }, basePoints, ignoreInf);
#endif
    }
  }

  std::vector<Multi_dimensional_barcode<T>> persistence_on_lines(
      const std::vector<std::pair<std::vector<T>, std::vector<T>>>& basePointsWithDirections,
      [[maybe_unused]] bool ignoreInf = true)
  {
    // TODO: Thread_safe has to use his own version of weak_copy(), so I had to decompose everything to factorize
    // but it is quite ugly. Does someone has a more elegant solution?
    if constexpr (Persistence::is_vine) {
      return _batch_persistence_on_lines_with_vine(
          complex_,
          [](const std::pair<std::vector<T>, std::vector<T>>& bpwd) { return Line<T>(bpwd.first, bpwd.second); },
          basePointsWithDirections);
    } else {
#ifdef GUDHI_USE_TBB
      return _batch_persistence_on_lines(
          weak_copy(),
          [](const std::pair<std::vector<T>, std::vector<T>>& bpwd) { return Line<T>(bpwd.first, bpwd.second); },
          basePointsWithDirections,
          ignoreInf);
#else
      return _batch_persistence_on_lines(
          complex_,
          [](const std::pair<std::vector<T>, std::vector<T>>& bpwd) { return Line<T>(bpwd.first, bpwd.second); },
          basePointsWithDirections,
          ignoreInf);
#endif
    }
  }

  std::vector<std::vector<Cycle>> get_representative_cycles(bool update = true)
  {
    return _get_representative_cycles(complex_, update);
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

  friend Slicer build_slicer_from_projective_cover_kernel(const Slicer& slicer, Dimension dim)
  {
    Projective_cover_kernel<Filtration_value> kernel(slicer.complex_, dim);
    return Slicer(kernel.create_complex());
  }

  friend void write_slicer_to_scc_file(const std::string& outFilePath,
                                       const Slicer& slicer,
                                       int numberOfParameters = -1,
                                       int degree = -1,
                                       bool rivetCompatible = false,
                                       bool ignoreLastGenerators = false,
                                       bool stripComments = false,
                                       bool reverse = false)
  {
    const Complex& cpx =
        slicer.complex_.is_ordered_by_dimension() ? slicer.complex_ : build_permuted_complex(slicer.complex_).first;
    write_complex_to_scc_file<typename Slicer::Filtration_value>(
        outFilePath, cpx, numberOfParameters, degree, rivetCompatible, ignoreLastGenerators, stripComments, reverse);
  };

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
  friend Thread_safe;

  // For ThreadSafe version
  Slicer(const std::vector<T>& slice, const std::vector<Index>& generatorOrder, const Persistence& persistence)
      : complex_(), slice_(slice), generatorOrder_(generatorOrder), persistence_(persistence, generatorOrder_)
  {}

  template <class Line>
  void _push_to(const Complex& complex, const Line& line)
  {
    const auto& filtrationValues = complex.get_filtration_values();
    for (Index i = 0u; i < filtrationValues.size(); i++) {
      slice_[i] = line.template compute_forward_intersection<T>(filtrationValues[i]);
    }
  }

  void _initialize_persistence_computation(const Complex& complex, const bool ignoreInf = true)
  {
    _initialize_order(complex, ignoreInf);
    persistence_.reinitialize(complex, generatorOrder_);
  }

  std::vector<std::vector<Cycle>> _get_representative_cycles(const Complex& complex, bool update = true)
  {
    static_assert(Persistence::has_rep_cycles, "Representative cycles not enabled.");

    const auto& dimensions = complex.get_dimensions();
    auto cycleKeys = persistence_.get_representative_cycles(update);
    auto numCycles = cycleKeys.size();
    std::vector<std::vector<Cycle>> out(complex.get_max_dimension() + 1);
    for (auto& cyclesDim : out) cyclesDim.reserve(numCycles);
    for (const auto& cycle : cycleKeys) {
      GUDHI_CHECK(!cycle.empty(), "A cycle should not be empty...");
      // assumes cycle to be never empty & all faces have same dimension
      out[dimensions[cycle[0]]].push_back(cycle);
    }
    return out;
  }

  template <typename F, typename F_arg>
  std::vector<Multi_dimensional_barcode<T>> _batch_persistence_on_lines_with_vine(const Complex& complex,
                                                                                  F&& get_line,
                                                                                  const std::vector<F_arg>& basePoints)
  {
    if (basePoints.size() == 0) return {};

    std::vector<Multi_dimensional_barcode<T>> out(basePoints.size());

    _push_to(complex, get_line(basePoints[0]));
    _initialize_persistence_computation(complex, false);
    out[0] = _get_barcode_by_dim<T>(complex.get_max_dimension());
    for (auto i = 1u; i < basePoints.size(); ++i) {
      _push_to(complex, get_line(basePoints[i]));
      vineyard_update();
      out[i] = _get_barcode_by_dim<T>(complex.get_max_dimension());
    }

    return out;
  }

#ifdef GUDHI_USE_TBB
  template <typename F, typename F_arg>
  std::vector<Multi_dimensional_barcode<T>> _batch_persistence_on_lines(const Thread_safe& localTemplate,
                                                                        F&& get_line,
                                                                        const std::vector<F_arg>& basePoints,
                                                                        const bool ignoreInf)
  {
    if (basePoints.size() == 0) return {};

    std::vector<Multi_dimensional_barcode<T>> out(basePoints.size());

    tbb::enumerable_thread_specific<Thread_safe> threadLocals(localTemplate);
    tbb::parallel_for(static_cast<std::size_t>(0), basePoints.size(), [&](const Index& i) {
      Thread_safe& s = threadLocals.local();
      s.push_to(get_line(basePoints[i]));
      s.initialize_persistence_computation(ignoreInf);
      out[i] = s.get_barcode();
    });

    return out;
  }
#else
  template <typename F, typename F_arg>
  std::vector<Multi_dimensional_barcode<T>> _batch_persistence_on_lines(const Complex& complex,
                                                                        F&& get_line,
                                                                        const std::vector<F_arg>& basePoints,
                                                                        const bool ignoreInf)
  {
    if (basePoints.size() == 0) return {};

    std::vector<Multi_dimensional_barcode<T>> out(basePoints.size());

    for (auto i = 0u; i < basePoints.size(); ++i) {
      _push_to(complex, get_line(basePoints[i]));
      _initialize_persistence_computation(complex, ignoreInf);
      out[i] = _get_barcode_by_dim<T>(complex.get_max_dimension());
    }

    return out;
  }
#endif

 private:
  Complex complex_;
  std::vector<T> slice_;  // filtration of the current slice
  std::vector<Index> generatorOrder_;
  Persistence persistence_;

  void _initialize_order(const Complex& complex, const bool ignoreInf = true)
  {
    const auto& dimensions = complex.get_dimensions();
    generatorOrder_.resize(complex.get_number_of_cycle_generators());
    std::iota(generatorOrder_.begin(), generatorOrder_.end(), 0);
    std::sort(generatorOrder_.begin(), generatorOrder_.end(), [&](Index i, Index j) {
      if (ignoreInf) {
        if (slice_[i] != Filtration_value::T_inf && slice_[j] == Filtration_value::T_inf) return true;
        // all elements at inf are considered equal
        if (slice_[i] == Filtration_value::T_inf) return false;
      }
      if (dimensions[i] > dimensions[j]) return false;
      if (dimensions[i] < dimensions[j]) return true;
      // if filtration values are equal, we don't care about order, so considered the same object
      return slice_[i] < slice_[j];
    });
    if (ignoreInf) {
      Index end = generatorOrder_.size();
      while (end > 0 && slice_[generatorOrder_[end - 1]] == Filtration_value::T_inf) --end;
      generatorOrder_.resize(end);
    }
  }

  template <class Interval, typename Value>
  void _retrieve_interval(const Interval& bar, Dimension& dim, Value& birth, Value& death)
  {
    const Value inf = Gudhi::multi_filtration::MF_T_inf<Value>;
    dim = bar.dim;
    birth = slice_[bar.birth];
    death = inf;
    if (bar.death != Persistence::nullDeath) death = slice_[bar.death];
    if (!(birth <= death)) {
      birth = inf;
      death = inf;
    }
  }

  template <typename Value>
  Barcode<Value> _get_barcode(int maxDim)
  {
    auto barcodeIndices = persistence_.get_barcode();
    Barcode<Value> out(barcodeIndices.size());
    Index i = 0;
    for (const auto& bar : barcodeIndices) {
      if (bar.dim <= maxDim) {
        _retrieve_interval(bar, out[i].dim, out[i].birth, out[i].death);
        ++i;
      }
    }
    out.resize(i);
    return out;
  }

  template <typename Value>
  Multi_dimensional_barcode<Value> _get_barcode_by_dim(int maxDim)
  {
    // TODO: This doesn't allow for negative dimensions
    // Hannah: not sure what this comment means ?
    Multi_dimensional_barcode<Value> out(maxDim + 1);
    Value birth, death;
    Dimension dim;
    for (const auto& bar : persistence_.get_barcode()) {
      if (bar.dim <= maxDim) {
        _retrieve_interval(bar, dim, birth, death);
        out[dim].emplace_back(birth, death, dim);
      }
    }
    return out;
  }

  template <typename Value>
  Flat_barcode<Value> _get_flat_barcode(int maxDim)
  {
    auto barcodeIndices = persistence_.get_barcode();
    Flat_barcode<Value> out(barcodeIndices.size() * 2);
    Index i = 0;
    Dimension dim;  // dummy
    for (const auto& bar : barcodeIndices) {
      if (bar.dim <= maxDim) {
        _retrieve_interval(bar, dim, out[i], out[i + 1]);
        i += 2;
      }
    }
    out.resize(i);
    return out;
  }

  template <typename Value>
  Multi_dimensional_flat_barcode<Value> _get_flat_barcode_by_dim(int maxDim)
  {
    Multi_dimensional_flat_barcode<Value> out(maxDim + 1);
    Value birth, death;
    Dimension dim;
    for (const auto& bar : persistence_.get_barcode()) {
      if (bar.dim <= maxDim) {
        _retrieve_interval(bar, dim, birth, death);
        out[dim].push_back(birth);
        out[dim].push_back(death);
      }
    }
    return out;
  }
};

}  // namespace multi_persistence
}  // namespace Gudhi

#endif  // MP_SLICER_H_INCLUDED
