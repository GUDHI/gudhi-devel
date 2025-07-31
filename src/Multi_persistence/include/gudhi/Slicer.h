/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       David Loiseaux
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - 2025/04 Hannah Schreiber: Reorganization + documentation.
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
#include <gudhi/Multi_persistence/Line.h>

namespace Gudhi {
namespace multi_persistence {

/**
 * @class Slicer Slicer.h gudhi/Slicer.h
 * @ingroup multi_persistence
 *
 * @brief Class slicing a multi-parameter persistence module. TODO: more details
 * 
 * @tparam MultiFiltrationValue Filtration value class respecting the @ref MultiFiltrationValue concept.
 * @tparam PersistenceAlgorithm Class respecting the @ref PersistenceAlgorithm concept. Used to compute persistence,
 * eventually vineyards and representative cycles.
 */
template <class MultiFiltrationValue, class PersistenceAlgorithm>
class Slicer
{
 public:
  using Persistence = PersistenceAlgorithm;                           /**< Persistence algorithm type. */
  using Filtration_value = MultiFiltrationValue;                      /**< Filtration value type. */
  using T = typename Filtration_value::value_type;                    /**< Numerical filtration value element type. */
  using Complex = Multi_parameter_filtered_complex<Filtration_value>; /**< Complex type. */
  using Index = typename Complex::Index;                              /**< Complex index type. */
  using Dimension = typename Complex::Dimension;                      /**< Dimension type. */
  template <typename Value = T>
  using Bar = Gudhi::persistence_matrix::Persistence_interval<Dimension, Value>;  /**< Bar type. */
  template <typename Value = T>
  using Barcode = std::vector<Bar<Value>>;                                 /**< Barcode type. */
  // TODO: replace by std::vector<std::array<Value,2> > to avoid double push_back for multi dim version?
  template <typename Value = T>
  using Flat_barcode = std::vector<Value>;                                 /**< Flat barcode type. */
  template <typename Value = T>
  using Multi_dimensional_barcode = std::vector<Barcode<Value>>;           /**< Barcode ordered by dimension type. */
  template <typename Value = T>
  using Multi_dimensional_flat_barcode = std::vector<Flat_barcode<Value>>; /**< Flat barcode ord. by dimension type. */
  using Cycle = std::vector<Index>;                                        /**< Cycle type. */
  using Thread_safe = Thread_safe_slicer<Slicer>;                          /**< Thread safe slicer type. */

  // CONSTRUCTORS

  /**
   * @brief Default constructor. Constructs an empty slicer.
   */
  Slicer() {}

  /**
   * @brief Constructs the slicer by copying the given complex. The current slice is not initialized to a default
   * value, it can be set with @ref set_slice or @ref push_to.
   *
   * It is recommended to use a complex which is ordered by dimension for better performance.
   */
  Slicer(const Complex& complex)
      : complex_(complex),
        slice_(complex.get_number_of_cycle_generators()),
        generatorOrder_(complex.get_number_of_cycle_generators()),
        persistence_()
  {}

  /**
   * @brief Constructs the slicer by moving the given complex. The current slice is not initialized to a default
   * value, it can be set with @ref set_slice or @ref push_to.
   *
   * It is recommended to use a complex which is ordered by dimension for better performance.
   */
  Slicer(Complex&& complex)
      : complex_(std::move(complex)),
        slice_(complex_.get_number_of_cycle_generators()),
        generatorOrder_(complex_.get_number_of_cycle_generators()),
        persistence_()
  {}

  /**
   * @brief Copy constructor.
   */
  Slicer(const Slicer& other)
      : complex_(other.complex_), slice_(other.slice_), generatorOrder_(other.generatorOrder_), persistence_()
  {}

  /**
   * @brief Move constructor.
   */
  Slicer(Slicer&& other)
      : complex_(std::move(other.complex_)),
        slice_(std::move(other.slice_)),
        generatorOrder_(std::move(other.generatorOrder_)),
        persistence_()
  {}

  // TODO: swap + assign operator?

  // ACCESS

  /**
   * @brief Returns a thread safe copy of this object. The copy is lighter than a real copy, as the complex is passed
   * by pointer and gives access to all const method and persistence related methods. But the returned object will be
   * invalidated if this object is destroyed.
   */
  Thread_safe weak_copy() const { return Thread_safe(*this); }

  /**
   * @brief Returns the number of generators in the stored module.
   */
  Index get_number_of_cycle_generators() const { return complex_.get_number_of_cycle_generators(); }

  /**
   * @brief Returns the number of parameters of the stored filtration. If the module is empty, the number returned is 0.
   */
  Index get_number_of_parameters() const { return complex_.get_number_of_parameters(); }

  // // only used for scc io for now
  // const Complex& get_chain_complex() const { return complex_; }

  /**
   * @brief Returns a const reference to the current permutation map, indicating in which order are the generators
   * with respect to the current slice (i.e., \$f order[i] \$f corresponds to the index in the complex of the
   * \$f i^{th} \$f generator in the filtration represented by the slice). It will be initialized with
   * @ref initialize_persistence_computation.
   *
   * If `ignoreInf` was true when calling @ref initialize_persistence_computation, indices of generators at infinity
   * are not stored in the container. That means that the size can be smaller than what
   * @ref get_number_of_cycle_generators returns.
   */
  const std::vector<Index>& get_current_order() const { return generatorOrder_; }

  /**
   * @brief Returns a const reference to the current slice. It can be initialized or updated with @ref set_slice
   * and @ref push_to.
   */
  const std::vector<T>& get_slice() const { return slice_; }

  /**
   * @brief Returns a const reference to the class computing the persistence of the current slice. It will be
   * initialized with @ref initialize_persistence_computation.
   */
  const Persistence& get_persistence_algorithm() const { return persistence_; }

  /**
   * @brief Returns two filtration values representing respectively the greatest common lower bound of all filtration
   * values in the filtration and the lowest common upper bound of them.
   */
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

  /**
   * @brief Returns a const reference to the filtration value container. A filtration value at index \$f i \$f 
   * correspond to the filtration value associated to the generators at index \$f i \$f.
   */
  const typename Complex::Filtration_value_container& get_filtration_values() const
  {
    return complex_.get_filtration_values();
  }

  /**
   * @brief Returns a reference to the filtration value container. A filtration value at index \$f i \$f 
   * correspond to the filtration value associated to the generators at index \$f i \$f.
   *
   * @warning The container is not const such that the user can easily modify/update a filtration value. But do not
   * modify the size of the container, nor the number of parameters.
   */
  typename Complex::Filtration_value_container& get_filtration_values() { return complex_.get_filtration_values(); }

  /**
   * @brief Returns a const reference to the filtration value associated to the generator at index \$f i \$f.
   */
  const Filtration_value& get_filtration_value(Index i) const { return complex_.get_filtration_values()[i]; }

  /**
   * @brief Returns a reference to the filtration value associated to the generator at index \$f i \$f.
   *
   * @warning The value is not const such that the user can easily modify/update the filtration value. But do not
   * modify the number of parameters.
   */
  Filtration_value& get_filtration_value(Index i) { return complex_.get_filtration_values()[i]; }

  /**
   * @brief Returns a const reference to the dimension container. A value at index \$f i \$f corresponds to the
   * dimension of the generator at index \$f i \$f.
   */
  const std::vector<Dimension>& get_dimensions() const { return complex_.get_dimensions(); }

  /**
   * @brief Returns the dimension of the generator at index \$f i \$f.
   */
  Dimension get_dimension(Index i) const { return complex_.get_dimensions()[i]; }

  /**
   * @brief Returns the maximal dimension of a generator in the module.
   */
  Dimension get_max_dimension() const { return complex_.get_max_dimension(); }

  /**
   * @brief Returns a const reference to the boundary container. The element at index \$f i \$f corresponds to the
   * boundary of the generator at index \$f i \$f.
   */
  const typename Complex::Boundary_container& get_boundaries() const { return complex_.get_boundaries(); }

  /**
   * @brief Returns the boundary of the generator at index \$f i \$f.
   */
  const typename Complex::Boundary& get_boundary(Index i) const { return complex_.get_boundaries()[i]; }

  // // TODO: only used to print info in python, so put in some interface instead
  // std::string to_str() const
  // {
  //   std::stringstream stream;
  //   stream << *this;
  //   return stream.str();
  // }

  // MODIFIERS

  /**
   * @brief Sets the current slice, that is the 1-parameter filtration values associated to each generator on that line.
   * The value at \$f slice[i] \$f has to corresponds to the value for the generator at index \$f i \$f.
   * One can also sets the slice directly from the line with @ref push_to.
   * 
   * @tparam Array Container which can be converted into a vector of `T`.
   */
  template <class Array = std::initializer_list<T>>
  void set_slice(const Array& slice)
  {
    slice_ = slice;
  }

  /**
   * @brief Sets the current slice by computing the 1-parameter filtration values fo each generator on the given line.
   * 
   * @tparam Line_like Any type convertible to a @ref Line class. Default value: `std::initializer_list<T>`.
   */
  template <class Line_like = std::initializer_list<T>>
  void push_to(const Line_like& line)
  {
    _push_to(complex_, Line<typename Line_like::value_type>(line));
  }

  /**
   * @brief Sets the current slice by computing the 1-parameter filtration values fo each generator on the given line.
   * 
   * @tparam U Template parameter of the given line.
   */
  template <class U>
  void push_to(const Line<U>& line)
  {
    _push_to(complex_, line);
  }

  /**
   * @brief Removes completely from the module all generator of dimension strictly higher than given. All
   * initializations are invalidated, so the slice has to be reset and @ref initialize_persistence_computation recalled.
   *
   * @warning If the internal complex was not ordered by dimension, the complex is sorted before pruning.
   * So, the indexing changes afterwards.
   * 
   * @param maxDim Maximal dimension to keep.
   */
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
                                       int degree = -1,
                                       bool rivetCompatible = false,
                                       bool ignoreLastGenerators = false,
                                       bool stripComments = false,
                                       bool reverse = false)
  {
    const Complex& cpx =
        slicer.complex_.is_ordered_by_dimension() ? slicer.complex_ : build_permuted_complex(slicer.complex_).first;
    write_complex_to_scc_file<typename Slicer::Filtration_value>(
        outFilePath, cpx, degree, rivetCompatible, ignoreLastGenerators, stripComments, reverse);
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

  template <class U>
  void _push_to(const Complex& complex, const Line<U>& line)
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
