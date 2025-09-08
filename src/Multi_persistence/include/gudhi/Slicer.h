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

#include <array>
#include <initializer_list>
#include <type_traits>
#include <numeric>  //std::iota
#include <utility>  //std::move
#include <cstdint>  //std::int32_t
#include <vector>
#include <ostream>
// #include <sstream>  //std::stringstream, to remove when to_str gets removed

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
  using Bar = Gudhi::persistence_matrix::Persistence_interval<Dimension, Value>; /**< Bar type. */
  /**
   * @brief Barcode type. A vector of @ref Bar, a tuple like structure containing birth, death and dimension of a bar.
   */
  template <typename Value = T>
  using Barcode = std::vector<Bar<Value>>;
  /**
   * @brief Flat barcode type. All bars are represented by a birth and a death value stored in arrays of size 2.
   */
  template <typename Value = T>
  using Flat_barcode = std::vector<std::array<Value,2> >;
  /**
   * @brief Barcode ordered by dimension type. A vector which has at index \f$ d \f$ the @ref Barcode of dimension
   * \f$ d \f$.
   */
  template <typename Value = T>
  using Multi_dimensional_barcode = std::vector<Barcode<Value>>;
  /**
   * @brief Flat barcode ordered by dimension type. A vector which has at index \f$ d \f$ the @ref Flat_barcode of
   * dimension \f$ d \f$.
   */
  template <typename Value = T>
  using Multi_dimensional_flat_barcode = std::vector<Flat_barcode<Value>>;
  using Cycle = std::vector<Index>;               /**< Cycle type. */
  using Thread_safe = Thread_safe_slicer<Slicer>; /**< Thread safe slicer type. */

  // CONSTRUCTORS

  /**
   * @brief Default constructor. Constructs an empty slicer.
   */
  Slicer() = default;

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
   * @brief Copy constructor. Persistence computation initialization is not updated.
   */
  Slicer(const Slicer& other)
      : complex_(other.complex_), slice_(other.slice_), generatorOrder_(other.generatorOrder_), persistence_()
  {}

  /**
   * @brief Move constructor. Persistence computation initialization is not updated.
   */
  Slicer(Slicer&& other) noexcept
      : complex_(std::move(other.complex_)),
        slice_(std::move(other.slice_)),
        generatorOrder_(std::move(other.generatorOrder_)),
        persistence_()
  {}

  ~Slicer() = default;

  /**
   * @brief Assign operator. Persistence computation initialization is not updated.
   */
  Slicer& operator=(const Slicer& other)
  {
    complex_ = other.complex_;
    slice_ = other.slice_;
    generatorOrder_ = other.generatorOrder_;
    persistence_.reset();

    return *this;
  }

  /**
   * @brief Move assign operator. Persistence computation initialization is not updated.
   */
  Slicer& operator=(Slicer&& other) noexcept
  {
    complex_ = std::move(other.complex_);
    slice_ = std::move(other.slice_);
    generatorOrder_ = std::move(other.generatorOrder_);
    persistence_.reset();

    return *this;
  }

  // TODO: swap ?

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
   * @brief Returns a reference to the current slice. It can also be initialized or updated with @ref set_slice
   * and @ref push_to.
   */
  std::vector<T>& get_slice() { return slice_; }

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
   * @tparam Array Container with a begin() and end() method and whose element can be converted into `T`.
   */
  template <class Array = std::initializer_list<T>>
  void set_slice(const Array& slice)
  {
    slice_ = std::vector<T>(slice.begin(), slice.end());
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
    persistence_.reset();
  }

  /**
   * @brief Projects all filtration values into the given grid. If @p coordinate is false, the entries are set to
   * the nearest upper bound value with the same parameter in the grid. Otherwise, the entries are set to the indices
   * of those nearest upper bound values.
   * An index \f$ i \f$ of the grid corresponds to the same parameter as the index \f$ i \f$ in a generator of the
   * filtration value. The internal vectors correspond to the possible values of the parameters, ordered by increasing
   * value, forming therefore all together a 2D grid.
   *
   * @param grid Vector of vector with size at least number of filtration parameters.
   * @param coordinate If true, the values are set to the coordinates of the projection in the grid. If false,
   * the values are set to the values at the coordinates of the projection. Default value: true.
   */
  void coarsen_on_grid(const std::vector<std::vector<T>>& grid, bool coordinate = true)
  {
    complex_.coarsen_on_grid(grid, coordinate);
  }

  // PERSISTENCE

  /**
   * @brief Returns true if and only if @ref initialize_persistence_computation was properly called.
   */
  [[nodiscard]] bool persistence_computation_is_initialized() const { return persistence_.is_initialized(); }

  /**
   * @brief Initializes the persistence computation of the current slice. If the slice was not set properly as
   * a valid 1-dimensional filtration, the behaviour is undefined.
   *
   * @param ignoreInf If true, all cells at infinity filtration values are ignored for the initialization, resulting
   * potentially in less storage use and better performance. But note that this can be problematic with the use of
   * @ref vineyard_update. Default value: true.
   */
  void initialize_persistence_computation(const bool ignoreInf = true)
  {
    _initialize_persistence_computation(complex_, ignoreInf);
  }

  /**
   * @brief After the persistence computation was initialized for a slice and the slice changes, this method can
   * update everything necessary for the barcode without re-computing everything from scratch (contrary to
   * @ref initialize_persistence_computation). Furthermore, it guarantees that the new barcode will "match" the
   * precedent one. TODO: explain exactly what it means and how to do the matching.
   * The method will have better performance if the complex is ordered by dimension.
   *
   * Only available if PersistenceAlgorithm::is_vine is true.
   *
   * @pre @ref initialize_persistence_computation has to be called at least once before.
   *
   * @warning If `ignoreInf` was set to true when initializing the persistence computation, any update of the slice has
   * to keep at infinity the boundaries which were before, otherwise the behaviour is undefined (it will throw with
   * high probability).
   */
  void vineyard_update()
  {
    static_assert(Persistence::is_vine, "vineyard_update() not enabled by the chosen PersistenceAlgorithm class.");

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

  /**
   * @brief Returns the barcode of the current slice. The barcode format will change depending on the template values.
   *
   * @pre @ref initialize_persistence_computation has to be called at some point before.
   *
   * @tparam byDim If true, the barcode is returned as @ref Multi_dimensional_barcode, otherwise as @ref Barcode.
   * @tparam Value Type of the birth and death values.
   * @param maxDim Maximal dimension to be included in the barcode. If negative, all dimensions are included.
   * Default value: -1.
   */
  template <bool byDim = true, typename Value = T, bool idx = false>
  std::conditional_t<byDim, Multi_dimensional_barcode<Value>, Barcode<Value>> get_barcode(int maxDim = -1)
  {
    if (maxDim < 0) maxDim = get_max_dimension();
    if constexpr (byDim) {
      return _get_barcode_by_dim<idx, Value>(maxDim);
    } else {
      return _get_barcode<idx, Value>(maxDim);
    }
  }

  /**
   * @brief Returns the barcode of the current slice. The barcode format will change depending on the template values.
   *
   * @pre @ref initialize_persistence_computation has to be called at some point before.
   *
   * @tparam byDim If true, the barcode is returned as @ref Multi_dimensional_flat_barcode, otherwise as
   * @ref Flat_barcode.
   * @tparam Value Type of the birth and death values.
   * @param maxDim Maximal dimension to be included in the barcode. If negative, all dimensions are included.
   * Default value: -1.
   */
  template <bool byDim = false, typename Value = T, bool idx = false>
  std::conditional_t<byDim, Multi_dimensional_flat_barcode<Value>, Flat_barcode<Value>> get_flat_barcode(
      int maxDim = -1)
  {
    if (maxDim < 0) maxDim = get_max_dimension();
    if constexpr (byDim) {
      return _get_flat_barcode_by_dim<idx, Value>(maxDim);
    } else {
      return _get_flat_barcode<idx, Value>(maxDim);
    }
  }

  /**
   * @brief Returns the representative cycles of the current slice. All cycles of dimension \f$ d \f$ are stored at
   * index \f$ d \f$ of the returned vector. A cycle is represented by a vector of boundary indices. That is, the index
   * \f$ i \f$ in a cycle represents the cell which boundary can be retrieved by @ref get_boundary "get_boundary(i)".
   *
   * Only available if PersistenceAlgorithm::has_rep_cycles is true.
   *
   * @pre @ref initialize_persistence_computation has to be called at least once before.
   *
   * @param update If true, updates the stored representative cycles, otherwise just returns the container in its
   * current state. So should be true at least the first time the method is used.
   */
  std::vector<std::vector<Cycle>> get_representative_cycles(bool update = true)
  {
    return _get_representative_cycles(complex_, update);
  }

  // FRIENDS

  /**
   * @brief Builds a new slicer by reordering the cells in the complex of the given slicer with the given permutation
   * map.
   */
  friend Slicer build_permuted_slicer(const Slicer& slicer, const std::vector<Index>& permutation)
  {
    GUDHI_CHECK(permutation.size() > slicer.get_number_of_cycle_generators(),
                "Too many elements in permutation vector.");
    return Slicer(build_permuted_complex(slicer.complex_, permutation));
  }

  /**
   * @brief Builds a new slicer by reordering the cells in the complex of the given slicer the same way than
   * @ref Multi_parameter_filtered_complex::sort_by_dimension_co_lexicographically. Returns a pair with the new slicer
   * as first element and the permutation map used as second element.
   */
  friend std::pair<Slicer, std::vector<Index>> build_permuted_slicer(const Slicer& slicer)
  {
    auto [complex, permutation] = build_permuted_complex(slicer.complex_);
    return std::make_pair(Slicer(std::move(complex)), std::move(permutation));
  }

  /**
   * @brief Builds a new slicer from the given one by projecting its filtration values on a grid.
   * See @ref coarsen_on_grid with the paramater `coordinate` at true.
   */
  friend auto build_slicer_coarsen_on_grid(const Slicer& slicer, const std::vector<std::vector<T>> grid)
  {
    using return_filtration_value = decltype(std::declval<Filtration_value>().template as_type<std::int32_t>());
    using return_complex = decltype(build_complex_coarsen_on_grid(slicer.complex_, grid));
    using return_pers = typename Persistence::template As_type<return_complex>;
    return Slicer<return_filtration_value, return_pers>(build_complex_coarsen_on_grid(slicer.complex_, grid));
  }

  /**
   * @brief Builds a new slicer using @ref Projective_cover_kernel. TODO: explain what that means.
   */
  friend Slicer build_slicer_from_projective_cover_kernel(const Slicer& slicer, Dimension dim)
  {
    Projective_cover_kernel<Filtration_value> kernel(slicer.complex_, dim);
    return Slicer(kernel.create_complex());
  }

  /**
   * @brief Writes the given slicer into a file with scc format. Assumes that every index appearing in a boundary of
   * the complex corresponds to an existing index in the complex (for example, the lowest dimension has always empty
   * boundaries).
   * See @ref build_slicer_from_scc_file to build a slicer from a scc format file.
   *
   * @param outFilePath Path with file name into which to write.
   * @param slicer Slicer to write. Every index appearing in a boundary of the complex has to correspond to an
   * existing index in the underlying complex.
   * @param degree TODO Default value: -1.
   * @param rivetCompatible Set to true if the written file has to be Rivet compatible. Note that Rivet only accepts
   * bi-filtrations. Default value: false.
   * @param ignoreLastGenerators Set to true, if the generators with last dimension in the list should be ignored
   * (maximal dimension by default, minimal dimension if `reverse` is true). Default value: false.
   * @param stripComments Set to true, if no comment should be written in the file (comments are lines starting with `#`
   * and which are ignored when read). Default value: false.
   * @param reverse Set to true if the generators should be written in increasing order of dimension instead of
   * decreasing. Default value: false.
   */
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

  /**
   * @brief Outstream operator.
   */
  friend std::ostream& operator<<(std::ostream& stream, Slicer& slicer)
  {
    stream << "-------------------- Slicer \n";

    stream << "--- Filtered complex \n";
    stream << slicer.complex_;

    stream << "--- Order \n";
    stream << "{";
    for (const auto& idx : slicer.generatorOrder_) stream << idx << ", ";
    stream << "}" << '\n';

    stream << "--- Current slice filtration\n";
    stream << "{";
    for (const auto& val : slicer.slice_) stream << val << ", ";
    if (!slicer.slice_.empty()) stream << "\b" << "\b";
    stream << "}" << '\n';

    stream << "--- PersBackend \n";
    stream << slicer.persistence_;

    return stream;
  }

 protected:
  friend Thread_safe;  // Thread_safe will use the "_*" methods below instead of "*".

  // For ThreadSafe version
  Slicer(const std::vector<T>& slice, const std::vector<Index>& generatorOrder, const Persistence& persistence)
      : complex_(), slice_(slice), generatorOrder_(generatorOrder), persistence_(persistence, generatorOrder_)
  {}

  Slicer(std::vector<T>&& slice, std::vector<Index>&& generatorOrder, Persistence&& persistence)
      : complex_(),
        slice_(std::move(slice)),
        generatorOrder_(std::move(generatorOrder)),
        persistence_(std::move(persistence), generatorOrder_)
  {}

  template <class U>
  void _push_to(const Complex& complex, const Line<U>& line)
  {
    const auto& filtrationValues = complex.get_filtration_values();
    for (Index i = 0U; i < filtrationValues.size(); i++) {
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
    static_assert(Persistence::has_rep_cycles,
                  "Representative cycles not enabled by the chosen PersistenceAlgorithm class.");

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

 private:
  Complex complex_;      /**< Complex storing all boundaries, filtration values and dimensions. */
  std::vector<T> slice_; /**< Filtration values of the current slice. The indices corresponds to those in complex_. */
  std::vector<Index> generatorOrder_; /**< Permutation map from current slice index to complex index. */
  Persistence persistence_;           /**< Class for persistence computations. */

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

  template <bool idx, class Interval, typename Value>
  void _retrieve_interval(const Interval& bar, Dimension& dim, Value& birth, Value& death)
  {
    const Value inf = Gudhi::multi_filtration::MF_T_inf<Value>;
    dim = bar.dim;
    if constexpr (idx) {
      birth = bar.birth;
      death = -1;
      if (bar.death != Persistence::nullDeath) death = bar.death;
    } else {
      birth = slice_[bar.birth];
      death = inf;
      if (bar.death != Persistence::nullDeath) death = slice_[bar.death];
      if (!(birth <= death)) {
        birth = inf;
        death = inf;
      }
    }
  }

  template <bool idx, typename Value>
  Barcode<Value> _get_barcode(int maxDim)
  {
    auto barcodeIndices = persistence_.get_barcode();
    Barcode<Value> out(barcodeIndices.size());
    Index i = 0;
    for (const auto& bar : barcodeIndices) {
      if (bar.dim <= maxDim) {
        _retrieve_interval<idx>(bar, out[i].dim, out[i].birth, out[i].death);
        ++i;
      }
    }
    out.resize(i);
    return out;
  }

  template <bool idx, typename Value>
  Multi_dimensional_barcode<Value> _get_barcode_by_dim(int maxDim)
  {
    // TODO: This doesn't allow for negative dimensions
    // Hannah: not sure what this comment means ?
    Multi_dimensional_barcode<Value> out(maxDim + 1);
    Value birth, death;
    Dimension dim;
    for (const auto& bar : persistence_.get_barcode()) {
      if (bar.dim <= maxDim) {
        _retrieve_interval<idx>(bar, dim, birth, death);
        out[dim].emplace_back(birth, death, dim);
      }
    }
    return out;
  }

  template <bool idx, typename Value>
  Flat_barcode<Value> _get_flat_barcode(int maxDim)
  {
    auto barcodeIndices = persistence_.get_barcode();
    Flat_barcode<Value> out(barcodeIndices.size());
    Index i = 0;
    Dimension dim;  // dummy
    for (const auto& bar : barcodeIndices) {
      if (bar.dim <= maxDim) {
        _retrieve_interval<idx>(bar, dim, out[i][0], out[i][1]);
        ++i;
      }
    }
    out.resize(i);
    return out;
  }

  template <bool idx, typename Value>
  Multi_dimensional_flat_barcode<Value> _get_flat_barcode_by_dim(int maxDim)
  {
    Multi_dimensional_flat_barcode<Value> out(maxDim + 1);
    Value birth, death;
    Dimension dim;
    for (const auto& bar : persistence_.get_barcode()) {
      if (bar.dim <= maxDim) {
        _retrieve_interval<idx>(bar, dim, birth, death);
        out[dim].emplace_back(std::array<Value, 2>{birth, death});
      }
    }
    return out;
  }
};

}  // namespace multi_persistence
}  // namespace Gudhi

#endif  // MP_SLICER_H_INCLUDED
