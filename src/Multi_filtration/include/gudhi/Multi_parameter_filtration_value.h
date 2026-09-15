/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber, David Loiseaux
 *
 *    Copyright (C) 2026 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file Multi_parameter_filtration_value.h
 * @author Hannah Schreiber, David Loiseaux
 * @brief Contains the @ref Gudhi::multi_filtration::Multi_parameter_filtration_value class.
 */

#ifndef MF_MULTI_PARAMETER_FILTRATION_VALUE_H_
#define MF_MULTI_PARAMETER_FILTRATION_VALUE_H_

#include <algorithm>    //std::lower_bound
#include <cmath>        //std::isnan, std::min, std::abs
#include <cstddef>      //std::size_t
#include <cstdint>      //std::int32_t, std::uint8_t
#include <cstring>      //memcpy
#include <iterator>     //std::distance
#include <ostream>      //std::ostream
#include <limits>       //std::numerical_limits
#include <stdexcept>    //std::logic_error
#include <type_traits>  //std::is_arithmetic
#include <utility>      //std::swap, std::move
#include <vector>
#include <initializer_list>

#ifdef GUDHI_USE_TBB
#include <oneapi/tbb/parallel_for.h>
#endif

#include <gudhi/Debug_utils.h>
#include <gudhi/serialization_utils.h>
#include <gudhi/Multi_filtration/multi_filtration_utils.h>
#include <gudhi/Multi_filtration/multi_filtration_conversions.h>

namespace Gudhi::multi_filtration {

// declaration needed pre C++20
template <typename OutValue, typename U>
void compute_coordinates_in_grid();
template <typename U, class RandomAccessArray>
void evaluate_coordinates_in_grid();
template <class RandomAccessArray>
void evaluate_coordinates_in_grid();
template <bool inverse>
bool is_strict_less_than_lexicographically();
template <bool inverse>
bool is_less_or_equal_than_lexicographically();
template <bool greatest>
bool factorize();

/**
 * @class Multi_parameter_filtration_value Multi_parameter_filtration_value.h gudhi/Multi_parameter_filtration_value.h
 * @ingroup multi_filtration
 *
 * @brief Class encoding the different generators, i.e., apparition times, of a \f$ k \f$-critical
 * \f$\mathbb R^n\f$-filtration value. E.g., the filtration value of a simplex, or, of the algebraic generator of a
 * module presentation.
 * Implements the concept @ref FiltrationValue of the @ref Gudhi::Simplex_tree and the concept
 * @ref Gudhi::multi_persistence::MultiFiltrationValue.
 *
 * @details
 * Multi-critical filtrations are filtrations such that the lifetime of each object is union of positive cones in
 * \f$\mathbb R^n\f$, e.g.,
 *  - \f$ \{ x \in \mathbb R^2 : x \ge (1,2)\} \cap \{ x \in \mathbb R^2 : x \ge (2,1)\} \f$ is finitely critical,
 *    and more particularly 2-critical, while
 *  - \f$ \{ x \in \mathbb R^2 : x \ge \mathrm{epigraph}(y \mapsto e^{-y})\} \f$ is not.
 *
 * Overloads `std::numeric_limits` can behave differently depending on the template parameters.
 *
 * @tparam StoragePolicy Has to follow the concept @ref StoragePolicy. Decides on how the underlying data is stored
 * and handled.
 * @tparam Co If `true`, reverses the poset order, i.e., the order \f$ \le \f$  in \f$ \mathbb R^n \f$ becomes
 * \f$ \ge \f$. That is, the positive cones representing a lifetime become all negative instead.
 * @tparam Ensure1Criticality If `true`, the methods ensure that the filtration value is always 1-critical by throwing
 * or refusing to compile if a modification increases the number of generators.
 */
template <class StoragePolicy, bool Co = false, bool Ensure1Criticality = false>
class Multi_parameter_filtration_value {
 public:
  /**
   * @brief Type of an element in the filtration value.
   * See documentation of corresponding StoragePolicy type or concept @ref StoragePolicy for more details.
   */
  using value_type = typename StoragePolicy::value_type;
  /**
   * @brief Underlying container for the elements in the filtration value.
   * See documentation of corresponding StoragePolicy type or concept @ref StoragePolicy for more details.
   */
  using Underlying_container = typename StoragePolicy::Underlying_container;
  /**
   * @brief Size type.
   * See documentation of corresponding StoragePolicy type or concept @ref StoragePolicy for more details.
   */
  using size_type = typename StoragePolicy::size_type;
  /**
   * @brief Reference type.
   * See documentation of corresponding StoragePolicy type or concept @ref StoragePolicy for more details.
   */
  using reference = typename StoragePolicy::reference;
  /**
   * @brief Const reference type. Note that this type could not be a reference.
   * See documentation of corresponding StoragePolicy type or concept @ref StoragePolicy for more details.
   */
  using const_reference = typename StoragePolicy::const_reference;
  /**
   * @brief Iterator type. At least LegacyForwardIterator.
   * See documentation of corresponding StoragePolicy type or concept @ref StoragePolicy for more details.
   */
  using iterator = typename StoragePolicy::iterator;
  /**
   * @brief Const iterator type. At least LegacyForwardIterator.
   * See documentation of corresponding StoragePolicy type or concept @ref StoragePolicy for more details.
   */
  using const_iterator = typename StoragePolicy::const_iterator;
  /**
   * @brief First template parameter, storage policy.
   */
  using Storage_policy = StoragePolicy;

  /**
   * @brief Plus infinity value of an entry of the filtration value.
   */
  constexpr static const value_type T_inf = StoragePolicy::template T_inf<>;

  /**
   * @brief Minus infinity value of an entry of the filtration value.
   */
  constexpr static const value_type T_m_inf = StoragePolicy::template T_m_inf<>;

  /**
   * @brief True if and only if @ref inf compiles and does not throw.
   */
  constexpr static const bool has_infinity = StoragePolicy::template has_infinity<Co>;

  /**
   * @brief True if and only if @ref minus_inf compiles and does not throw.
   */
  constexpr static const bool has_minus_infinity = StoragePolicy::template has_minus_infinity<Co>;

  /**
   * @brief True if and only if @ref nan compiles and does not throw.
   */
  constexpr static const bool has_quiet_NaN = StoragePolicy::has_quiet_NaN;

  // CONSTRUCTORS

  /**
   * @brief Default constructor. Builds filtration value with one generator and given number of parameters.
   * If Co is false, all values are at +inf, if Co is true, all values are at -inf.
   *
   * @param numberOfParameters If negative, takes the default value instead. Default value: 2.
   */
  Multi_parameter_filtration_value(int numberOfParameters = 2)
      : generators_(numberOfParameters < 0 ? 2 : numberOfParameters, _get_default_value()) {}

  /**
   * @brief Builds filtration value with one generator and given number of parameters.
   * All values are initialized at the given value.
   *
   * @param numberOfParameters If negative, is set to 2 instead.
   * @param value Initialization value for every value in the generator.
   */
  Multi_parameter_filtration_value(int numberOfParameters, value_type value)
      : generators_(numberOfParameters < 0 ? 2 : numberOfParameters, value) {}

  /**
   * @brief Builds filtration value with one generator that is initialized with the given range. The number of
   * parameters are therefore deduced from the length of the range.
   *
   * @param range Values of the generator.
   */
  Multi_parameter_filtration_value(std::initializer_list<value_type> range) : generators_(range.begin(), range.end()) {}

  /**
   * @brief Builds filtration value with one generator that is initialized with the given range. The number of
   * parameters are therefore deduced from the length of the range.
   *
   * @tparam ValueRange Range of types convertible to `T`. Should have a begin() and end() method.
   * @param range Values of the generator.
   */
  template <class ValueRange, class = std::enable_if_t<detail::RangeTraits<ValueRange>::has_begin>>
  Multi_parameter_filtration_value(const ValueRange &range) : generators_(range.begin(), range.end()) {}

  /**
   * @brief Builds filtration value with one generator that is initialized with the given range. The range is
   * determined from the two given iterators. The number of parameters are therefore deduced from the distance
   * between the two.
   *
   * @tparam Iterator Iterator type that has to satisfy the requirements of standard LegacyInputIterator and
   * dereferenced elements have to be convertible to `T`.
   * @param itBegin Iterator pointing to the start of the range.
   * @param itEnd Iterator pointing to the end of the range.
   */
  template <class Iterator, class = std::enable_if_t<!std::is_arithmetic_v<Iterator>>>
  Multi_parameter_filtration_value(Iterator itBegin, Iterator itEnd) : generators_(itBegin, itEnd) {}

  /**
   * @brief Builds filtration value with given number of parameters and values from the given range. Lets \f$ p \f$
   * be the number of parameters. The \f$ p \f$ first elements of the range have to correspond to the first generator,
   * the \f$ p \f$ next elements to the second generator and so on... So the length of the range has to be a multiple
   * of \f$ p \f$ and the number of generators will be \f$ length / p \f$. The range is represented by two iterators.
   *
   * @tparam Iterator Iterator type that has to satisfy the requirements of standard LegacyInputIterator and
   * dereferenced elements have to be convertible to `T`.
   * @param itBegin Iterator pointing to the start of the range.
   * @param itEnd Iterator pointing to the end of the range.
   * @param numberOfParameters Negative values are associated to 0.
   */
  template <class Iterator, class = std::enable_if_t<!std::is_arithmetic_v<Iterator>>>
  Multi_parameter_filtration_value(Iterator itBegin, Iterator itEnd, int numberOfParameters)
      : generators_(itBegin, itEnd, numberOfParameters < 0 ? 0 : numberOfParameters) {
    if constexpr (Ensure1Criticality) {
      if (generators_.num_generators() != 1)
        throw std::invalid_argument("Multiparameter filtration value is not 1-critical.");
    }
  }

  /**
   * @brief Builds filtration value with given number of parameters and by copying the given container as is. The
   * format of that container depends on the template parameter @ref StoragePolicy.
   *
   * @param generators Values.
   * @param numberOfParameters Negative values are associated to 0.
   */
  Multi_parameter_filtration_value(const Underlying_container &generators, int numberOfParameters)
      : generators_(generators, numberOfParameters < 0 ? 0 : numberOfParameters) {
    if constexpr (Ensure1Criticality) {
      if (generators_.num_generators() != 1)
        throw std::invalid_argument("Multiparameter filtration value is not 1-critical.");
    }
  }

  /**
   * @brief Builds filtration value with given number of parameters and by moving the given container into it. The
   * format of that container depends on the template parameter @ref StoragePolicy.
   *
   * @param generators Values to move.
   * @param numberOfParameters Negative values are associated to 0.
   */
  Multi_parameter_filtration_value(Underlying_container &&generators, int numberOfParameters)
      : generators_(std::move(generators), numberOfParameters < 0 ? 0 : numberOfParameters) {
    if constexpr (Ensure1Criticality) {
      if (generators_.num_generators() != 1)
        throw std::invalid_argument("Multiparameter filtration value is not 1-critical.");
    }
  }

  /**
   * @brief Builds filtration value by copying given @ref StoragePolicy.
   */
  Multi_parameter_filtration_value(const StoragePolicy &generators) : generators_(generators) {}

  /**
   * @brief Builds filtration value by moving given @ref StoragePolicy.
   */
  Multi_parameter_filtration_value(StoragePolicy &&generators) : generators_(std::move(generators)) {}

  /**
   * @brief Copy constructor.
   */
  Multi_parameter_filtration_value(const Multi_parameter_filtration_value &other) : generators_(other.generators_) {}

  /**
   * @brief Move constructor.
   */
  Multi_parameter_filtration_value(Multi_parameter_filtration_value &&other) noexcept
      : generators_(std::move(other.generators_)) {}

  ~Multi_parameter_filtration_value() = default;

  /**
   * @brief Assign operator.
   */
  Multi_parameter_filtration_value &operator=(const Multi_parameter_filtration_value &other) = default;

  /**
   * @brief Move assign operator.
   */
  Multi_parameter_filtration_value &operator=(Multi_parameter_filtration_value &&other) noexcept {
    generators_ = std::move(other.generators_);
    return *this;
  }

  /**
   * @brief Swap operator.
   */
  friend void swap(Multi_parameter_filtration_value &f1, Multi_parameter_filtration_value &f2) noexcept {
    swap(f1.generators_, f2.generators_);
  }

  // CONVERTERS

  /**
   * @brief Returns a copy as a @ref Multi_parameter_filtration_value with given template arguments.
   *
   * @tparam OtherStoragePolicy New `StoragePolicy` template of @ref Multi_parameter_filtration_value.
   * @tparam OtherCo New `Co` template of @ref Multi_parameter_filtration_value.
   * @tparam OtherEnsure1Criticality New `Ensure1Criticality` template of @ref Multi_parameter_filtration_value.
   *
   * @throw std::logic_error When trying to convert a filtration value with more than 1 generator into a
   * @ref Multi_parameter_filtration_value with `OtherEnsure1Criticality` template argument at `true`.
   */
  template <class OtherStoragePolicy, bool OtherCo = Co, bool OtherEnsure1Criticality = Ensure1Criticality,
            class = std::enable_if_t<detail::RangeTraits<OtherStoragePolicy>::is_storage_policy>>
  Multi_parameter_filtration_value<OtherStoragePolicy, OtherCo, OtherEnsure1Criticality> as_type() const {
    if constexpr (OtherEnsure1Criticality) {
      if (num_generators() > 1)
        throw std::logic_error("Cannot convert a k-critical filtration value into a 1-critical one.");
    }
    OtherStoragePolicy gens = Gudhi::multi_filtration::as_type<OtherStoragePolicy, OtherCo>(generators_);
    return Multi_parameter_filtration_value<OtherStoragePolicy, OtherCo, OtherEnsure1Criticality>(std::move(gens));
  }

  /**
   * @brief If `U` is a native arithmetic type, returns a copy by casting every filtration value element in that type.
   * Otherwise, `U` has to be @ref Multi_parameter_filtration_value with the desired template arguments, and
   * the method returns a copy into that new format.
   */
  template <typename U,
            class = std::enable_if_t<std::is_arithmetic_v<U> || detail::RangeTraits<U>::is_multi_filtration>>
  auto as_type() const {
    if constexpr (std::is_arithmetic_v<U>) {
      return as_type<typename StoragePolicy::template As_type<U>, Co, Ensure1Criticality>();
    } else {
      return as_type<typename U::Storage_policy, U::has_negative_cones(), U::ensures_1_criticality()>();
    }
  }

  /**
   * @brief Returns a copy of this with given number of parameters and generators. If the given number of parameters
   * \f$ p \f$ is smaller than the current one, only the \f$ p \f$ first parameters are copied and if \f$ p \f$ is
   * greater than the current one, then additional parameters are added at the end using the given default value.
   * Same for the given number of generators (similar to @ref set_num_generators "").
   *
   * @warning All new generators or parameters will be set to -infinity (`Co` is true) or infinity (`Co` is false).
   * That is, the new filtration value is potentially not minimal anymore. So,
   * if @ref StoragePolicy::has_minimal_set_representation is true, make sure to fill them with real generators or to
   * call @ref simplify before using other methods as most methods will have an undefined behaviour for that case,
   * if the set of generators is not minimal or sorted.
   * 
   * @tparam U Arithmetic type @ref value_type can convert into. Default: @ref value_type.
   */
  template <typename U = value_type>
  auto copy(size_type numberOfParameters, size_type numberOfGenerators) const {
    using SP = typename StoragePolicy::template As_type<U>;
    SP sp = generators_.template copy<U>(numberOfParameters, numberOfGenerators,
                                         Co ? SP::template T_m_inf<> : SP::template T_inf<>);
    return Multi_parameter_filtration_value<SP, Co, Ensure1Criticality>(std::move(sp));
  }

  // ACCESS

  /**
   * @brief Returns reference to underlying @ref StoragePolicy.
   */
  StoragePolicy &get_underlying_policy() { return generators_; }

  /**
   * @brief Returns const reference to underlying @ref StoragePolicy.
   */
  const StoragePolicy &get_underlying_policy() const { return generators_; }

  /**
   * @brief Returns reference to underlying value container.
   */
  Underlying_container &get_underlying_container() { return generators_.get_underlying_container(); }

  /**
   * @brief Returns const reference to underlying value container.
   */
  const Underlying_container &get_underlying_container() const { return generators_.get_underlying_container(); }

  /**
   * @brief Returns reference to value of parameter `p` of generator `g`.
   */
  reference operator()(size_type g, size_type p) { return generators_(g, p); }

  /**
   * @brief Returns const reference or copy to value of parameter `p` of generator `g`.
   */
  const_reference operator()(size_type g, size_type p) const { return generators_(g, p); }

  /**
   * @brief Let \f$ g \f$ be the first value in `indices` and \f$ p \f$ the second value.
   * Returns reference to value of parameter \f$ p \f$ of generator \f$ g \f$.
   *
   * @tparam IndexRange Range with a begin() and size() method.
   * @param indices Range with at least two elements. The first element should correspond to the generator number and
   * the second element to the parameter number.
   */
  template <class IndexRange = std::initializer_list<size_type>,
            class = std::enable_if_t<detail::RangeTraits<IndexRange>::has_begin>>
  reference operator[](const IndexRange &indices) {
    GUDHI_CHECK(indices.size() >= 2,
                std::invalid_argument(
                    "Exactly 2 indices allowed only: first the generator number, second the parameter number."));
    auto it = indices.begin();
    size_type g = *it;
    return this->operator()(g, *(++it));
  }

  /**
   * @brief Let \f$ g \f$ be the first value in `indices` and \f$ p \f$ the second value.
   * Returns const reference to or copy of value of parameter \f$ p \f$ of generator \f$ g \f$.
   *
   * @tparam IndexRange Range with a begin() and size() method.
   * @param indices Range with at least two elements. The first element should correspond to the generator number and
   * the second element to the parameter number.
   */
  template <class IndexRange = std::initializer_list<size_type>,
            class = std::enable_if_t<detail::RangeTraits<IndexRange>::has_begin>>
  const_reference operator[](const IndexRange &indices) const {
    GUDHI_CHECK(indices.size() >= 2,
                std::invalid_argument(
                    "Exactly 2 indices allowed only: first the generator number, second the parameter number."));
    auto it = indices.begin();
    size_type g = *it;
    return this->operator()(g, *(++it));
  }

  /**
   * @brief Returns begin iterator to generator `g`. See @ref iterator.
   */
  iterator begin(size_type g) { return generators_.begin(g); }

  /**
   * @brief Returns end iterator to generator `g`. See @ref iterator.
   */
  iterator end(size_type g) { return generators_.end(g); }

  /**
   * @brief Returns begin const iterator to generator `g`. See @ref const_iterator.
   */
  const_iterator begin(size_type g) const { return generators_.begin(g); }

  /**
   * @brief Returns end const iterator to generator `g`. See @ref const_iterator.
   */
  const_iterator end(size_type g) const { return generators_.end(g); }

  /**
   * @brief Returns the number of parameters in the filtration value.
   */
  size_type num_parameters() const { return generators_.num_parameters(); }

  /**
   * @brief Returns the number of generators in the filtration value, i.e. the criticality of the element.
   */
  constexpr size_type num_generators() const {
    if constexpr (Ensure1Criticality) {
      return 1;  // for possible optimizations
    } else {
      return generators_.num_generators();
    }
  }

  /**
   * @brief Returns the total number of values in the filtration value, that is,
   * @ref num_parameters() * @ref num_generators().
   */
  size_type num_entries() const { return generators_.num_entries(); }

  /**
   * @brief Returns a filtration value with given number of parameters for which @ref is_plus_inf() returns `true`.
   * @throw std::logic_error If @ref has_infinity is false.
   */
  static Multi_parameter_filtration_value inf(int numberOfParameters) {
    StoragePolicy out = StoragePolicy::template inf<Co>(numberOfParameters < 0 ? 0 : numberOfParameters);
    return Multi_parameter_filtration_value(std::move(out));
  }

  /**
   * @brief Returns a filtration value with given number of parameters for which @ref is_minus_inf() returns `true`.
   * @throw std::logic_error If @ref has_minus_infinity is false.
   */
  static Multi_parameter_filtration_value minus_inf(int numberOfParameters) {
    StoragePolicy out = StoragePolicy::template minus_inf<Co>(numberOfParameters < 0 ? 0 : numberOfParameters);
    return Multi_parameter_filtration_value(std::move(out));
  }

  /**
   * @brief Returns a filtration value with given number of parameters for which @ref is_nan() returns `true`.
   * @throw std::logic_error If @ref has_quiet_NaN is false.
   */
  static Multi_parameter_filtration_value nan(int numberOfParameters) {
    StoragePolicy out = StoragePolicy::nan(numberOfParameters < 0 ? 0 : numberOfParameters);
    return Multi_parameter_filtration_value(std::move(out));
  }

  // DESCRIPTORS

  /**
   * @brief Returns value of `Ensure1Criticality`.
   */
  static constexpr bool ensures_1_criticality() { return Ensure1Criticality; }

  /**
   * @brief Returns value of `Co`.
   */
  static constexpr bool has_negative_cones() { return Co; }

  /**
   * @brief Returns `true` if and only if the filtration value is considered as plus infinity.
   */
  [[nodiscard]] bool is_plus_inf() const { return generators_.template is_plus_inf<Co>(); }

  /**
   * @brief Returns `true` if and only if the filtration value is considered as minus infinity.
   */
  [[nodiscard]] bool is_minus_inf() const { return generators_.template is_minus_inf<Co>(); }

  /**
   * @brief Returns `true` if and only if the filtration value is considered as NaN.
   */
  [[nodiscard]] bool is_nan() const { return generators_.is_nan(); }

  /**
   * @brief Returns `true` if and only if the filtration value is not considered as plus infinity,
   * minus infinity or NaN.
   */
  [[nodiscard]] bool is_finite() const { return generators_.template is_finite<Co>(); }

  // COMPARISON OPERATORS

  /**
   * @brief Returns `true` if and only if the first argument is lexicographically strictly less than the second
   * argument. The "words" considered for the lexicographical order are all the generators concatenated together
   * in order of generator index and then in order of parameter index. Different from @ref operator< "", this order
   * is total. To make the order total, a NaN value will be considered higher than any other value possible.
   *
   * @tparam inverse If true, the parameter index and generator index order is inverted.
   */
  template <bool inverse = false>
  friend bool is_strict_less_than_lexicographically(const Multi_parameter_filtration_value &a,
                                                    const Multi_parameter_filtration_value &b) {
    if (a.is_nan()) return false;
    if (b.is_nan()) return true;
    if (&a == &b) return false;

    GUDHI_CHECK(a.num_parameters() == b.num_parameters(),
                std::invalid_argument("Only filtration values with same number of parameters can be compared."));

    if (a.is_minus_inf()) return !b.is_minus_inf();
    if (b.is_plus_inf()) return !a.is_plus_inf();
    if (b.is_minus_inf() || a.is_plus_inf()) return false;

    for (size_type g = 0U; g < std::min(a.num_generators(), b.num_generators()); ++g) {
      size_type gA = g;
      size_type gB = g;
      if constexpr (inverse) {
        gA = a.num_generators() - 1 - g;
        gB = b.num_generators() - 1 - g;
      }
      for (size_type p = 0U; p < a.num_parameters(); ++p) {
        if constexpr (inverse) p = a.num_parameters() - 1 - p;
        if (detail::_is_nan(a(gA, p)) && !detail::_is_nan(b(gB, p))) return false;
        if (detail::_is_nan(b(gB, p))) return true;
        if (a(gA, p) < b(gB, p)) return true;
        if (b(gB, p) < a(gA, p)) return false;
        if constexpr (inverse) p = a.num_parameters() - 1 - p;
      }
    }
    return a.num_generators() < b.num_generators();
  }

  /**
   * @brief Returns `true` if and only if the first argument is lexicographically less than or equal to the second
   * argument. The "words" considered for the lexicographical order are all the generators concatenated together
   * in order of generator index and then in order of parameter index. Different from @ref operator<= "", this order
   * is total. To make the order total, a NaN value will be considered higher than any other value possible.
   *
   * @tparam inverse If true, the parameter index and generator index order is inverted.
   */
  template <bool inverse = false>
  friend bool is_less_or_equal_than_lexicographically(const Multi_parameter_filtration_value &a,
                                                      const Multi_parameter_filtration_value &b) {
    if (b.is_nan()) return true;
    if (a.is_nan()) return false;
    if (&a == &b) return true;

    GUDHI_CHECK(a.num_parameters() == b.num_parameters(),
                std::invalid_argument("Only filtration values with same number of parameters can be compared."));

    if (a.is_plus_inf()) return b.is_plus_inf();
    if (b.is_minus_inf()) return a.is_minus_inf();
    if (b.is_plus_inf() || a.is_minus_inf()) return true;

    for (size_type g = 0U; g < std::min(a.num_generators(), b.num_generators()); ++g) {
      size_type gA = g;
      size_type gB = g;
      if constexpr (inverse) {
        gA = a.num_generators() - 1 - g;
        gB = b.num_generators() - 1 - g;
      }
      for (size_type p = 0U; p < a.num_parameters(); ++p) {
        if constexpr (inverse) p = a.num_parameters() - 1 - p;
        if (detail::_is_nan(a(gA, p)) && !detail::_is_nan(b(gB, p))) return false;
        if (detail::_is_nan(b(gB, p))) return true;
        if (a(gA, p) < b(gB, p)) return true;
        if (b(gB, p) < a(gA, p)) return false;
        if constexpr (inverse) p = a.num_parameters() - 1 - p;
      }
    }
    return a.num_generators() <= b.num_generators();
  }

  /**
   * @brief Returns `true` if and only if the cones generated by @p b are strictly contained in the
   * cones generated by @p a (recall that the cones are positive if `Co` is false and negative if `Co` is true).
   * Both @p a and @p b must have to have the same number of parameters.
   *
   * Note that not all filtration values are comparable. That is, \f$ a < b \f$ and \f$ b < a \f$ returning both false
   * does **not** imply \f$ a == b \f$. If a total order is needed, use @ref is_strict_less_than_lexicographically
   * instead.
   */
  friend bool operator<(const Multi_parameter_filtration_value &a, const Multi_parameter_filtration_value &b) {
    return _is_lower_than<true>(a, b);
  }

  /**
   * @brief Returns `true` if and only if the cones generated by @p b are contained in or are (partially)
   * equal to the cones generated by @p a (recall that the cones are positive if `Co` is false and negative if `Co` is
   * true).
   * Both @p a and @p b must have to have the same number of parameters.
   *
   * Note that not all filtration values are comparable. That is, \f$ a \le b \f$ and \f$ b \le a \f$ can both return
   * `false`. If a total order is needed, use @ref is_less_or_equal_than_lexicographically instead.
   */
  friend bool operator<=(const Multi_parameter_filtration_value &a, const Multi_parameter_filtration_value &b) {
    return _is_lower_than<false>(a, b);
  }

  /**
   * @brief Returns `true` if and only if the cones generated by @p a are strictly contained in the
   * cones generated by @p b (recall that the cones are positive if `Co` is false and negative if `Co` is true).
   * Both @p a and @p b must have to have the same number of parameters.
   *
   * Note that not all filtration values are comparable. That is, \f$ a > b \f$ and \f$ b > a \f$ returning both false
   * does **not** imply \f$ a == b \f$. If a total order is needed, use @ref is_strict_less_than_lexicographically
   * instead.
   */
  friend bool operator>(const Multi_parameter_filtration_value &a, const Multi_parameter_filtration_value &b) {
    return b < a;
  }

  /**
   * @brief Returns `true` if and only if the cones generated by @p a are contained in or are (partially)
   * equal to the cones generated by @p b (recall that the cones are positive if `Co` is false and negative if `Co` is
   * true).
   * Both @p a and @p b must have to have the same number of parameters.
   *
   * Note that not all filtration values are comparable. That is, \f$ a \ge b \f$ and \f$ b \ge a \f$ can both return
   * `false`. If a total order is needed, use @ref is_less_or_equal_than_lexicographically instead.
   */
  friend bool operator>=(const Multi_parameter_filtration_value &a, const Multi_parameter_filtration_value &b) {
    return b <= a;
  }

  /**
   * @brief Returns `true` if and only if for each \f$ i,j \f$, \f$ a(i,j) \f$ is equal to \f$ b(i,j) \f$.
   * In particular, if both do not have the same number of generators or parameters or one of them is NaN,
   * returns `false`.
   */
  friend bool operator==(const Multi_parameter_filtration_value &a, const Multi_parameter_filtration_value &b) {
    if (a.is_nan() || b.is_nan()) return false;
    if (&a == &b) return true;
    if (a.num_generators() != b.num_generators()) return false;
    if (a.num_parameters() != b.num_parameters()) return false;
    // assumes that for same StoragePolicy, the internal order is the same,
    // i.e. if both are equal, the generators and parameters are ordered the same
    for (size_type g = 0; g < a.num_generators(); ++g) {
      for (size_type p = 0; p < a.num_parameters(); ++p) {
        value_type va = a(g, p);
        value_type vb = b(g, p);
        if (detail::_is_nan(va) || detail::_is_nan(vb)) return false;
        if (va != vb) return false;
      }
    }
    return true;
  }

  /**
   * @brief Returns `true` if and only if \f$ a == b \f$ returns `false`.
   */
  friend bool operator!=(const Multi_parameter_filtration_value &a, const Multi_parameter_filtration_value &b) {
    return !(a == b);
  }

  // ARITHMETIC OPERATORS

  // opposite
  /**
   * @brief Returns a filtration value such that an entry at index \f$ i,j \f$ is equal to \f$ -f(i,j) \f$.
   * If `StoragePolicy` has a fixed parameter, this one will not change value.
   *
   * Used conventions:
   * - \f$ -NaN = NaN \f$.
   *
   * @param f Value to opposite.
   * @return The opposite of @p f.
   */
  friend Multi_parameter_filtration_value operator-(const Multi_parameter_filtration_value &f) {
    using F = Multi_parameter_filtration_value;

    if (f.is_nan()) return f;

    Multi_parameter_filtration_value out = f;

    for (size_type g = 0; g < f.num_generators(); ++g) {
      for (size_type p = 0; p < f.num_parameters(); ++p) {
        if constexpr (StoragePolicy::has_an_implicit_axis) {
          if (StoragePolicy::implicit_axis() == p) continue;
        }
        auto &v = out(g, p);
        if (v == F::T_inf)
          v = F::T_m_inf;
        else if (v == F::T_m_inf)
          v = F::T_inf;
        else
          v = -v;
      }
    }

    return out;
  }

  // subtraction
  /**
   * @brief Returns a filtration value such that an entry at index \f$ (g,p) \f$, with
   * \f$ 0 \leq g \leq num_generators \f$ and \f$ 0 \leq p \leq num_parameters \f$ is equal to \f$ f(g,p) - r(p) \f$
   * if \f$ p < length_r \f$ and to \f$ f(g,p) \f$ otherwise.
   * If `StoragePolicy` has a fixed parameter, this one will not change value.
   *
   * Used conventions:
   * - \f$ inf - inf = NaN \f$,
   * - \f$ -inf - (-inf) = NaN \f$,
   * - \f$ NaN - b = NaN \f$,
   * - \f$ a - NaN = NaN \f$.
   *
   * All NaN values are represented by `std::numeric_limits<T>::quiet_NaN()` independently if
   * `std::numeric_limits<T>::has_quiet_NaN` is true or not.
   *
   * @tparam ValueRange Range with a begin() and end() method.
   * @param f First element of the subtraction.
   * @param r Second element of the subtraction.
   */
  template <class ValueRange, class = std::enable_if_t<detail::RangeTraits<ValueRange>::has_begin ||
                                                       detail::RangeTraits<ValueRange>::is_multi_filtration>>
  friend Multi_parameter_filtration_value operator-(Multi_parameter_filtration_value f, const ValueRange &r) {
    f -= r;
    return f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ (g,p) \f$, with
   * \f$ 0 \leq g \leq num_generators \f$ and \f$ 0 \leq p \leq num_parameters \f$ is equal to \f$ r(p) - f(g,p) \f$
   * if \f$ p < length_r \f$ and to \f$ -f(g,p) \f$ otherwise.
   * If `StoragePolicy` has a fixed parameter, this one will not change value.
   *
   * Used conventions:
   * - \f$ inf - inf = NaN \f$,
   * - \f$ -inf - (-inf) = NaN \f$,
   * - \f$ NaN - b = NaN \f$,
   * - \f$ a - NaN = NaN \f$.
   *
   * All NaN values are represented by `std::numeric_limits<T>::quiet_NaN()` independently if
   * `std::numeric_limits<T>::has_quiet_NaN` is true or not.
   *
   * @tparam ValueRange Range with a begin() and end() method.
   * @param r First element of the subtraction.
   * @param f Second element of the subtraction.
   */
  template <class ValueRange, class = std::enable_if_t<detail::RangeTraits<ValueRange>::has_begin &&
                                                       !detail::RangeTraits<ValueRange>::is_multi_filtration>>
  friend Multi_parameter_filtration_value operator-(const ValueRange &r, Multi_parameter_filtration_value f) {
    f._apply_operation(r, [](value_type &valF, const value_type &valR) -> bool {
      valF = -valF;
      return detail::_add(valF, valR);
    });
    return f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ (g,p) \f$ is equal to \f$ f(g,p) - val \f$.
   * If `StoragePolicy` has a fixed parameter, this one will not change value.
   *
   * Used conventions:
   * - \f$ inf - inf = NaN \f$,
   * - \f$ -inf - (-inf) = NaN \f$,
   * - \f$ NaN - b = NaN \f$,
   * - \f$ a - NaN = NaN \f$.
   *
   * All NaN values are represented by `std::numeric_limits<T>::quiet_NaN()` independently if
   * `std::numeric_limits<T>::has_quiet_NaN` is true or not.
   *
   * @param f First element of the subtraction.
   * @param val Second element of the subtraction.
   */
  friend Multi_parameter_filtration_value operator-(Multi_parameter_filtration_value f, const value_type &val) {
    f -= val;
    return f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ (g,p) \f$ is equal to \f$ val - f(g,p) \f$.
   * If `StoragePolicy` has a fixed parameter, this one will not change value.
   *
   * Used conventions:
   * - \f$ inf - inf = NaN \f$,
   * - \f$ -inf - (-inf) = NaN \f$,
   * - \f$ NaN - b = NaN \f$,
   * - \f$ a - NaN = NaN \f$.
   *
   * All NaN values are represented by `std::numeric_limits<T>::quiet_NaN()` independently if
   * `std::numeric_limits<T>::has_quiet_NaN` is true or not.
   *
   * @param val First element of the subtraction.
   * @param f Second element of the subtraction.
   */
  friend Multi_parameter_filtration_value operator-(const value_type &val, Multi_parameter_filtration_value f) {
    f._apply_operation([val](value_type &valF) -> bool {
      valF = -valF;
      return detail::_add(valF, val);
    });
    return f;
  }

  /**
   * @brief Modifies the first parameter such that an entry at index \f$ (g,p) \f$, with
   * \f$ 0 \leq g \leq num_generators \f$ and \f$ 0 \leq p \leq num_parameters \f$ is equal to \f$ f(g,p) - r(p) \f$
   * if \f$ p < length_r \f$ and to \f$ f(g,p) \f$ otherwise.
   * If `StoragePolicy` has a fixed parameter, this one will not change value.
   *
   * Used conventions:
   * - \f$ inf - inf = NaN \f$,
   * - \f$ -inf - (-inf) = NaN \f$,
   * - \f$ NaN - b = NaN \f$,
   * - \f$ a - NaN = NaN \f$.
   *
   * All NaN values are represented by `std::numeric_limits<T>::quiet_NaN()` independently if
   * `std::numeric_limits<T>::has_quiet_NaN` is true or not.
   *
   * @tparam ValueRange Range with a begin() and end() method.
   * @param f First element of the subtraction.
   * @param r Second element of the subtraction.
   */
  template <class ValueRange, class = std::enable_if_t<detail::RangeTraits<ValueRange>::has_begin ||
                                                       detail::RangeTraits<ValueRange>::is_multi_filtration>>
  friend Multi_parameter_filtration_value &operator-=(Multi_parameter_filtration_value &f, const ValueRange &r) {
    f._apply_operation(r,
                       [](value_type &valF, const value_type &valR) -> bool { return detail::_subtract(valF, valR); });
    return f;
  }

  /**
   * @brief Modifies the first parameter such that an entry at index \f$ (g,p) \f$ is equal to \f$ f(g,p) - val \f$.
   * If `StoragePolicy` has a fixed parameter, this one will not change value.
   *
   * Used conventions:
   * - \f$ inf - inf = NaN \f$,
   * - \f$ -inf - (-inf) = NaN \f$,
   * - \f$ NaN - b = NaN \f$,
   * - \f$ a - NaN = NaN \f$.
   *
   * All NaN values are represented by `std::numeric_limits<T>::quiet_NaN()` independently if
   * `std::numeric_limits<T>::has_quiet_NaN` is true or not.
   *
   * @param f First element of the subtraction.
   * @param val Second element of the subtraction.
   */
  friend Multi_parameter_filtration_value &operator-=(Multi_parameter_filtration_value &f, const value_type &val) {
    f._apply_operation([val](value_type &valF) -> bool { return detail::_subtract(valF, val); });
    return f;
  }

  // addition
  /**
   * @brief Returns a filtration value such that an entry at index \f$ (g,p) \f$, with
   * \f$ 0 \leq g \leq num_generators \f$ and \f$ 0 \leq p \leq num_parameters \f$ is equal to \f$ f(g,p) + r(p) \f$
   * if \f$ p < length_r \f$ and to \f$ f(g,p) \f$ otherwise.
   * If `StoragePolicy` has a fixed parameter, this one will not change value.
   *
   * Used conventions:
   * - \f$ inf + (-inf) = NaN \f$,
   * - \f$ -inf + inf = NaN \f$,
   * - \f$ NaN + b = NaN \f$,
   * - \f$ a + NaN = NaN \f$.
   *
   * All NaN values are represented by `std::numeric_limits<T>::quiet_NaN()` independently if
   * `std::numeric_limits<T>::has_quiet_NaN` is true or not.
   *
   * @tparam ValueRange Range with a begin() and end() method.
   * @param f First element of the addition.
   * @param r Second element of the addition.
   */
  template <class ValueRange, class = std::enable_if_t<detail::RangeTraits<ValueRange>::has_begin ||
                                                       detail::RangeTraits<ValueRange>::is_multi_filtration>>
  friend Multi_parameter_filtration_value operator+(Multi_parameter_filtration_value f, const ValueRange &r) {
    f += r;
    return f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ (g,p) \f$, with
   * \f$ 0 \leq g \leq num_generators \f$ and \f$ 0 \leq p \leq num_parameters \f$ is equal to \f$ r(p) + f(g,p) \f$
   * if \f$ p < length_r \f$ and to \f$ f(g,p) \f$ otherwise.
   * If `StoragePolicy` has a fixed parameter, this one will not change value.
   *
   * Used conventions:
   * - \f$ inf + (-inf) = NaN \f$,
   * - \f$ -inf + inf = NaN \f$,
   * - \f$ NaN + b = NaN \f$,
   * - \f$ a + NaN = NaN \f$.
   *
   * All NaN values are represented by `std::numeric_limits<T>::quiet_NaN()` independently if
   * `std::numeric_limits<T>::has_quiet_NaN` is true or not.
   *
   * @tparam ValueRange Range with a begin() and end() method.
   * @param r First element of the addition.
   * @param f Second element of the addition.
   */
  template <class ValueRange, class = std::enable_if_t<detail::RangeTraits<ValueRange>::has_begin &&
                                                       !detail::RangeTraits<ValueRange>::is_multi_filtration>>
  friend Multi_parameter_filtration_value operator+(const ValueRange &r, Multi_parameter_filtration_value f) {
    f += r;
    return f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ (g,p) \f$ is equal to \f$ f(g,p) + val \f$.
   * If `StoragePolicy` has a fixed parameter, this one will not change value.
   *
   * Used conventions:
   * - \f$ inf + (-inf) = NaN \f$,
   * - \f$ -inf + inf = NaN \f$,
   * - \f$ NaN + b = NaN \f$,
   * - \f$ a + NaN = NaN \f$.
   *
   * All NaN values are represented by `std::numeric_limits<T>::quiet_NaN()` independently if
   * `std::numeric_limits<T>::has_quiet_NaN` is true or not.
   *
   * @param f First element of the addition.
   * @param val Second element of the addition.
   */
  friend Multi_parameter_filtration_value operator+(Multi_parameter_filtration_value f, const value_type &val) {
    f += val;
    return f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ (g,p) \f$ is equal to \f$ val + f(g,p) \f$.
   * If `StoragePolicy` has a fixed parameter, this one will not change value.
   *
   * Used conventions:
   * - \f$ inf + (-inf) = NaN \f$,
   * - \f$ -inf + inf = NaN \f$,
   * - \f$ NaN + b = NaN \f$,
   * - \f$ a + NaN = NaN \f$.
   *
   * All NaN values are represented by `std::numeric_limits<T>::quiet_NaN()` independently if
   * `std::numeric_limits<T>::has_quiet_NaN` is true or not.
   *
   * @param val First element of the addition.
   * @param f Second element of the addition.
   */
  friend Multi_parameter_filtration_value operator+(const value_type &val, Multi_parameter_filtration_value f) {
    f += val;
    return f;
  }

  /**
   * @brief Modifies the first parameter such that an entry at index \f$ (g,p) \f$, with
   * \f$ 0 \leq g \leq num_generators \f$ and \f$ 0 \leq p \leq num_parameters \f$ is equal to \f$ f(g,p) + r(p) \f$
   * if \f$ p < length_r \f$ and to \f$ f(g,p) \f$ otherwise.
   * If `StoragePolicy` has a fixed parameter, this one will not change value.
   *
   * Used conventions:
   * - \f$ inf + (-inf) = NaN \f$,
   * - \f$ -inf + inf = NaN \f$,
   * - \f$ NaN + b = NaN \f$,
   * - \f$ a + NaN = NaN \f$.
   *
   * All NaN values are represented by `std::numeric_limits<T>::quiet_NaN()` independently if
   * `std::numeric_limits<T>::has_quiet_NaN` is true or not.
   *
   * @tparam ValueRange Range with a begin() and end() method.
   * @param f First element of the addition.
   * @param r Second element of the addition.
   */
  template <class ValueRange, class = std::enable_if_t<detail::RangeTraits<ValueRange>::has_begin ||
                                                       detail::RangeTraits<ValueRange>::is_multi_filtration>>
  friend Multi_parameter_filtration_value &operator+=(Multi_parameter_filtration_value &f, const ValueRange &r) {
    f._apply_operation(r, [](value_type &valF, const value_type &valR) -> bool { return detail::_add(valF, valR); });
    return f;
  }

  /**
   * @brief Modifies the first parameter such that an entry at index \f$ (g,p) \f$ is equal to \f$ f(g,p) + val \f$.
   * If `StoragePolicy` has a fixed parameter, this one will not change value.
   *
   * Used conventions:
   * - \f$ inf + (-inf) = NaN \f$,
   * - \f$ -inf + inf = NaN \f$,
   * - \f$ NaN + b = NaN \f$,
   * - \f$ a + NaN = NaN \f$.
   *
   * All NaN values are represented by `std::numeric_limits<T>::quiet_NaN()` independently if
   * `std::numeric_limits<T>::has_quiet_NaN` is true or not.
   *
   * @param f First element of the addition.
   * @param val Second element of the addition.
   */
  friend Multi_parameter_filtration_value &operator+=(Multi_parameter_filtration_value &f, const value_type &val) {
    f._apply_operation([val](value_type &valF) -> bool { return detail::_add(valF, val); });
    return f;
  }

  // multiplication
  /**
   * @brief Returns a filtration value such that an entry at index \f$ (g,p) \f$, with
   * \f$ 0 \leq g \leq num_generators \f$ and \f$ 0 \leq p \leq num_parameters \f$ is equal to \f$ f(g,p) * r(p) \f$
   * if \f$ p < length_r \f$ and to \f$ f(g,p) \f$ otherwise.
   * If `StoragePolicy` has a fixed parameter, this one will not change value.
   *
   * Used conventions:
   * - \f$ inf * 0 = NaN \f$,
   * - \f$ 0 * inf = NaN \f$,
   * - \f$ -inf * 0 = NaN \f$,
   * - \f$ 0 * (-inf) = NaN \f$,
   * - \f$ NaN * b = NaN \f$,
   * - \f$ a * NaN = NaN \f$.
   *
   * All NaN values are represented by `std::numeric_limits<T>::quiet_NaN()` independently if
   * `std::numeric_limits<T>::has_quiet_NaN` is true or not.
   *
   * @tparam ValueRange Range with a begin() and end() method.
   * @param f First element of the multiplication.
   * @param r Second element of the multiplication.
   */
  template <class ValueRange, class = std::enable_if_t<detail::RangeTraits<ValueRange>::has_begin ||
                                                       detail::RangeTraits<ValueRange>::is_multi_filtration>>
  friend Multi_parameter_filtration_value operator*(Multi_parameter_filtration_value f, const ValueRange &r) {
    f *= r;
    return f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ (g,p) \f$, with
   * \f$ 0 \leq g \leq num_generators \f$ and \f$ 0 \leq p \leq num_parameters \f$ is equal to \f$ r(p) * f(g,p) \f$
   * if \f$ p < length_r \f$ and to \f$ f(g,p) \f$ otherwise.
   * If `StoragePolicy` has a fixed parameter, this one will not change value.
   *
   * Used conventions:
   * - \f$ inf * 0 = NaN \f$,
   * - \f$ 0 * inf = NaN \f$,
   * - \f$ -inf * 0 = NaN \f$,
   * - \f$ 0 * (-inf) = NaN \f$,
   * - \f$ NaN * b = NaN \f$,
   * - \f$ a * NaN = NaN \f$.
   *
   * All NaN values are represented by `std::numeric_limits<T>::quiet_NaN()` independently if
   * `std::numeric_limits<T>::has_quiet_NaN` is true or not.
   *
   * @tparam ValueRange Range with a begin() and end() method.
   * @param r First element of the multiplication.
   * @param f Second element of the multiplication.
   */
  template <class ValueRange, class = std::enable_if_t<detail::RangeTraits<ValueRange>::has_begin &&
                                                       !detail::RangeTraits<ValueRange>::is_multi_filtration>>
  friend Multi_parameter_filtration_value operator*(const ValueRange &r, Multi_parameter_filtration_value f) {
    f *= r;
    return f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ (g,p) \f$ is equal to \f$ f(g,p) * val \f$.
   * If `StoragePolicy` has a fixed parameter, this one will not change value.
   *
   * Used conventions:
   * - \f$ inf * 0 = NaN \f$,
   * - \f$ 0 * inf = NaN \f$,
   * - \f$ -inf * 0 = NaN \f$,
   * - \f$ 0 * (-inf) = NaN \f$,
   * - \f$ NaN * b = NaN \f$,
   * - \f$ a * NaN = NaN \f$.
   *
   * All NaN values are represented by `std::numeric_limits<T>::quiet_NaN()` independently if
   * `std::numeric_limits<T>::has_quiet_NaN` is true or not.
   *
   * @param f First element of the multiplication.
   * @param val Second element of the multiplication.
   */
  friend Multi_parameter_filtration_value operator*(Multi_parameter_filtration_value f, const value_type &val) {
    f *= val;
    return f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ (g,p) \f$ is equal to \f$ val * f(g,p) \f$.
   * If `StoragePolicy` has a fixed parameter, this one will not change value.
   *
   * Used conventions:
   * - \f$ inf * 0 = NaN \f$,
   * - \f$ 0 * inf = NaN \f$,
   * - \f$ -inf * 0 = NaN \f$,
   * - \f$ 0 * (-inf) = NaN \f$,
   * - \f$ NaN * b = NaN \f$,
   * - \f$ a * NaN = NaN \f$.
   *
   * All NaN values are represented by `std::numeric_limits<T>::quiet_NaN()` independently if
   * `std::numeric_limits<T>::has_quiet_NaN` is true or not.
   *
   * @param val First element of the multiplication.
   * @param f Second element of the multiplication.
   */
  friend Multi_parameter_filtration_value operator*(const value_type &val, Multi_parameter_filtration_value f) {
    f *= val;
    return f;
  }

  /**
   * @brief Modifies the first parameter such that an entry at index \f$ (g,p) \f$, with
   * \f$ 0 \leq g \leq num_generators \f$ and \f$ 0 \leq p \leq num_parameters \f$ is equal to \f$ f(g,p) * r(p) \f$
   * if \f$ p < length_r \f$ and to \f$ f(g,p) \f$ otherwise.
   * If `StoragePolicy` has a fixed parameter, this one will not change value.
   *
   * Used conventions:
   * - \f$ inf * 0 = NaN \f$,
   * - \f$ 0 * inf = NaN \f$,
   * - \f$ -inf * 0 = NaN \f$,
   * - \f$ 0 * (-inf) = NaN \f$,
   * - \f$ NaN * b = NaN \f$,
   * - \f$ a * NaN = NaN \f$.
   *
   * All NaN values are represented by `std::numeric_limits<T>::quiet_NaN()` independently if
   * `std::numeric_limits<T>::has_quiet_NaN` is true or not.
   *
   * @tparam ValueRange Range with a begin() and end() method.
   * @param f First element of the multiplication.
   * @param r Second element of the multiplication.
   */
  template <class ValueRange, class = std::enable_if_t<detail::RangeTraits<ValueRange>::has_begin ||
                                                       detail::RangeTraits<ValueRange>::is_multi_filtration>>
  friend Multi_parameter_filtration_value &operator*=(Multi_parameter_filtration_value &f, const ValueRange &r) {
    f._apply_operation(r,
                       [](value_type &valF, const value_type &valR) -> bool { return detail::_multiply(valF, valR); });
    return f;
  }

  /**
   * @brief Modifies the first parameter such that an entry at index \f$ (g,p) \f$ is equal to \f$ f(g,p) * val \f$.
   * If `StoragePolicy` has a fixed parameter, this one will not change value.
   *
   * Used conventions:
   * - \f$ inf * 0 = NaN \f$,
   * - \f$ 0 * inf = NaN \f$,
   * - \f$ -inf * 0 = NaN \f$,
   * - \f$ 0 * (-inf) = NaN \f$,
   * - \f$ NaN * b = NaN \f$,
   * - \f$ a * NaN = NaN \f$.
   *
   * All NaN values are represented by `std::numeric_limits<T>::quiet_NaN()` independently if
   * `std::numeric_limits<T>::has_quiet_NaN` is true or not.
   *
   * @param f First element of the multiplication.
   * @param val Second element of the multiplication.
   */
  friend Multi_parameter_filtration_value &operator*=(Multi_parameter_filtration_value &f, const value_type &val) {
    f._apply_operation([val](value_type &valF) -> bool { return detail::_multiply(valF, val); });
    return f;
  }

  // division
  /**
   * @brief Returns a filtration value such that an entry at index \f$ (g,p) \f$, with
   * \f$ 0 \leq g \leq num_generators \f$ and \f$ 0 \leq p \leq num_parameters \f$ is equal to \f$ f(g,p) / r(p) \f$
   * if \f$ p < length_r \f$ and to \f$ f(g,p) \f$ otherwise.
   * If `StoragePolicy` has a fixed parameter, this one will not change value.
   *
   * Used conventions:
   * - \f$ a / 0 = NaN \f$,
   * - \f$ inf / inf = NaN \f$,
   * - \f$ -inf / inf = NaN \f$,
   * - \f$ inf / -inf = NaN \f$,
   * - \f$ -inf / -inf = NaN \f$,
   * - \f$ NaN / b = NaN \f$,
   * - \f$ a / NaN = NaN \f$,
   * - \f$ a / inf = 0 \f$,
   * - \f$ a / -inf = 0 \f$.
   *
   * All NaN values are represented by `std::numeric_limits<T>::quiet_NaN()` independently if
   * `std::numeric_limits<T>::has_quiet_NaN` is true or not.
   *
   * @tparam ValueRange Range with a begin() and end() method.
   * @param f First element of the division.
   * @param r Second element of the division.
   */
  template <class ValueRange, class = std::enable_if_t<detail::RangeTraits<ValueRange>::has_begin ||
                                                       detail::RangeTraits<ValueRange>::is_multi_filtration>>
  friend Multi_parameter_filtration_value operator/(Multi_parameter_filtration_value f, const ValueRange &r) {
    f /= r;
    return f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ (g,p) \f$, with
   * \f$ 0 \leq g \leq num_generators \f$ and \f$ 0 \leq p \leq num_parameters \f$ is equal to \f$ r(p) / f(g,p) \f$
   * if \f$ p < length_r \f$ and to \f$ f(g,p) \f$ otherwise.
   * If `StoragePolicy` has a fixed parameter, this one will not change value.
   *
   * Used conventions:
   * - \f$ a / 0 = NaN \f$,
   * - \f$ inf / inf = NaN \f$,
   * - \f$ -inf / inf = NaN \f$,
   * - \f$ inf / -inf = NaN \f$,
   * - \f$ -inf / -inf = NaN \f$,
   * - \f$ NaN / b = NaN \f$,
   * - \f$ a / NaN = NaN \f$,
   * - \f$ a / inf = 0 \f$,
   * - \f$ a / -inf = 0 \f$.
   *
   * All NaN values are represented by `std::numeric_limits<T>::quiet_NaN()` independently if
   * `std::numeric_limits<T>::has_quiet_NaN` is true or not.
   *
   * @tparam ValueRange Range with a begin() and end() method.
   * @param r First element of the division.
   * @param f Second element of the division.
   */
  template <class ValueRange, class = std::enable_if_t<detail::RangeTraits<ValueRange>::has_begin &&
                                                       !detail::RangeTraits<ValueRange>::is_multi_filtration>>
  friend Multi_parameter_filtration_value operator/(const ValueRange &r, Multi_parameter_filtration_value f) {
    f._apply_operation(r, [](value_type &valF, const value_type &valR) -> bool {
      value_type tmp = valF;
      valF = valR;
      return detail::_divide(valF, tmp);
    });
    return f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ (g,p) \f$ is equal to \f$ f(g,p) / val \f$.
   * If `StoragePolicy` has a fixed parameter, this one will not change value.
   *
   * Used conventions:
   * - \f$ a / 0 = NaN \f$,
   * - \f$ inf / inf = NaN \f$,
   * - \f$ -inf / inf = NaN \f$,
   * - \f$ inf / -inf = NaN \f$,
   * - \f$ -inf / -inf = NaN \f$,
   * - \f$ NaN / b = NaN \f$,
   * - \f$ a / NaN = NaN \f$,
   * - \f$ a / inf = 0 \f$,
   * - \f$ a / -inf = 0 \f$.
   *
   * All NaN values are represented by `std::numeric_limits<T>::quiet_NaN()` independently if
   * `std::numeric_limits<T>::has_quiet_NaN` is true or not.
   *
   * @param f First element of the division.
   * @param val Second element of the division.
   */
  friend Multi_parameter_filtration_value operator/(Multi_parameter_filtration_value f, const value_type &val) {
    f /= val;
    return f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ (g,p) \f$ is equal to \f$ val / f(g,p) \f$.
   * If `StoragePolicy` has a fixed parameter, this one will not change value.
   *
   * Used conventions:
   * - \f$ a / 0 = NaN \f$,
   * - \f$ inf / inf = NaN \f$,
   * - \f$ -inf / inf = NaN \f$,
   * - \f$ inf / -inf = NaN \f$,
   * - \f$ -inf / -inf = NaN \f$,
   * - \f$ NaN / b = NaN \f$,
   * - \f$ a / NaN = NaN \f$,
   * - \f$ a / inf = 0 \f$,
   * - \f$ a / -inf = 0 \f$.
   *
   * All NaN values are represented by `std::numeric_limits<T>::quiet_NaN()` independently if
   * `std::numeric_limits<T>::has_quiet_NaN` is true or not.
   *
   * @param val First element of the division.
   * @param f Second element of the division.
   */
  friend Multi_parameter_filtration_value operator/(const value_type &val, Multi_parameter_filtration_value f) {
    f._apply_operation([val](value_type &valF) -> bool {
      value_type tmp = valF;
      valF = val;
      return detail::_divide(valF, tmp);
    });
    return f;
  }

  /**
   * @brief Modifies the first parameter such that an entry at index \f$ (g,p) \f$, with
   * \f$ 0 \leq g \leq num_generators \f$ and \f$ 0 \leq p \leq num_parameters \f$ is equal to \f$ f(g,p) / r(p) \f$
   * if \f$ p < length_r \f$ and to \f$ f(g,p) \f$ otherwise.
   * If `StoragePolicy` has a fixed parameter, this one will not change value.
   *
   * Used conventions:
   * - \f$ a / 0 = NaN \f$,
   * - \f$ inf / inf = NaN \f$,
   * - \f$ -inf / inf = NaN \f$,
   * - \f$ inf / -inf = NaN \f$,
   * - \f$ -inf / -inf = NaN \f$,
   * - \f$ NaN / b = NaN \f$,
   * - \f$ a / NaN = NaN \f$,
   * - \f$ a / inf = 0 \f$,
   * - \f$ a / -inf = 0 \f$.
   *
   * All NaN values are represented by `std::numeric_limits<T>::quiet_NaN()` independently if
   * `std::numeric_limits<T>::has_quiet_NaN` is true or not.
   *
   * @tparam ValueRange Range with a begin() and end() method.
   * @param f First element of the division.
   * @param r Second element of the division.
   */
  template <class ValueRange, class = std::enable_if_t<detail::RangeTraits<ValueRange>::has_begin ||
                                                       detail::RangeTraits<ValueRange>::is_multi_filtration>>
  friend Multi_parameter_filtration_value &operator/=(Multi_parameter_filtration_value &f, const ValueRange &r) {
    f._apply_operation(r, [](value_type &valF, const value_type &valR) -> bool { return detail::_divide(valF, valR); });
    return f;
  }

  /**
   * @brief Modifies the first parameter such that an entry at index \f$ (g,p) \f$ is equal to \f$ f(g,p) / val \f$.
   * If `StoragePolicy` has a fixed parameter, this one will not change value.
   *
   * Used conventions:
   * - \f$ a / 0 = NaN \f$,
   * - \f$ inf / inf = NaN \f$,
   * - \f$ -inf / inf = NaN \f$,
   * - \f$ inf / -inf = NaN \f$,
   * - \f$ -inf / -inf = NaN \f$,
   * - \f$ NaN / b = NaN \f$,
   * - \f$ a / NaN = NaN \f$,
   * - \f$ a / inf = 0 \f$,
   * - \f$ a / -inf = 0 \f$.
   *
   * All NaN values are represented by `std::numeric_limits<T>::quiet_NaN()` independently if
   * `std::numeric_limits<T>::has_quiet_NaN` is true or not.
   *
   * @param f First element of the division.
   * @param val Second element of the division.
   */
  friend Multi_parameter_filtration_value &operator/=(Multi_parameter_filtration_value &f, const value_type &val) {
    f._apply_operation([val](value_type &valF) -> bool { return detail::_divide(valF, val); });
    return f;
  }

  // MODIFIERS

  /**
   * @brief Reserves space for the given number of generators in the underlying container. Does nothing if
   * `Ensure1Criticality` is true.
   */
  void reserve([[maybe_unused]] size_type number_of_generators) {
    if constexpr (!Ensure1Criticality) {
      generators_.reserve(number_of_generators);
    }
  }

  /**
   * @brief Sets the number of generators. If there were less generators before, new generators with default values are
   * constructed. If there were more generators before, the exceed of generators is destroyed (any generator with index
   * higher or equal to @p g to be more precise). If @p g is zero, the methods does nothing.
   *
   * Fails to compile if `Ensure1Criticality` is true.
   *
   * @warning All new generators will be set to -infinity (`Co` is true) or infinity (`Co` is false). That is, the new
   * filtration value is not minimal anymore. So, if @ref StoragePolicy::has_minimal_set_representation is true, make
   * sure to fill them with real generators or to remove them before using other methods.
   *
   * @warning Be sure to call @ref simplify if necessary after initializing all the generators.
   * If @ref StoragePolicy::has_minimal_set_representation is true, most methods will have
   * an undefined behaviour if the set of generators is not minimal or sorted.
   *
   * @param g New number of generators.
   */
  void set_num_generators(size_type g) {
    static_assert(!Ensure1Criticality, "Number of generators cannot be set for a 1-critical only filtration value.");

    if (g <= 0) return;
    generators_.set_num_generators(g, _get_default_null_value());
  }

  /**
   * @brief Adds the given generator to the filtration value such that the set remains sorted.
   * If @ref StoragePolicy::has_minimal_set_representation is true, it also ensures that it remains minimal. In that
   * case, it is therefore possible that the generator is ignored if it does not generated any new lifetime or that
   * old generators disappear if they are overshadowed by the new one. Note that this can also happen if
   * @ref StoragePolicy::has_minimal_set_representation is false: with true, it is just guaranteed.
   *
   * @tparam GeneratorRange Range of elements convertible to `T`. Must have a begin(), end() method and the iterator
   * type should satisfy the requirements of the standard `LegacyForwardIterator`.
   * @param x New generator to add. Has to have the same number of parameters than @ref num_parameters().
   * @return true If and only if the generator is actually added to the set of generators.
   * @return false Otherwise.
   */
  template <class GeneratorRange = std::initializer_list<value_type>,
            class = std::enable_if_t<detail::RangeTraits<GeneratorRange>::has_begin>>
  bool add_generator(const GeneratorRange &x) {
    return add_generator(x.begin(), x.end());
  }

  /**
   * @brief Adds the given generator to the filtration value such that the set remains sorted.
   * If @ref StoragePolicy::has_minimal_set_representation is true, it also ensures that it remains minimal. In that
   * case, it is therefore possible that the generator is ignored if it does not generated any new lifetime or that
   * old generators disappear if they are overshadowed by the new one. Note that this can also happen if
   * @ref StoragePolicy::has_minimal_set_representation is false: with true, it is just guaranteed.
   *
   * @tparam Iterator Iterator class satisfying the requirements of the standard `LegacyForwardIterator`.
   * The dereferenced type has to be convertible to `T`.
   * @param genStart Iterator pointing to the begining of the range.
   * @param genEnd Iterator pointing to the end of the range.
   * @return true If and only if the generator is actually added to the set of generators.
   * @return false Otherwise.
   */
  template <class Iterator>
  bool add_generator(Iterator genStart, Iterator genEnd) {
    GUDHI_CHECK(std::distance(genStart, genEnd) == static_cast<int>(num_parameters()),
                std::invalid_argument("Wrong range size. Should correspond to the number of parameters."));

    if constexpr (StoragePolicy::has_minimal_set_representation) {
      size_type end = num_generators();

      if (_generator_can_be_added(genStart, 0, end)) {
        generators_.set_num_generators(end, _get_default_null_value());  // never sizes up
        generators_.emplace_back(genStart, genEnd);
        if constexpr (Ensure1Criticality) {
          if (generators_.num_generators() != 1)
            throw std::logic_error("Multiparameter filtration value is not 1-critical anymore.");
        }
        generators_.sort();
        return true;
      }

      return false;
    } else {
      bool res = generators_.template emplace<Co>(genStart, genEnd, _get_default_null_value());
      if constexpr (Ensure1Criticality) {
        if (generators_.num_generators() != 1)
          throw std::logic_error("Multiparameter filtration value is not 1-critical anymore.");
      }
      if (res) generators_.sort();
      return res;
    }
  }

  /**
   * @brief Adds the given generator to the filtration value without any verifications or simplifications. It will be
   * added either at the end of the set or at the right position (w.r.t. the internal order) if `StoragePolicy` can do
   * so trivially.
   *
   * Fails to compile if `Ensure1Criticality` is true.
   *
   * @warning If @ref StoragePolicy::has_minimal_set_representation is true and the resulting set of generators is
   * not minimal or sorted after this modification, some methods will have an undefined behaviour. Be sure to
   * call @ref simplify() before using them. Same if the internal order is not preserved.
   *
   * @tparam GeneratorRange Range of elements convertible to `T`. Must have a begin(), end() and size() method.
   * @param x New generator to add. Must have the same number of parameters than @ref num_parameters().
   */
  template <class GeneratorRange = std::initializer_list<value_type>,
            class = std::enable_if_t<detail::RangeTraits<GeneratorRange>::has_begin>>
  void add_guaranteed_generator(const GeneratorRange &x) {
    static_assert(!Ensure1Criticality, "Cannot add additional generator to a 1-critical only filtration value.");
    add_guaranteed_generator(x.begin(), x.end());
  }

  /**
   * @brief Adds the given generator to the filtration value without any verifications or simplifications. It will be
   * added either at the end of the set or at the right position (w.r.t. the internal order) if `StoragePolicy` can do
   * so trivially.
   *
   * Fails to compile if `Ensure1Criticality` is true.
   *
   * @warning If @ref StoragePolicy::has_minimal_set_representation is true and the resulting set of generators is
   * not minimal or sorted after this modification, some methods will have an undefined behaviour. Be sure to
   * call @ref simplify() before using them. Same if the internal order is not preserved.
   *
   * @tparam Iterator Iterator class satisfying the requirements of the standard `LegacyForwardIterator`.
   * The dereferenced type has to be convertible to `T`.
   * @param genStart Iterator pointing to the begining of the range.
   * @param genEnd Iterator pointing to the end of the range.
   */
  template <class Iterator>
  void add_guaranteed_generator(Iterator genStart, Iterator genEnd) {
    static_assert(!Ensure1Criticality, "Cannot add additional generator to a 1-critical only filtration value.");

    GUDHI_CHECK(std::distance(genStart, genEnd) == static_cast<int>(num_parameters()),
                std::invalid_argument("Wrong range size. Should correspond to the number of parameters."));

    if constexpr (StoragePolicy::has_minimal_set_representation) {
      generators_.emplace_back(genStart, genEnd);
    } else {
      generators_.template emplace<Co>(genStart, genEnd, _get_default_null_value());
    }
  }

  /**
   * @brief If @ref StoragePolicy::has_minimal_set_representation is true, simplifies the current set of generators
   * such that it becomes minimal. Also orders it in increasing lexicographical order. Only necessary if generators
   * were added "by hand" without verification either through the constructor or with @ref add_guaranteed_generator "",
   * etc. If @ref StoragePolicy::has_minimal_set_representation is false, does nothing.
   */
  void simplify() {
    if constexpr (Ensure1Criticality || !StoragePolicy::has_minimal_set_representation) {
      return;
    } else {
      size_type end = 0;

      for (size_type curr = 0; curr < num_generators(); ++curr) {
        if (_generator_can_be_added(generators_.begin(curr), 0, end)) {
          generators_.swap_generators(end, curr);
          ++end;
        }
      }

      generators_.set_num_generators(end, _get_default_null_value());
      generators_.sort();
    }
  }

  /**
   * @brief Sets each generator of the filtration value to the least common upper bound between it and the given value.
   *
   * More formally, it pushes the current generator to the cone \f$ \{ y \in \mathbb R^n : y \ge x \} \f$
   * originating in \f$ x \f$. The resulting value corresponds to the intersection of both
   * cones: \f$ \mathrm{this} = \min \{ y \in \mathbb R^n : y \ge this \} \cap \{ y \in \mathbb R^n : y \ge x \} \f$.
   *
   * @tparam GeneratorRange Range of elements convertible to `T`. Must have a begin(), end() and size() method.
   * @param x Range towards to push. Has to have as many elements than @ref num_parameters().
   * @param exclude_infinite_values If true, values in x at infinity or minus infinity are ignored.
   * @return true If the filtration value was actually modified.
   * @return false Otherwise.
   */
  template <class GeneratorRange = std::initializer_list<value_type>,
            class = std::enable_if_t<detail::RangeTraits<GeneratorRange>::has_begin ||
                                     detail::RangeTraits<GeneratorRange>::is_multi_filtration>>
  bool push_to_least_common_upper_bound(const GeneratorRange &x, bool exclude_infinite_values = false) {
    std::vector<value_type> pushTo;
    if constexpr (detail::RangeTraits<GeneratorRange>::is_multi_filtration) {
      if (x.num_generators() == 0) return false;
      GUDHI_CHECK(x.num_parameters() == num_parameters(),
                  std::invalid_argument("Wrong range size. Should correspond to the number of parameters."));

      if (x.is_nan()) return false;
      if (x.is_plus_inf()) {
        if constexpr (StoragePolicy::template has_infinity<Co>) {
          if (is_plus_inf()) return false;
          generators_ = StoragePolicy::template inf<Co>(num_parameters());
        } else {
          throw std::invalid_argument("Cannot represent an infinite implicit axis with the current StoragePolicy.");
        }
        return true;
      }
      if (x.is_minus_inf()) return false;
      pushTo.resize(x.num_parameters());
      for (size_type p = 0; p < x.num_parameters(); ++p) pushTo[p] = x(0, p);
    } else {
      GUDHI_CHECK(x.size() == num_parameters(),
                  std::invalid_argument("Wrong range size. Should correspond to the number of parameters."));

      pushTo.insert(pushTo.end(), x.begin(), x.end());
    }

    if constexpr (StoragePolicy::has_an_implicit_axis) {
      if (is_nan()) return false;

      auto implicitAxisGen = generators_.get_generator_of_implicit_value(pushTo[StoragePolicy::implicit_axis()]);
      using I = decltype(implicitAxisGen);
      GUDHI_CHECK(implicitAxisGen != StoragePolicy::template null_value<I>(),
                  std::invalid_argument("Value at implicit axis is not valid."));

      if (implicitAxisGen == StoragePolicy::template T_inf<I>) {
        if constexpr (StoragePolicy::template has_infinity<Co>) {
          if (is_plus_inf()) return false;
          generators_ = StoragePolicy::template inf<Co>(num_parameters());
        } else {
          throw std::invalid_argument("Cannot represent an infinite implicit axis with the current StoragePolicy.");
        }
        return true;
      }

      if constexpr (Ensure1Criticality) {
        if (implicitAxisGen != 0)
          throw std::invalid_argument(
              "Pushing to an index other than 0 or +inf is not permitted with Ensure1Criticality at true.");
      }

      bool modified = false;
      auto implGen = static_cast<size_type>(implicitAxisGen);
      for (size_type p = 0; p < num_parameters(); ++p) {
        if (p != StoragePolicy::implicit_axis() &&
            (!exclude_infinite_values || (pushTo[p] != T_inf && pushTo[p] != T_m_inf)) && !detail::_is_nan(pushTo[p])) {
          value_type newVal = _get_default_null_value();

          for (size_type i = 0; i < std::min(implGen, generators_.num_generators()); ++i) {
            auto &v = generators_(i, p);
            if (!detail::_is_nan(v)) {
              value_type threshold = std::max(v, pushTo[p]);
              if (_dominates(newVal, threshold)) newVal = threshold;
              if (v != _get_default_null_value()) {
                modified = true;
                v = _get_default_null_value();
              }
            }
          }

          if (implGen < generators_.num_generators()) {
            auto &vs = generators_(implGen, p);
            value_type threshold = std::max(vs, pushTo[p]);
            if (_dominates(newVal, threshold)) newVal = threshold;
            if (newVal != vs) {
              modified = true;
              vs = newVal;
            }
            for (size_type i = implGen + 1; i < generators_.num_generators(); ++i) {
              auto &v = generators_(i, p);
              if (!detail::_is_nan(v) && v < pushTo[p]) {
                modified = true;
                v = pushTo[p];
              }
            }
          } else {
            pushTo[p] = newVal;
          }
        }
      }
      if (implGen >= generators_.num_generators()) {
        modified = true;
        generators_.template emplace<Co>(pushTo.begin(), pushTo.end(), _get_default_null_value());
      }

      // to potentially simplify
      if constexpr (StoragePolicy::template has_infinity<Co>) {
        if (is_plus_inf()) generators_ = StoragePolicy::template inf<Co>(num_parameters());
      }
      return modified;
    } else {
      bool xIsInf = true, xIsMinusInf = true, xIsNaN = true;
      bool thisIsInf = false, thisIsMinusInf = false, thisIsNaN = false;

      // if one is not finite, we can avoid the heavy simplification process
      _get_infinity_statuses(pushTo, xIsInf, xIsMinusInf, xIsNaN);
      _get_infinity_statuses(thisIsInf, thisIsMinusInf, thisIsNaN);

      if (thisIsInf || thisIsNaN || xIsNaN || xIsMinusInf || (xIsInf && exclude_infinite_values)) return false;

      if (xIsInf) {
        generators_ = StoragePolicy::template inf<Co>(num_parameters());
        return true;
      }

      if (thisIsMinusInf) {
        generators_ = StoragePolicy(pushTo.begin(), pushTo.end());
        return true;
      }

      bool modified = false;

      for (size_type g = 0; g < num_generators(); ++g) {
        for (size_type p = 0; p < num_parameters(); ++p) {
          value_type valX = pushTo[p];
          auto &val = generators_(g, p);
          if (!exclude_infinite_values || (valX != T_inf && valX != T_m_inf)) {
            modified |= val < valX;
            val = valX > val ? valX : val;
          }
        }
      }

      if (modified && num_generators() > 1) {
        if constexpr (StoragePolicy::has_minimal_set_representation)
          simplify();
        else
          generators_.sort();
      }

      return modified;
    }
  }

  /**
   * @brief Sets each generator of the filtration value to the greatest common lower bound between it and the given
   * value.
   *
   * More formally, it pulls the current generator to the cone \f$ \{ y \in \mathbb R^n : y \le x \} \f$
   * originating in \f$ x \f$. The resulting value corresponds to the intersection of both
   * cones: \f$ \mathrm{this} = \min \{ y \in \mathbb R^n : y \le this \} \cap \{ y \in \mathbb R^n : y \le x \} \f$.
   *
   * @tparam GeneratorRange Range of elements convertible to `T`. Must have a begin(), end() and size() method.
   * @param x Range towards to pull. Has to have as many elements than @ref num_parameters().
   * @param exclude_infinite_values If true, values at infinity or minus infinity are not affected.
   * @return true If the filtration value was actually modified.
   * @return false Otherwise.
   */
  template <class GeneratorRange = std::initializer_list<value_type>,
            class = std::enable_if_t<detail::RangeTraits<GeneratorRange>::has_begin ||
                                     detail::RangeTraits<GeneratorRange>::is_multi_filtration>>
  bool pull_to_greatest_common_lower_bound(const GeneratorRange &x, bool exclude_infinite_values = false) {
    std::vector<value_type> pushTo;
    if constexpr (detail::RangeTraits<GeneratorRange>::is_multi_filtration) {
      if (x.num_generators() == 0) return false;
      GUDHI_CHECK(x.num_parameters() == num_parameters(),
                  std::invalid_argument("Wrong range size. Should correspond to the number of parameters."));

      if (x.is_nan()) return false;
      if (x.is_minus_inf()) {
        if (is_minus_inf()) return false;
        generators_ = StoragePolicy::template minus_inf<Co>(num_parameters());
        return true;
      }
      if (x.is_plus_inf()) return false;
      pushTo.resize(x.num_parameters());
      for (size_type p = 0; p < x.num_parameters(); ++p) pushTo[p] = x(0, p);
    } else {
      GUDHI_CHECK(x.size() == num_parameters(),
                  std::invalid_argument("Wrong range size. Should correspond to the number of parameters."));

      pushTo.insert(pushTo.end(), x.begin(), x.end());
    }

    if constexpr (StoragePolicy::has_an_implicit_axis) {
      if (is_nan()) return false;

      auto implicitAxisGen = generators_.get_generator_of_implicit_value(pushTo[StoragePolicy::implicit_axis()]);
      GUDHI_CHECK(implicitAxisGen != StoragePolicy::template null_value<decltype(implicitAxisGen)>(),
                  std::invalid_argument("Value at implicit axis is not valid."));

      bool modified = false;
      auto implGen = static_cast<size_type>(implicitAxisGen);
      for (size_type p = 0; p < num_parameters(); ++p) {
        if (p != StoragePolicy::implicit_axis() &&
            (!exclude_infinite_values || (pushTo[p] != T_m_inf && pushTo[p] != T_inf)) && !detail::_is_nan(pushTo[p])) {
          value_type newVal = _get_default_null_value();

          for (size_type i = 0; i < std::min(implGen, generators_.num_generators()); ++i) {
            auto &v = generators_(i, p);
            if (!detail::_is_nan(v) && v > pushTo[p]) {
              modified = true;
              v = pushTo[p];
            }
          }

          if (implGen < generators_.num_generators()) {
            for (size_type i = implGen; i < generators_.num_generators(); ++i) {
              auto v = generators_(i, p);
              if (!detail::_is_nan(v)) {
                value_type threshold = std::min(v, pushTo[p]);
                if (_dominates(newVal, threshold)) newVal = threshold;
              }
            }
            modified |= generators_(implGen, p) != newVal;
            generators_(implGen, p) = newVal;
          }
        }
      }
      if (implGen < generators_.num_generators()) {
        modified |= (implGen + 1) != generators_.num_generators();
        generators_.set_num_generators(implGen + 1, _get_default_null_value());
      }

      return modified;
    } else {
      bool xIsInf = true, xIsMinusInf = true, xIsNaN = true;
      bool thisIsInf = false, thisIsMinusInf = false, thisIsNaN = false;

      // if one is not finite, we can avoid the heavy simplification process
      _get_infinity_statuses(pushTo, xIsInf, xIsMinusInf, xIsNaN);
      _get_infinity_statuses(thisIsInf, thisIsMinusInf, thisIsNaN);

      if (thisIsMinusInf || thisIsNaN || xIsNaN || xIsInf || (xIsMinusInf && exclude_infinite_values)) return false;

      if (xIsMinusInf) {
        generators_ = StoragePolicy::template minus_inf<Co>(num_parameters());
        return true;
      }

      if (thisIsInf) {
        generators_ = StoragePolicy(pushTo.begin(), pushTo.end());
        return true;
      }

      bool modified = false;

      for (size_type g = 0; g < num_generators(); ++g) {
        for (size_type p = 0; p < num_parameters(); ++p) {
          value_type valX = pushTo[p];
          auto &val = generators_(g, p);
          if (!exclude_infinite_values || (valX != T_inf && valX != T_m_inf)) {
            modified |= val > valX;
            val = valX < val ? valX : val;
          }
        }
      }

      if (modified && num_generators() > 1) {
        if constexpr (StoragePolicy::has_minimal_set_representation)
          simplify();
        else
          generators_.sort();
      }

      return modified;
    }
  }

  /**
   * @brief Projects the generator into the given grid. If @p coordinate is false, the entries are set to
   * the nearest value with the same parameter in the grid. Otherwise, the entries are set to the indices
   * of those nearest values. If an entry in the generator is higher than any value in the grid, this entry
   * is set to infinity if @p coordinate is false and to the grids size at the corresponding parameter otherwise.
   * The grid has to be represented as a vector of ordered ranges of values convertible into `T`. An index
   * \f$ i \f$ of the vector corresponds to the same parameter as the index \f$ i \f$ in a generator of the filtration
   * value. The ranges correspond to the possible values of the parameters, ordered by increasing value, forming
   * therefore all together a 2D grid.
   *
   * @tparam OneDimArray A range of values \f$ U \f$ convertible into `T`. Has to implement
   * a begin, end and operator[] method and a `value_type` definition equal to \f$ U \f$.
   * @param grid Vector of @p OneDimArray with size at least number of filtration parameters. Each array has to be
   * ordered by increasing value.
   * @param coordinate If true, the values are set to the coordinates of the projection in the grid. If false,
   * the values are set to the values at the coordinates of the projection.
   */
  template <typename OneDimArray>
  void project_onto_grid(const std::vector<OneDimArray> &grid, bool coordinate = true) {
    GUDHI_CHECK(
        grid.size() >= num_parameters(),
        std::invalid_argument("The grid should not be smaller than the number of parameters in the filtration value."));

    auto project_generator_value = [&](value_type &val, const OneDimArray &filtration) -> void {
      auto v = static_cast<typename OneDimArray::value_type>(val);
      std::size_t d = std::distance(filtration.begin(), std::lower_bound(filtration.begin(), filtration.end(), v));
      if (d == filtration.size()) {
        val = coordinate ? static_cast<value_type>(d) : T_inf;
      } else {
        if (d != 0 && std::abs(v - filtration[d]) > std::abs(v - filtration[d - 1])) {
          --d;
        }
        val = coordinate ? static_cast<value_type>(d) : static_cast<value_type>(filtration[d]);
      }
    };

#ifdef GUDHI_USE_TBB
    if constexpr (Ensure1Criticality) {
      [[maybe_unused]] value_type &dummy = generators_(0, 0);  // triggers eventual resize
      tbb::parallel_for(size_type(0), num_parameters(), [&](size_type p) -> void {
        if constexpr (StoragePolicy::has_an_implicit_axis) {
          if (p == StoragePolicy::implicit_axis()) return;
        }
        project_generator_value(generators_(0, p), grid[p]);
      });
    } else {
      tbb::parallel_for(size_type(0), num_generators(), [&](size_type g) -> void {
        [[maybe_unused]] value_type &dummy = generators_(g, 0);  // triggers eventual resize
        tbb::parallel_for(size_type(0), num_parameters(), [&](size_type p) -> void {
          if constexpr (StoragePolicy::has_an_implicit_axis) {
            if (p == StoragePolicy::implicit_axis()) return;
          }
          project_generator_value(generators_(g, p), grid[p]);
        });
      });
      // TODO: simplify also in case of coordinate == true ?
      if (!coordinate && num_generators() > 1) {
        if constexpr (StoragePolicy::has_minimal_set_representation)
          simplify();
        else
          generators_.sort();
      }
    }
#else
    for (size_type g = 0; g < num_generators(); ++g) {
      for (size_type p = 0; p < num_parameters(); ++p) {
        if constexpr (StoragePolicy::has_an_implicit_axis) {
          if (p == StoragePolicy::implicit_axis()) continue;
        }
        project_generator_value(generators_(g, p), grid[p]);
      }
    }
    if (!coordinate && num_generators() > 1) {
      if constexpr (StoragePolicy::has_minimal_set_representation)
        simplify();
      else
        generators_.sort();
    }
#endif
  }

  // FUNCTIONALITIES

  /**
   * @private
   */
  template <bool greatest>
  friend Multi_parameter_filtration_value factorize(const Multi_parameter_filtration_value &f) {
    if (f.num_generators() <= 1) return f;

    if constexpr (StoragePolicy::has_an_implicit_axis) {
      std::vector<value_type> gen(f.num_parameters(), greatest ? T_m_inf : T_inf);
      bool isInf = true;
      bool isNaN = true;
      for (size_type p = 0; p < f.num_parameters(); ++p) {
        bool nan = true;
        for (size_type g = 0; g < f.num_generators() && gen[p] != (greatest ? T_inf : T_m_inf); ++g) {
          value_type val = f(g, p);
          if (!detail::_is_nan(val)) {
            nan = false;
            if (val != _get_default_null_value()) {
              isInf = false;
              if constexpr (greatest) {
                gen[p] = val > gen[p] ? val : gen[p];
              } else {
                gen[p] = val < gen[p] ? val : gen[p];
              }
            }
          } else {
            isInf = false;
          }
        }
        if (nan)
          gen[p] = std::numeric_limits<value_type>::quiet_NaN();
        else
          isNaN = false;
      }
      if (isNaN) return nan(f.num_parameters());
      if (isInf) return Co ? minus_inf(f.num_parameters()) : inf(f.num_parameters());
      StoragePolicy result = f.generators_.copy(f.num_parameters(), 0, _get_default_null_value());
      result.template emplace<Co>(gen.begin(), gen.end(), _get_default_null_value());
      return Multi_parameter_filtration_value(std::move(result));
    } else {
      StoragePolicy result(f.num_parameters(), greatest ? T_m_inf : T_inf);
      for (size_type p = 0; p < f.num_parameters(); ++p) {
        bool nan = true;
        for (size_type g = 0; g < f.num_generators() && result(0, p) != (greatest ? T_inf : T_m_inf); ++g) {
          value_type val = f(g, p);
          if (!detail::_is_nan(val)) {
            nan = false;
            if constexpr (greatest) {
              result(0, p) = val > result(0, p) ? val : result(0, p);
            } else {
              result(0, p) = val < result(0, p) ? val : result(0, p);
            }
          }
        }
        if (nan) result(0, p) = std::numeric_limits<value_type>::quiet_NaN();
      }
      return Multi_parameter_filtration_value(std::move(result));
    }
  }

  /**
   * @brief Returns a generator with the minimal values of all parameters in any generator of the given filtration
   * value. That is, the greatest lower bound of all generators.
   */
  friend Multi_parameter_filtration_value factorize_below(const Multi_parameter_filtration_value &f) {
    return factorize<false>(f);
  }

  /**
   * @brief Returns a generator with the maximal values of all parameters in any generator of the given filtration
   * value. That is, the least upper bound of all generators.
   */
  friend Multi_parameter_filtration_value factorize_above(const Multi_parameter_filtration_value &f) {
    return factorize<true>(f);
  }

  /**
   * @brief Computes the coordinates in the given grid, corresponding to the nearest values of the entries
   * in the given filtration value.
   * The grid has to be represented as a 2-dimensional array of ordered values convertible into `OutValue`. An index
   * \f$ i \f$ of the array corresponds to the same parameter as the index \f$ i \f$ in a generator of the filtration
   * value. The ranges correspond to the possible values of the parameters, ordered by increasing value, forming
   * therefore all together a 2D grid.
   *
   * @tparam OutValue Signed arithmetic type. Default value: std::int32_t.
   * @tparam RandomAccessArray A range of values \f$ U \f$ convertible into `T`. Has to implement
   * a begin, end and operator[] method and a `value_type` definition equal to \f$ U \f$.
   * @param f Filtration value to project.
   * @param grid Vector of @p RandomAccessArray with size at least number of filtration parameters. Each array
   * has to be ordered by increasing value.
   * @return Filtration value \f$ out \f$ whose entry correspond to the indices of the projected values. That is,
   * the projection of \f$ f(g,p) \f$ is \f$ grid[p][out(g,p)] \f$.
   */
  template <typename OutValue = std::int32_t, class RandomAccessArray = std::vector<value_type>>
  friend auto compute_coordinates_in_grid(Multi_parameter_filtration_value f,
                                          const std::vector<RandomAccessArray> &grid) {
    // TODO: by replicating the code of "project_onto_grid", this could be done with just one copy
    // instead of two. But it is not clear if it is really worth it, i.e., how much the change in type is really
    // necessary in the use cases. To see later.
    f.project_onto_grid(grid);
    if constexpr (std::is_same_v<OutValue, value_type>) {
      return f;
    } else {
      return f.as_type<OutValue>();
    }
  }

  /**
   * @brief Computes the values in the given grid corresponding to the coordinates given by the given filtration
   * value. That is, if \f$ out \f$ is the result, \f$ out(g,p) = grid[p][f(g,p)] \f$. Assumes therefore, that the
   * values stored in the filtration value corresponds to indices existing in the given grid.
   *
   * @tparam RandomAccessArray A range of values convertible into `U`. Has to implement
   * a size and operator[] method and a `value_type` definition.
   * @tparam U Signed arithmetic type. Default: `RandomAccessArray::value_type`.
   * @param f Filtration value storing coordinates compatible with `grid`.
   * @param grid Vector of @p RandomAccessArray with size at least number of filtration parameters. Each array
   * has to be ordered by increasing value.
   * @return Filtration value \f$ out \f$ whose entry correspond to \f$ out(g,p) = grid[p][f(g,p)] \f$.
   */
  template <typename U, class RandomAccessArray>
  friend auto evaluate_coordinates_in_grid(const Multi_parameter_filtration_value &f,
                                           const std::vector<RandomAccessArray> &grid) {
    GUDHI_CHECK(grid.size() >= f.num_parameters(),
                std::invalid_argument(
                    "The size of the grid should correspond to the number of parameters in the filtration value."));

    using SP = typename StoragePolicy::template As_type<U>;
    U grid_inf = SP::template T_inf<>;
    U grid_m_inf = SP::template T_m_inf<>;
    auto get_filtration_value = [grid_inf](const Multi_parameter_filtration_value &f) {
      if constexpr (StoragePolicy::has_an_implicit_axis) {
        return f.copy<U>(f.num_parameters(), 0);
      } else {
        return Multi_parameter_filtration_value<SP, Co, Ensure1Criticality>(f.num_parameters(), grid_inf);
      }
    };

    Multi_parameter_filtration_value<SP, Co, Ensure1Criticality> out = get_filtration_value(f);
    out.reserve(f.num_generators());
    std::vector<U> tmpGen(f.num_parameters());

    for (size_type g = 0; g < f.num_generators(); ++g) {
      for (size_type p = 0; p < f.num_parameters(); ++p) {
        const RandomAccessArray &filtration = grid[p];
        const value_type &c = f(g, p);
        if (detail::_is_nan(c))
          tmpGen[p] = std::numeric_limits<U>::quiet_NaN();
        else if (c < 0)
          tmpGen[p] = grid_m_inf;
        else if (c == T_inf || static_cast<std::size_t>(c) >= filtration.size())
          tmpGen[p] = grid_inf;
        else
          tmpGen[p] = static_cast<U>(filtration[c]);
      }
      out.add_generator(tmpGen);
    }
    return out;
  }

  /**
   * @brief Computes the values in the given grid corresponding to the coordinates given by the given filtration
   * value. That is, if \f$ out \f$ is the result, \f$ out(g,p) = grid[p][f(g,p)] \f$. Assumes therefore, that the
   * values stored in the filtration value corresponds to indices existing in the given grid.
   *
   * It is equivalent to the other version when `U` corresponds to `RandomAccessArray::value_type`.
   *
   * @tparam RandomAccessArray A range with size and operator[] method and a `value_type` definition.
   * @param f Filtration value storing coordinates compatible with `grid`.
   * @param grid Vector of @p RandomAccessArray with size at least number of filtration parameters. Each array
   * has to be ordered by increasing value.
   * @return Filtration value \f$ out \f$ whose entry correspond to \f$ out(g,p) = grid[p][f(g,p)] \f$.
   */
  template <class RandomAccessArray>
  friend auto evaluate_coordinates_in_grid(const Multi_parameter_filtration_value &f,
                                           const std::vector<RandomAccessArray> &grid) {
    return evaluate_coordinates_in_grid<typename RandomAccessArray::value_type, RandomAccessArray>(f, grid);
  }

  // UTILITIES

  /**
   * @brief Outstream operator.
   */
  friend std::ostream &operator<<(std::ostream &stream, const Multi_parameter_filtration_value &f) {
    const size_type num_gen = f.num_generators();
    const size_type num_param = f.num_parameters();

    stream << "( k = " << num_gen << " ) ( p = " << num_param << " ) [ ";
    for (size_type g = 0; g < num_gen; ++g) {
      stream << "[";
      for (size_type p = 0; p < num_param; ++p) {
        stream << f(g, p);
        if (p < num_param - 1) stream << ", ";
      }
      stream << "]";
      if (g < num_gen - 1) stream << "; ";
    }
    stream << " ]";

    return stream;
  }

  /**
   * @brief Instream operator.
   */
  friend std::istream &operator>>(std::istream &stream, Multi_parameter_filtration_value &f) {
    size_type num_gen;
    size_type num_param;
    char delimiter;
    stream >> delimiter;  // (
    stream >> delimiter;  // k
    stream >> delimiter;  // =
    if (delimiter != '=')
      throw std::invalid_argument("Invalid incoming stream format for Multi_parameter_filtration_value.");
    stream >> num_gen;
    if (!stream.good())
      throw std::invalid_argument("Invalid incoming stream format for Multi_parameter_filtration_value.");
    stream >> delimiter;  // )
    stream >> delimiter;  // (
    stream >> delimiter;  // p
    stream >> delimiter;  // =
    stream >> num_param;
    if (!stream.good())
      throw std::invalid_argument("Invalid incoming stream format for Multi_parameter_filtration_value.");
    f = Multi_parameter_filtration_value(num_param);
    f.generators_.set_num_generators(num_gen, _get_default_null_value());
    stream >> delimiter;  // )
    stream >> delimiter;  // [
    if (delimiter != '[')
      throw std::invalid_argument("Invalid incoming stream format for Multi_parameter_filtration_value.");
    if (num_gen == 0) return stream;
    for (size_type i = 0; i < num_gen; ++i) {
      stream >> delimiter;  // [
      for (size_type j = 0; j < num_param; ++j) {
        f(i, j) = detail::_get_value<value_type>(stream);
        if (!stream.good())
          throw std::invalid_argument("Invalid incoming stream format for Multi_parameter_filtration_value.");
        stream >> delimiter;  // , or last ]
      }
      stream >> delimiter;  // ; or last ]
    }
    if (delimiter != ']')
      throw std::invalid_argument("Invalid incoming stream format for Multi_parameter_filtration_value.");

    return stream;
  }

  /**
   * @brief Returns true if and only if the given filtration value is at plus infinity.
   */
  friend bool is_positive_infinity(const Multi_parameter_filtration_value &f) { return f.is_plus_inf(); }

  /**
   * @brief Adds the generators of the second argument to the first argument. If `Ensure1Criticality` is true,
   * the method assumes that the two filtration values are comparable, that is, that the result of the union is also
   * 1-critical. A check for this is only done in Debug Mode, as it is costly.
   *
   * @param f1 Filtration value to modify.
   * @param f2 Filtration value to merge with the first one. Should have the same number of parameters than the other.
   * @return true If the first argument was actually modified.
   * @return false Otherwise.
   */
  friend bool unify_lifetimes(Multi_parameter_filtration_value &f1, const Multi_parameter_filtration_value &f2) {
    GUDHI_CHECK(f1.num_parameters() == f2.num_parameters(),
                std::invalid_argument("Cannot unify two filtration values with different number of parameters."));

    // TODO: verify if this really makes a differences in the 1-critical case, otherwise just keep the general case
    // if general case is kept: add (num_gen == 1) test to throw if unification is not 1-critical anymore.
    if constexpr (Ensure1Criticality) {
      // WARNING: costly check
      GUDHI_CHECK(
          f1 <= f2 || f2 <= f1,
          std::invalid_argument("When 1-critical only, two non-comparable filtration values cannot be unified."));

      if constexpr (Co) {
        return f1.push_to_least_common_upper_bound(f2);
      } else {
        return f1.pull_to_greatest_common_lower_bound(f2);
      }
    } else {
      bool modified = false;
      for (size_type g = 0; g < f2.num_generators(); ++g) {
        auto start = f2.generators_.begin(g);
        auto end = f2.generators_.end(g);
        modified |= f1.add_generator(start, end);
      }
      return modified;
    }
  }

  /**
   * @brief Stores in the first argument the origins of the cones in the intersection of the positive
   * (negative if `Co` is true) cones generated by the two arguments.
   *
   * If @ref StoragePolicy::has_an_implicit_axis is true, then both values must have same implicit values,
   * otherwise the implicit values in the result will likely be wrong.
   *
   * @param f1 First set of cones which will be modified.
   * @param f2 Second set of cones. Should have the same number of parameters than the first one.
   * @return true If the first argument was actually modified.
   * @return false Otherwise.
   */
  friend bool intersect_lifetimes(Multi_parameter_filtration_value &f1, const Multi_parameter_filtration_value &f2) {
    GUDHI_CHECK(f1.num_parameters() == f2.num_parameters(),
                std::invalid_argument("Cannot intersect two filtration values with different number of parameters."));

    if constexpr (StoragePolicy::has_an_implicit_axis) {
      GUDHI_CHECK(f1.generators_.has_same_implicit_values(f2.generators_),
                  std::invalid_argument("Cannot intersect two filtration values with non-aligned implicit axes."));
    }

    bool f2IsInf = false, f2IsMinusInf = false, f2IsNaN = false;
    bool f1IsInf = false, f1IsMinusInf = false, f1IsNaN = false;

    // if one is not finite, we can avoid the heavy simplification process
    f1._get_infinity_statuses(f1IsInf, f1IsMinusInf, f1IsNaN);
    f2._get_infinity_statuses(f2IsInf, f2IsMinusInf, f2IsNaN);

    if (f1IsNaN || f2IsNaN) return false;

    // inf cases first to avoid costly g1 * g2 check
    if constexpr (Co) {
      if (f1IsInf) {
        if (f2IsInf) return false;
        f1 = f2;
        return true;
      }
      if (f1IsMinusInf) {
        return false;
      }
    } else {
      if (f1IsMinusInf) {
        if (f2IsMinusInf) return false;
        f1 = f2;
        return true;
      }
      if (f1IsInf) {
        return false;
      }
    }

    if constexpr (Ensure1Criticality) {
      if constexpr (Co) {
        return f1.pull_to_greatest_common_lower_bound(f2);
      } else {
        return f1.push_to_least_common_upper_bound(f2);
      }
    } else if constexpr (StoragePolicy::has_an_implicit_axis) {
      bool modified = false;
      for (size_type p = 0; p < f1.num_parameters(); ++p) {
        if (p != StoragePolicy::implicit_axis()) {
          value_type threshold1 = f1(0, p);
          value_type threshold2 = f2(0, p);
          for (size_type g = 0; g < std::max(f1.num_generators(), f2.num_generators()); ++g) {
            if (g < f1.num_generators())
              threshold1 = _strictly_dominates(threshold1, f1(g, p)) ? f1(g, p) : threshold1;
            else {
              f1.set_num_generators(g + 1);
              modified = true;
            }
            if (g < f2.num_generators()) threshold2 = _strictly_dominates(threshold2, f2(g, p)) ? f2(g, p) : threshold2;
            if (_strictly_dominates(threshold2, threshold1)) {
              if (f1(g, p) != threshold2) modified = true;
              f1(g, p) = threshold2;
            } else {
              if (f1(g, p) != threshold1) modified = true;
              f1(g, p) = threshold1;
            }
          }
        }
      }
      return modified;
    } else {
      const size_type num_param = f1.num_parameters();
      Multi_parameter_filtration_value res = Co ? minus_inf(num_param) : inf(num_param);
      std::vector<value_type> newGen(num_param);
      // TODO: see if the order can be used to avoid g1 * g2 add_generator and
      // perhaps even to replace add_generator by add_guaranteed_generator
      for (size_type g1 = 0; g1 < f1.num_generators(); ++g1) {
        for (size_type g2 = 0; g2 < f2.num_generators(); ++g2) {
          for (size_type p = 0; p < num_param; ++p) {
            if constexpr (Co) {
              newGen[p] = std::min(f1(g1, p), f2(g2, p));
            } else {
              newGen[p] = std::max(f1(g1, p), f2(g2, p));
            }
          }
          res.add_generator(newGen);
        }
      }
      swap(f1, res);

      return f1 != res;
    }
  }

  /**
   * @brief Serialize given value into the buffer at given pointer.
   *
   * @param value Value to serialize.
   * @param start Pointer to the start of the space in the buffer where to store the serialization.
   * @return End position of the serialization in the buffer.
   */
  friend char *serialize_value_to_char_buffer(const Multi_parameter_filtration_value &value, char *start) {
    return serialize_value_to_char_buffer(value.generators_, start);
  }

  /**
   * @brief Deserialize the value from a buffer at given pointer and stores it in given value.
   *
   * @param value Value to fill with the deserialized filtration value.
   * @param start Pointer to the start of the space in the buffer where the serialization is stored.
   * @return End position of the serialization in the buffer.
   */
  friend const char *deserialize_value_from_char_buffer(Multi_parameter_filtration_value &value, const char *start) {
    return deserialize_value_from_char_buffer(value.generators_, start);
  }

  /**
   * @brief Returns the serialization size of the given filtration value.
   */
  friend std::size_t get_serialization_size_of(const Multi_parameter_filtration_value &value) {
    return get_serialization_size_of(value.generators_);
  }

 private:
  StoragePolicy generators_;

  constexpr static value_type _get_default_value() { return Co ? T_inf : T_m_inf; }

  constexpr static value_type _get_default_null_value() { return Co ? T_m_inf : T_inf; }

  template <bool strict>
  static bool _is_lower_than_with_implicit_axis(const Multi_parameter_filtration_value &a,
                                                const Multi_parameter_filtration_value &b) {
    if constexpr (StoragePolicy::has_an_implicit_axis) {
      // works because values at implicit axis are unique and the generators are ordered
      // by increasing value at that axis
      GUDHI_CHECK(StoragePolicy::implicit_axis() < b.num_parameters(),
                  std::runtime_error("Implicit axis out of bound??"));

      if (a.generators_.has_same_implicit_values(b.generators_)) {
        for (size_type p = 0; p < b.num_parameters(); ++p) {
          if (p != StoragePolicy::implicit_axis()) {
            value_type threshold = a(0, p);
            for (size_type g = 0; g < b.num_generators(); ++g) {
              if (detail::_is_nan(b(g, p))) return false;
              if (g >= a.num_generators() && _strictly_dominates(threshold, b(g, p))) return false;
              if (g < a.num_generators()) {
                if (detail::_is_nan(a(g, p))) return false;
                threshold = _strictly_dominates(threshold, a(g, p)) ? a(g, p) : threshold;
                if constexpr (strict)
                  if (a(g, p) == threshold && b(g, p) == threshold) return false;
              }
              if (_strictly_dominates(threshold, b(g, p))) return false;
            }
          }
        }
        return true;
      }

      // as the implicit values are not the same, we have to align them here.
      // still assumes that the implicit value is always an integer and unique.
      // there is probably also a way to do this with floating point, by just using the relative order, but it is a
      // bit of a headache for something we don't need for now.
      // but if we end up implementing this, we could probably also use it for
      // StoragePolicy::has_lexicographical_storage (replacing g x g x p complexity with p x g)
      struct Viewer {
        StoragePolicy const *f;
        size_type p;
        int curr;
        int end;
        size_type currG;
        value_type max;

        Viewer(const StoragePolicy &s, size_type p_, value_type m)
            : f(&s),
              p(p_),
              curr(s(0, StoragePolicy::implicit_axis())),
              end(s(s.num_generators() - 1, StoragePolicy::implicit_axis()) + 1),
              currG(0),
              max(m) {
          while (currG < s.num_generators() && s(currG, p) == max) {
            ++currG;
            if (currG < s.num_generators())
              curr = s(currG, StoragePolicy::implicit_axis());
            else
              curr = end;
          }
        }
        value_type get_value() const {
          if (static_cast<int>((*f)(currG, StoragePolicy::implicit_axis())) == curr) return (*f)(currG, p);
          return max;
        }
        void next() {
          if (static_cast<int>((*f)(currG, StoragePolicy::implicit_axis())) == curr) ++currG;
          ++curr;
        }
        [[nodiscard]] bool is_end() const { return curr >= end; }
      };

      for (size_type p = 0; p < b.num_parameters(); ++p) {
        if (p != StoragePolicy::implicit_axis()) {
          Viewer viewA(a.get_underlying_policy(), p, _get_default_null_value());
          Viewer viewB(b.get_underlying_policy(), p, _get_default_null_value());

          if (viewA.is_end() && viewB.is_end()) return !strict;  // both are inf or both -inf
          if (viewA.is_end()) return Co;
          if (viewB.is_end()) return !Co;

          if (viewA.curr > viewB.curr) return false;

          value_type threshold = viewA.get_value();  // can't be max
          while (!viewA.is_end() && viewA.curr < viewB.curr) {
            auto v = viewA.get_value();
            threshold = _strictly_dominates(threshold, v) ? v : threshold;
            viewA.next();
          }

          while (!viewB.is_end()) {
            auto vB = viewB.get_value();
            if (detail::_is_nan(vB)) return false;
            if (viewA.is_end() && _strictly_dominates(threshold, vB)) return false;
            if (!viewA.is_end()) {
              auto vA = viewA.get_value();
              if (detail::_is_nan(vA)) return false;
              if (vA == viewA.max && vB == viewB.max) {
                viewA.next();
                viewB.next();
                continue;
              }
              threshold = _strictly_dominates(threshold, vA) ? vA : threshold;
              if constexpr (strict)
                if (vA == threshold && vB == threshold) return false;
              viewA.next();
              if (_strictly_dominates(threshold, vB)) return false;
            }
            viewB.next();
          }
        }
      }
      return true;
    }
  }

  template <bool strict>
  static bool _is_lower_than_other(const Multi_parameter_filtration_value &a,
                                   const Multi_parameter_filtration_value &b) {
    if constexpr (!StoragePolicy::has_an_implicit_axis) {
      // TODO: try optimization above for StoragePolicy::has_lexicographical_storage
      for (size_type i = 0; i < b.num_generators(); ++i) {
        // for each generator in b, verify if it is in the cone of at least one generator of a
        bool isContained = false;
        for (size_type j = 0; j < a.num_generators() && !isContained; ++j) {
          // lexicographical order, so if a[j][0] strictly dom b[j][0], than a[j'] can never contain b[i]
          // for all j' > j.
          if constexpr (strict) {
            if constexpr (StoragePolicy::has_lexicographical_storage)
              if (_dominates(a.generators_(j, 0), b.generators_(i, 0))) return false;
            isContained = _strictly_contains(a.generators_, j, b.generators_, i);
          } else {
            if constexpr (StoragePolicy::has_lexicographical_storage)
              if (_strictly_dominates(a.generators_(j, 0), b.generators_(i, 0))) return false;
            isContained = _contains(a.generators_, j, b.generators_, i);
          }
        }
        if (!isContained) return false;
      }
      return true;
    }
  }

  template <bool strict>
  static bool _is_lower_than(const Multi_parameter_filtration_value &a, const Multi_parameter_filtration_value &b) {
    GUDHI_CHECK(a.num_parameters() == b.num_parameters(),
                std::invalid_argument("Only filtration values with same number of parameters can be compared."));

    if (a.num_generators() == 0 || b.num_generators() == 0) return false;
    if (a.is_nan() || b.is_nan()) return false;
    if (&a == &b) return !strict;

    if (b.num_parameters() == 1) {
      if constexpr (strict)
        return a(0, 0) < b(0, 0);
      else
        return a(0, 0) <= b(0, 0);
    }

    if constexpr (StoragePolicy::has_an_implicit_axis) {
      return _is_lower_than_with_implicit_axis<strict>(a, b);
    } else {
      return _is_lower_than_other<strict>(a, b);
    }
  }

  /**
   * @brief Verifies if @p b is strictly contained in the positive cone originating in `a`.
   */
  static bool _strictly_contains(const StoragePolicy &a, size_type g_a, const StoragePolicy &b, size_type g_b) {
    bool isSame = true;
    for (auto i = 0U; i < a.num_parameters(); ++i) {
      value_type a_i, b_i;
      if constexpr (Co) {
        a_i = b(g_b, i);
        b_i = a(g_a, i);
      } else {
        a_i = a(g_a, i);
        b_i = b(g_b, i);
      }
      if (a_i > b_i || detail::_is_nan(a_i) || detail::_is_nan(b_i)) return false;
      if (isSame && a_i != b_i) isSame = false;
    }
    return !isSame;
  }

  /**
   * @brief Verifies if @p b is contained in the positive cone originating in `a`.
   */
  static bool _contains(const StoragePolicy &a, size_type g_a, const StoragePolicy &b, size_type g_b) {
    for (size_type i = 0U; i < a.num_parameters(); ++i) {
      value_type a_i, b_i;
      if constexpr (Co) {
        a_i = b(g_b, i);
        b_i = a(g_a, i);
      } else {
        a_i = a(g_a, i);
        b_i = b(g_b, i);
      }
      if (a_i > b_i || (!detail::_is_nan(a_i) && detail::_is_nan(b_i)) ||
          (detail::_is_nan(a_i) && !detail::_is_nan(b_i)))
        return false;
    }
    return true;
  }

  /**
   * @brief Verifies if the first element of @p b strictly dominates the first element of `a`.
   */
  static bool _strictly_dominates(value_type a, value_type b) {
    if constexpr (Co) {
      return a < b;
    } else {
      return a > b;
    }
  }

  /**
   * @brief Verifies if the first element of @p b dominates the first element of `a`.
   */
  static bool _dominates(value_type a, value_type b) {
    if constexpr (Co) {
      return a <= b;
    } else {
      return a >= b;
    }
  }

  /**
   * @brief Applies operation on the elements of the filtration value.
   */
  template <class ValueRange, class F, class = std::enable_if_t<detail::RangeTraits<ValueRange>::has_begin>>
  void _apply_operation(const ValueRange &range, F &&operate) {
    if (is_nan()) return;
    bool isNaN = true;
    for (size_type g = 0; g < num_generators(); ++g) {
      auto it = range.begin();
      for (size_type p = 0; p < num_parameters() && it != range.end(); ++p) {
        if constexpr (StoragePolicy::has_an_implicit_axis) {
          if (StoragePolicy::implicit_axis() == p) {
            ++it;
            continue;
          }
        }
        if (std::forward<F>(operate)(generators_(g, p), *it)) isNaN = false;
        ++it;
      }
    }
    if (isNaN && StoragePolicy::has_quiet_NaN) generators_ = StoragePolicy::nan(num_parameters());
  }

  template <class OtherStoragePolicy, bool OtherCo, bool OtherEnsure1Criticality, class F>
  void _apply_operation(
      const Multi_parameter_filtration_value<OtherStoragePolicy, OtherCo, OtherEnsure1Criticality> &range,
      F &&operate) {
    if (is_nan()) return;
    if (range.is_nan()) {
      if constexpr (StoragePolicy::has_quiet_NaN) generators_ = StoragePolicy::nan(num_parameters());
      return;
    }
    if (range.num_generators() == 0) return;
    bool isNaN = true;
    for (size_type g = 0; g < num_generators(); ++g) {
      for (size_type p = 0; p < std::min(num_parameters(), range.num_parameters()); ++p) {
        if constexpr (StoragePolicy::has_an_implicit_axis) {
          if (StoragePolicy::implicit_axis() == p) continue;
        }
        if (std::forward<F>(operate)(generators_(g, p), range(0, p))) isNaN = false;
      }
    }
    if (isNaN && StoragePolicy::has_quiet_NaN) generators_ = StoragePolicy::nan(num_parameters());
  }

  /**
   * @brief Applies operation on the elements of the filtration value.
   */
  template <class F>
  void _apply_operation(F &&operate) {
    if (is_nan()) return;
    bool isNaN = true;
    for (size_type g = 0; g < num_generators(); ++g) {
      for (size_type p = 0; p < num_parameters(); ++p) {
        if constexpr (StoragePolicy::has_an_implicit_axis) {
          if (StoragePolicy::implicit_axis() == p) continue;
        }
        if (std::forward<F>(operate)(generators_(g, p))) isNaN = false;
      }
    }
    if (isNaN && StoragePolicy::has_quiet_NaN) generators_ = StoragePolicy::nan(num_parameters());
  }

  template <class GeneratorRange>
  static void _get_infinity_statuses(const GeneratorRange &r, bool &isInf, bool &isMinusInf, bool &isNaN) {
    for (auto itB = r.begin(); itB != r.end(); ++itB) {
      value_type v = *itB;
      if (v != T_inf) isInf = false;
      if (v != T_m_inf) isMinusInf = false;
      if (!detail::_is_nan(v)) isNaN = false;
      if (!isInf && !isMinusInf && !isNaN) return;
    }
  }

  void _get_infinity_statuses(bool &isInf, bool &isMinusInf, bool &isNaN) const {
    if (is_nan()) {
      isNaN = true;
      return;
    }
    if (is_plus_inf()) {
      isInf = true;
      return;
    }
    if (is_minus_inf()) {
      isMinusInf = true;
      return;
    }
  }

  enum class Rel : std::uint8_t { EQUAL, DOMINATES, IS_DOMINATED, NONE };

  template <class Iterator>
  Rel _get_domination_relation(size_type g, Iterator it) const {
    bool equal = true;
    bool allGreater = true;
    bool allSmaller = true;
    bool allNaNA = true;
    bool allNaNB = true;
    for (size_type p = 0; p < num_parameters(); ++p) {
      if (generators_(g, p) < static_cast<value_type>(*it)) {
        if (!allSmaller) return Rel::NONE;
        equal = false;
        allGreater = false;
      } else if (generators_(g, p) > static_cast<value_type>(*it)) {
        if (!allGreater) return Rel::NONE;
        equal = false;
        allSmaller = false;
      }
      if (!detail::_is_nan(generators_(g, p))) allNaNA = false;
      if (!detail::_is_nan(*it)) allNaNB = false;
      ++it;
    }
    if (allNaNA || allNaNB) return Rel::IS_DOMINATED;
    if (equal) return Rel::EQUAL;

    if constexpr (Co) {
      if (allSmaller) return Rel::DOMINATES;
      return Rel::IS_DOMINATED;
    } else {
      if (allGreater) return Rel::DOMINATES;
      return Rel::IS_DOMINATED;
    }
  }

  /**
   * @brief Verifies how x can be added as a new generator with respect to an already existing generator, represented
   * by `indices[curr]`. If x is dominated by or is equal to `indices[curr]`, it cannot be added. If it dominates
   * `indices[curr]`, it has to replace `indices[curr]`. If there is no relation between both, `indices[curr]` has
   * no influence on the addition of x.
   *
   * Assumes between 'curr' and 'end' everything is simplified:
   * no nan values and if there is an inf/-inf, then 'end - curr == 1'.
   */
  template <class Iterator>
  bool _generator_can_be_added(Iterator x, size_type curr, size_type &end) {
    while (curr != end) {
      Rel res = _get_domination_relation(curr, x);
      if (res == Rel::IS_DOMINATED || res == Rel::EQUAL) return false;  // x dominates or is equal
      if (res == Rel::DOMINATES) {                                      // x is dominated
        --end;
        generators_.swap_generators(curr, end);
      } else {  // no relation
        ++curr;
      }
    }
    return true;
  }
};

}  // namespace Gudhi::multi_filtration

namespace std {

template <class StoragePolicy, bool Co, bool Ensure1Criticality>
class numeric_limits<Gudhi::multi_filtration::Multi_parameter_filtration_value<StoragePolicy, Co, Ensure1Criticality>> {
 public:
  using T = typename StoragePolicy::value_type;
  using Filtration_value =
      Gudhi::multi_filtration::Multi_parameter_filtration_value<StoragePolicy, Co, Ensure1Criticality>;

  static constexpr bool has_infinity = StoragePolicy::template has_infinity<Co>;
  static constexpr bool has_quiet_NaN = StoragePolicy::has_quiet_NaN;

  static constexpr Filtration_value infinity(std::size_t p = 1) { return Filtration_value::inf(p); };

  // non-standard
  static constexpr Filtration_value minus_infinity(std::size_t p = 1) { return Filtration_value::minus_inf(p); };

  static constexpr Filtration_value max() noexcept(false) {
    throw std::logic_error(
        "The max value cannot be represented with no finite numbers of parameters."
        "Use `max(numberOfParameters)` instead");
  };

  static constexpr Filtration_value max(std::size_t p) {
    if constexpr (has_infinity) {
      return Filtration_value(p, std::numeric_limits<T>::max());
    } else {
      throw std::logic_error("No biggest value possible for Co-filtrations yet.");
    }
  };

  static constexpr Filtration_value lowest(std::size_t p = 1) { return Filtration_value::minus_inf(p); };

  static constexpr Filtration_value quiet_NaN(std::size_t p = 1) {
    if constexpr (has_quiet_NaN) {
      return Filtration_value::nan(p);
    } else {
      throw std::logic_error("Does not have a NaN value.");
    }
  };
};

}  // namespace std

#endif  // MF_MULTI_PARAMETER_FILTRATION_VALUE_H_
