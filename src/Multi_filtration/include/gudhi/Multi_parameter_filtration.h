/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber, David Loiseaux
 *
 *    Copyright (C) 2024-25 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file Multi_parameter_filtration.h
 * @author Hannah Schreiber, David Loiseaux
 * @brief Contains the @ref Gudhi::multi_filtration::Multi_parameter_filtration class.
 */

#ifndef MF_MULTI_PARAMETER_FILTRATION_H_
#define MF_MULTI_PARAMETER_FILTRATION_H_

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
#include <numeric>      //std::iota
#include <vector>
#include <initializer_list>

#include <gudhi/Debug_utils.h>
#include <gudhi/simple_mdspan.h>
#include <gudhi/Multi_filtration/multi_filtration_utils.h>

namespace Gudhi::multi_filtration {

// declaration needed pre C++20 for friends with templates defined inside a class
template <typename U>
U compute_linear_projection();
template <typename U>
U compute_euclidean_distance_to();
template <typename U>
U compute_norm();
template <typename OutValue, typename U>
void compute_coordinates_in_grid();
template <typename U>
void evaluate_coordinates_in_grid();
template <bool inverse>
bool is_strict_less_than_lexicographically();
template <bool inverse>
bool is_less_or_equal_than_lexicographically();

/**
 * @class Multi_parameter_filtration Multi_parameter_filtration.h gudhi/Multi_parameter_filtration.h
 * @ingroup multi_filtration
 *
 * @brief Class encoding the different generators, i.e., apparition times, of a \f$ k \f$-critical
 * \f$\mathbb R^n\f$-filtration value. E.g., the filtration value of a simplex, or, of the algebraic generator of a
 * module presentation. The encoding is compacted into a single vector, so if a lot of non trivial modifications are
 * done (that not only consists of simply adding new generators at the end of the vector), it is probably preferable
 * to use @ref Dynamic_multi_parameter_filtration instead. Implements the concept @ref FiltrationValue of the
 * @ref Gudhi::Simplex_tree and the concept @ref Gudhi::multi_persistence::MultiFiltrationValue.
 *
 * @details Overloads `std::numeric_limits` such that:
 * - `std::numeric_limits<Multi_parameter_filtration>::has_infinity` returns `true`,
 * - `std::numeric_limits<Multi_parameter_filtration>::has_quiet_NaN` returns `std::numeric_limits<T>::has_quiet_NaN`,
 * - `std::numeric_limits<Multi_parameter_filtration>::infinity(int)` returns
 * @ref Multi_parameter_filtration::inf(int) "",
 * - `std::numeric_limits<Multi_parameter_filtration>::minus_infinity(int)` returns
 * @ref Multi_parameter_filtration::minus_inf(int) "",
 * - `std::numeric_limits<Multi_parameter_filtration>::max(int num_param)` returns a @ref Multi_parameter_filtration
 * with one generator of `num_param` parameters evaluated at value `std::numeric_limits<T>::max()`,
 * - `std::numeric_limits<Multi_parameter_filtration>::quiet_NaN(int)` returns
 * @ref Multi_parameter_filtration::nan(int) if `std::numeric_limits<Multi_parameter_filtration>::has_quiet_NaN`
 * and throws otherwise.
 *
 * Multi-critical filtrations are filtrations such that the lifetime of each object is union of positive cones in
 * \f$\mathbb R^n\f$, e.g.,
 *  - \f$ \{ x \in \mathbb R^2 : x \ge (1,2)\} \cap \{ x \in \mathbb R^2 : x \ge (2,1)\} \f$ is finitely critical,
 *    and more particularly 2-critical, while
 *  - \f$ \{ x \in \mathbb R^2 : x \ge \mathrm{epigraph}(y \mapsto e^{-y})\} \f$ is not.
 *
 * @tparam T Arithmetic type of an entry for one parameter of a filtration value. Has to be **signed** and
 * to implement `std::isnan(T)`, `std::numeric_limits<T>::has_quiet_NaN`, `std::numeric_limits<T>::quiet_NaN()`,
 * `std::numeric_limits<T>::has_infinity`, `std::numeric_limits<T>::infinity()` and `std::numeric_limits<T>::max()`.
 * If `std::numeric_limits<T>::has_infinity` returns `false`, a call to `std::numeric_limits<T>::infinity()`
 * can simply throw. Examples are the native types `double`, `float` and `int`.
 * @tparam Co If `true`, reverses the poset order, i.e., the order \f$ \le \f$  in \f$ \mathbb R^n \f$ becomes
 * \f$ \ge \f$. That is, the positive cones representing a lifetime become all negative instead.
 * @tparam Ensure1Criticality If `true`, the methods ensure that the filtration value is always 1-critical by throwing
 * or refusing to compile if a modification increases the number of generators.
 */
template <typename T, bool Co = false, bool Ensure1Criticality = false>
class Multi_parameter_filtration
{
 private:
  using view_extents = extents<std::size_t, Gudhi::dynamic_extent, Gudhi::dynamic_extent>;
  using Viewer = Gudhi::Simple_mdspan<T, view_extents>;

 public:
  using Underlying_container = std::vector<T>; /**< Underlying container for values. */

  // CONSTRUCTORS

  /**
   * @brief Default constructor. Builds filtration value with one generator and given number of parameters.
   * If Co is false, all values are at -inf, if Co is true, all values are at +inf.
   *
   * @param number_of_parameters If negative, takes the default value instead. Default value: 2.
   */
  Multi_parameter_filtration(int number_of_parameters = 2)
      : generators_(number_of_parameters < 0 ? 2 : number_of_parameters, _get_default_value()),
        generator_view_(generators_.data(), generators_.empty() ? 0 : 1, generators_.size())
  {}

  /**
   * @brief Builds filtration value with one generator and given number of parameters.
   * All values are initialized at the given value.
   *
   * @param number_of_parameters If negative, is set to 2 instead.
   * @param value Initialization value for every value in the generator.
   */
  Multi_parameter_filtration(int number_of_parameters, T value)
      : generators_(number_of_parameters < 0 ? 2 : number_of_parameters, value),
        generator_view_(generators_.data(), generators_.empty() ? 0 : 1, generators_.size())
  {}

  /**
   * @brief Builds filtration value with one generator that is initialized with the given range. The number of
   * parameters are therefore deduced from the length of the range.
   *
   * @tparam ValueRange Range of types convertible to `T`. Should have a begin() and end() method.
   * @param range Values of the generator.
   */
  template <class ValueRange = std::initializer_list<T>, class = std::enable_if_t<RangeTraits<ValueRange>::has_begin> >
  Multi_parameter_filtration(const ValueRange &range)
      : generators_(range.begin(), range.end()),
        generator_view_(generators_.data(), generators_.empty() ? 0 : 1, generators_.size())
  {}

  /**
   * @brief Builds filtration value with one generator that is initialized with the given range. The range is
   * determined from the two given iterators. The number of parameters are therefore deduced from the distance
   * between the two.
   *
   * @tparam Iterator Iterator type that has to satisfy the requirements of standard LegacyInputIterator and
   * dereferenced elements have to be convertible to `T`.
   * @param it_begin Iterator pointing to the start of the range.
   * @param it_end Iterator pointing to the end of the range.
   */
  template <class Iterator>
  Multi_parameter_filtration(Iterator it_begin, Iterator it_end)
      : generators_(it_begin, it_end),
        generator_view_(generators_.data(), generators_.empty() ? 0 : 1, generators_.size())
  {}

  /**
   * @brief Builds filtration value with given number of parameters and values from the given range. Lets \f$ p \f$
   * be the number of parameters. The \f$ p \f$ first elements of the range have to correspond to the first generator,
   * the \f$ p \f$ next elements to the second generator and so on... So the length of the range has to be a multiple
   * of \f$ p \f$ and the number of generators will be \f$ length / p \f$. The range is represented by two iterators.
   *
   * @tparam Iterator Iterator type that has to satisfy the requirements of standard LegacyInputIterator and
   * dereferenced elements have to be convertible to `T`.
   * @param it_begin Iterator pointing to the start of the range.
   * @param it_end Iterator pointing to the end of the range.
   * @param number_of_parameters Negative values are associated to 0.
   */
  template <class Iterator, class = std::enable_if_t<!std::is_arithmetic_v<Iterator> > >
  Multi_parameter_filtration(Iterator it_begin, Iterator it_end, int number_of_parameters)
      : generators_(it_begin, it_end),
        generator_view_(generators_.data(),
                        number_of_parameters <= 0 ? 0 : generators_.size() / number_of_parameters,
                        number_of_parameters)
  {
    if constexpr (Ensure1Criticality) {
      if (generator_view_.extent(0) != 1) throw std::logic_error("Multiparameter filtration value is not 1-critical.");
    }
  }

  /**
   * @brief Builds filtration value with given number of parameters and values from the given range. Lets \f$ p \f$
   * be the number of parameters. The \f$ p \f$ first elements of the range have to correspond to the first generator,
   * the \f$ p \f$ next elements to the second generator and so on... So the length of the range has to be a multiple
   * of \f$ p \f$ and the number of generators will be \f$ length / p \f$. The range is represented by
   * @ref Multi_parameter_filtration::Underlying_container "" and copied into the underlying container of the class.
   *
   * @param generators Values.
   * @param number_of_parameters Negative values are associated to 0.
   */
  Multi_parameter_filtration(const Underlying_container &generators, int number_of_parameters)
      : generators_(generators),
        generator_view_(generators_.data(),
                        number_of_parameters <= 0 ? 0 : generators_.size() / number_of_parameters,
                        number_of_parameters)
  {
    GUDHI_CHECK(number_of_parameters > 0 || generators_.empty(),
                std::invalid_argument("Number of parameters cannot be 0 if the container is not empty."));

    if constexpr (Ensure1Criticality) {
      if (generator_view_.extent(0) != 1) throw std::logic_error("Multiparameter filtration value is not 1-critical.");
    }
  }

  /**
   * @brief Builds filtration value with given number of parameters and values from the given range. Lets \f$ p \f$
   * be the number of parameters. The \f$ p \f$ first elements of the range have to correspond to the first generator,
   * the \f$ p \f$ next elements to the second generator and so on... So the length of the range has to be a multiple
   * of \f$ p \f$ and the number of generators will be \f$ length / p \f$. The range is represented by
   * @ref Multi_parameter_filtration::Underlying_container "" and **moved** into the underlying container of the class.
   *
   * @param generators Values to move.
   * @param number_of_parameters Negative values are associated to 0.
   */
  Multi_parameter_filtration(Underlying_container &&generators, int number_of_parameters)
      : generators_(std::move(generators)),
        generator_view_(generators_.data(),
                        number_of_parameters <= 0 ? 0 : generators_.size() / number_of_parameters,
                        number_of_parameters)
  {
    GUDHI_CHECK(number_of_parameters > 0 || generators_.empty(),
                std::invalid_argument("Number of parameters cannot be 0 if the container is not empty."));

    if constexpr (Ensure1Criticality) {
      if (generator_view_.extent(0) != 1) throw std::logic_error("Multiparameter filtration value is not 1-critical.");
    }
  }

  /**
   * @brief Copy constructor.
   */
  Multi_parameter_filtration(const Multi_parameter_filtration &other)
      : generators_(other.generators_),
        generator_view_(generators_.data(), other.num_generators(), other.num_parameters())
  {}

  /**
   * @brief Copy constructor.
   *
   * @tparam U Type convertible into `T`.
   */
  template <typename U, bool OtherCo, bool OtherEnsure1Criticality>
  Multi_parameter_filtration(const Multi_parameter_filtration<U, OtherCo, OtherEnsure1Criticality> &other)
      : generators_(other.begin(), other.end()),
        generator_view_(generators_.data(), other.num_generators(), other.num_parameters())
  {
    if constexpr (Ensure1Criticality && !OtherEnsure1Criticality) {
      if (generator_view_.extent(0) != 1) throw std::logic_error("Multiparameter filtration value is not 1-critical.");
    }
  }

  /**
   * @brief Move constructor.
   */
  Multi_parameter_filtration(Multi_parameter_filtration &&other) noexcept
      : generators_(std::move(other.generators_)),
        generator_view_(generators_.data(), other.num_generators(), other.num_parameters())
  {}

  ~Multi_parameter_filtration() = default;

  /**
   * @brief Assign operator.
   */
  Multi_parameter_filtration &operator=(const Multi_parameter_filtration &other)
  {
    generators_ = other.generators_;
    generator_view_ = Viewer(generators_.data(), other.num_generators(), other.num_parameters());
    return *this;
  }

  /**
   * @brief Assign operator.
   *
   * @tparam U Type convertible into `T`.
   */
  template <typename U, bool OtherCo, bool OtherEnsure1Criticality>
  Multi_parameter_filtration &operator=(const Multi_parameter_filtration<U, OtherCo, OtherEnsure1Criticality> &other)
  {
    if constexpr (Ensure1Criticality && !OtherEnsure1Criticality) {
      if (other.num_generators() != 1) throw std::logic_error("Multiparameter filtration value is not 1-critical.");
    }
    generators_ = Underlying_container(other.begin(), other.end());
    generator_view_ = Viewer(generators_.data(), other.num_generators(), other.num_parameters());
    return *this;
  }

  /**
   * @brief Move assign operator.
   */
  Multi_parameter_filtration &operator=(Multi_parameter_filtration &&other) noexcept
  {
    generators_ = std::move(other.generators_);
    generator_view_ = Viewer(generators_.data(), other.num_generators(), other.num_parameters());
    return *this;
  }

  /**
   * @brief Swap operator.
   */
  friend void swap(Multi_parameter_filtration &f1, Multi_parameter_filtration &f2) noexcept
  {
    f1.generators_.swap(f2.generators_);
    swap(f1.generator_view_, f2.generator_view_);
  }

  // VECTOR-LIKE

  using value_type = T;                                                                 /**< Value type. */
  using size_type = typename Underlying_container::size_type;                           /**< Size type. */
  using difference_type = typename Underlying_container::difference_type;               /**< Difference type. */
  using reference = value_type &;                                                       /**< Reference type. */
  using const_reference = const value_type &;                                           /**< Const reference type. */
  using pointer = typename Underlying_container::pointer;                               /**< Pointer type. */
  using const_pointer = typename Underlying_container::const_pointer;                   /**< Const pointer type. */
  using iterator = typename Underlying_container::iterator;                             /**< Iterator type. */
  using const_iterator = typename Underlying_container::const_iterator;                 /**< Const iterator type. */
  using reverse_iterator = typename Underlying_container::reverse_iterator;             /**< Reverse iterator type. */
  using const_reverse_iterator = typename Underlying_container::const_reverse_iterator; /**< Const reverse iterator. */

  /**
   * @brief Returns reference to value of parameter `p` of generator `g`.
   */
  reference operator()(size_type g, size_type p) { return generator_view_(g, p); }

  /**
   * @brief Returns const reference to value of parameter `p` of generator `g`.
   */
  const_reference operator()(size_type g, size_type p) const { return generator_view_(g, p); }

  /**
   * @brief Let \f$ g \f$ be the first value in `indices` and \f$ p \f$ the second value.
   * Returns reference to value of parameter \f$ p \f$ of generator \f$ g \f$.
   *
   * @tparam IndexRange Range with a begin() and size() method.
   * @param indices Range with at least two elements. The first element should correspond to the generator number and
   * the second element to the parameter number.
   */
  template <class IndexRange = std::initializer_list<size_type>,
            class = std::enable_if_t<RangeTraits<IndexRange>::has_begin> >
  reference operator[](const IndexRange &indices)
  {
    GUDHI_CHECK(indices.size() == 2,
                std::invalid_argument(
                    "Exactly 2 indices allowed only: first the generator number, second the parameter number."));
    auto it = indices.begin();
    size_type g = *it;
    return generator_view_(g, *(++it));
  }

  /**
   * @brief Let \f$ g \f$ be the first value in `indices` and \f$ p \f$ the second value.
   * Returns reference to value of parameter \f$ p \f$ of generator \f$ g \f$.
   *
   * @tparam IndexRange Range with a begin() and size() method.
   * @param indices Range with at least two elements. The first element should correspond to the generator number and
   * the second element to the parameter number.
   */
  template <class IndexRange = std::initializer_list<size_type>,
            class = std::enable_if_t<RangeTraits<IndexRange>::has_begin> >
  const_reference operator[](const IndexRange &indices) const
  {
    GUDHI_CHECK(indices.size() == 2,
                std::invalid_argument(
                    "Exactly 2 indices allowed only: first the generator number, second the parameter number."));
    auto it = indices.begin();
    size_type g = *it;
    return generator_view_(g, *(++it));
  }

  /**
   * @brief Returns an iterator pointing the begining of the underlying container. The @ref num_parameters() first
   * elements corresponds to the first generator, the @ref num_parameters() next to the second and so on.
   *
   * @warning If a generator is modified and the new set of generators is not minimal or not sorted, the behaviour
   * of most methods is undefined. It is possible to call @ref simplify() after construction if there is a doubt to
   * ensure this property.
   */
  iterator begin() noexcept { return generators_.begin(); }

  /**
   * @brief Returns an iterator pointing the begining of the underlying container. The @ref num_parameters() first
   * elements corresponds to the first generator, the @ref num_parameters() next to the second and so on.
   */
  const_iterator begin() const noexcept { return generators_.begin(); }

  /**
   * @brief Returns an iterator pointing the begining of the underlying container. The @ref num_parameters() first
   * elements corresponds to the first generator, the @ref num_parameters() next to the second and so on.
   */
  const_iterator cbegin() const noexcept { return generators_.cbegin(); }

  /**
   * @brief Returns an iterator pointing the end of the underlying container.
   */
  iterator end() noexcept { return generators_.end(); }

  /**
   * @brief Returns an iterator pointing the end of the underlying container.
   */
  const_iterator end() const noexcept { return generators_.end(); }

  /**
   * @brief Returns an iterator pointing the end of the underlying container.
   */
  const_iterator cend() const noexcept { return generators_.cend(); }

  /**
   * @brief Returns a reverse iterator pointing to the first element from the back of the underlying container.
   * The @ref num_parameters() first elements corresponds to the last generator (in parameter reverse order), the
   * @ref num_parameters() next to the second to last and so on.
   *
   * @warning If a generator is modified and the new set of generators is not minimal or not sorted, the behaviour
   * of most methods is undefined. It is possible to call @ref simplify() after construction if there is a doubt to
   * ensure this property.
   */
  reverse_iterator rbegin() noexcept { return generators_.rbegin(); }

  /**
   * @brief Returns a reverse iterator pointing to the first element from the back of the underlying container.
   * The @ref num_parameters() first elements corresponds to the last generator (in parameter reverse order), the
   * @ref num_parameters() next to the second to last and so on.
   */
  const_reverse_iterator rbegin() const noexcept { return generators_.rbegin(); }

  /**
   * @brief Returns a reverse iterator pointing to the first element from the back of the underlying container.
   * The @ref num_parameters() first elements corresponds to the last generator (in parameter reverse order), the
   * @ref num_parameters() next to the second to last and so on.
   */
  const_reverse_iterator crbegin() const noexcept { return generators_.crbegin(); }

  /**
   * @brief Returns a reverse iterator pointing to the end of the reversed underlying container.
   */
  reverse_iterator rend() noexcept { return generators_.rend(); }

  /**
   * @brief Returns a reverse iterator pointing to the end of the reversed underlying container.
   */
  const_reverse_iterator rend() const noexcept { return generators_.rend(); }

  /**
   * @brief Returns a reverse iterator pointing to the end of the reversed underlying container.
   */
  const_reverse_iterator crend() const noexcept { return generators_.crend(); }

  /**
   * @brief Returns the size of the underlying container. Corresponds exactly to @ref num_entries(), but enables to use
   * the class as a classic range with a `begin`, `end` and `size` method.
   */
  size_type size() const noexcept { return generators_.size(); }

  /**
   * @brief Reserves space for the given number of generators in the underlying container. Does nothing if
   * `Ensure1Criticality` is true.
   */
  void reserve([[maybe_unused]] size_type number_of_generators)
  {
    if constexpr (Ensure1Criticality) {
      return;
    } else {
      generators_.reserve(num_parameters() * number_of_generators);
    }
  }

  // CONVERTERS

  // like numpy
  /**
   * @brief Returns a copy with entries casted into the type given as template parameter.
   *
   * @tparam U New type for the entries.
   * @tparam OCo New value for `Co`. Default value: `Co`.
   * @tparam OEns New value for `Ensure1Criticality`. Note that if `OEns` is set to true and the value is not
   * 1-critical, the method will throw. Default value: `Ensure1Criticality`.
   * @return Copy with new entry type.
   */
  template <typename U, bool OCo = Co, bool OEns = Ensure1Criticality>
  Multi_parameter_filtration<U, OCo, OEns> as_type() const
  {
    std::vector<U> out(generators_.begin(), generators_.end());
    return Multi_parameter_filtration<U, OCo, OEns>(std::move(out), num_parameters());
  }

  // ACCESS

  /**
   * @brief Returns the number of parameters in the filtration value.
   */
  size_type num_parameters() const { return generator_view_.extent(1); }

  /**
   * @brief Returns the number of generators in the filtration value, i.e. the criticality of the element.
   */
  size_type num_generators() const
  {
    if constexpr (Ensure1Criticality) {
      return 1;  // for possible optimizations? If there is none, we can just keep the other version
    } else {
      return generator_view_.extent(0);
    }
  }

  /**
   * @brief Returns the total number of values in the filtration value, that is,
   * @ref num_parameters() * @ref num_generators().
   */
  size_type num_entries() const { return generators_.size(); }

  /**
   * @brief Returns a filtration value with given number of parameters for which @ref is_plus_inf() returns `true`.
   */
  static Multi_parameter_filtration inf(int number_of_parameters)
  {
    return Multi_parameter_filtration(number_of_parameters, T_inf);
  }

  /**
   * @brief Returns a filtration value with given number of parameters for which @ref is_minus_inf() returns `true`.
   */
  static Multi_parameter_filtration minus_inf(int number_of_parameters)
  {
    return Multi_parameter_filtration(number_of_parameters, T_m_inf);
  }

  /**
   * @brief If `std::numeric_limits<T>::has_quiet_NaN` is true, returns a filtration value with given number of
   * parameters for which @ref is_nan() returns `true`. Otherwise, throws.
   */
  static Multi_parameter_filtration nan(int number_of_parameters)
  {
    if constexpr (std::numeric_limits<T>::has_quiet_NaN) {
      return Multi_parameter_filtration(number_of_parameters, std::numeric_limits<T>::quiet_NaN());
    } else {
      throw std::logic_error("No NaN value exists.");
    }
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
  [[nodiscard]] bool is_plus_inf() const
  {
    for (const T &v : generators_) {
      if (v != T_inf) return false;
    }
    return true;
  }

  /**
   * @brief Returns `true` if and only if the filtration value is considered as minus infinity.
   */
  [[nodiscard]] bool is_minus_inf() const
  {
    for (const T &v : generators_) {
      if (v != T_m_inf) return false;
    }
    return true;
  }

  /**
   * @brief Returns `true` if and only if the filtration value is considered as NaN.
   */
  [[nodiscard]] bool is_nan() const
  {
    if constexpr (std::numeric_limits<T>::has_quiet_NaN) {
      for (const auto &v : generators_) {
        if (!std::isnan(v)) return false;
      }
      return true;
    } else {
      return false;
    }
  }

  /**
   * @brief Returns `true` if and only if the filtration value is non-empty and is not considered as plus infinity,
   * minus infinity or NaN.
   */
  [[nodiscard]] bool is_finite() const
  {
    bool isInf = true, isMinusInf = true, isNan = true;
    for (const auto &v : generators_) {
      if (v != T_inf) isInf = false;
      if (v != T_m_inf) isMinusInf = false;
      if (!_is_nan(v)) isNan = false;
      if (!isInf && !isMinusInf && !isNan) return true;
    }
    return false;
  }

  // COMPARAISON OPERATORS

  /**
   * @brief Returns `true` if and only if the first argument is lexicographically strictly less than the second
   * argument. The "words" considered for the lexicographical order are all the generators concatenated together
   * in order of generator index and then in order of parameter index. Different from @ref operator< "", this order
   * is total.
   *
   * @tparam inverse If true, the parameter index and generator index order is inverted.
   */
  template <bool inverse = false>
  friend bool is_strict_less_than_lexicographically(const Multi_parameter_filtration &a,
                                                    const Multi_parameter_filtration &b)
  {
    if (&a == &b) return false;

    GUDHI_CHECK(a.num_parameters() == b.num_parameters(),
                std::invalid_argument("Only filtration values with same number of parameters can be compared."));

    for (std::size_t i = 0U; i < a.num_parameters() * std::min(a.num_generators(), b.num_generators()); ++i) {
      std::size_t iA = i;
      std::size_t iB = i;
      if constexpr (inverse) {
        iA = a.generators_.size() - 1 - i;
        iB = b.generators_.size() - 1 - i;
      }
      if (_is_nan(a.generators_[iA]) && !_is_nan(b.generators_[iB])) return false;
      if (_is_nan(b.generators_[iB])) return true;
      if (a.generators_[iA] < b.generators_[iB]) return true;
      if (b.generators_[iB] < a.generators_[iA]) return false;
    }
    return a.num_generators() < b.num_generators();
  }

  /**
   * @brief Returns `true` if and only if the first argument is lexicographically less than or equal to the second
   * argument. The "words" considered for the lexicographical order are all the generators concatenated together
   * in order of generator index and then in order of parameter index. Different from @ref operator<= "", this order
   * is total.
   *
   * @tparam inverse If true, the parameter index and generator index order is inverted.
   */
  template <bool inverse = false>
  friend bool is_less_or_equal_than_lexicographically(const Multi_parameter_filtration &a,
                                                      const Multi_parameter_filtration &b)
  {
    if (&a == &b) return true;

    GUDHI_CHECK(a.num_parameters() == b.num_parameters(),
                std::invalid_argument("Only filtration values with same number of parameters can be compared."));

    for (std::size_t i = 0U; i < a.num_parameters() * std::min(a.num_generators(), b.num_generators()); ++i) {
      std::size_t iA = i;
      std::size_t iB = i;
      if constexpr (inverse) {
        iA = a.generators_.size() - 1 - i;
        iB = b.generators_.size() - 1 - i;
      }
      if (_is_nan(a.generators_[iA]) && !_is_nan(b.generators_[iB])) return false;
      if (_is_nan(b.generators_[iB])) return true;
      if (a.generators_[iA] < b.generators_[iB]) return true;
      if (b.generators_[iB] < a.generators_[iA]) return false;
    }
    return a.num_generators() <= b.num_generators();
  }

  /**
   * @brief Returns `true` if and only if the cones generated by @p b are strictly contained in the
   * cones generated by @p a (recall that the cones are positive if `Co` is false and negative if `Co` is true).
   * Both @p a and @p b  have to have the same number of parameters.
   *
   * Note that not all filtration values are comparable. That is, \f$ a < b \f$ and \f$ b < a \f$ returning both false
   * does **not** imply \f$ a == b \f$. If a total order is needed, use @ref is_strict_less_than_lexicographically
   * instead.
   */
  friend bool operator<(const Multi_parameter_filtration &a, const Multi_parameter_filtration &b)
  {
    if (&a == &b) return false;

    GUDHI_CHECK(a.num_parameters() == b.num_parameters(),
                std::invalid_argument("Only filtration values with same number of parameters can be compared."));

    if (a.num_generators() == 0 || b.num_generators() == 0) return false;

    const auto &view_a = a.generator_view_;
    const auto &view_b = b.generator_view_;
    // TODO: verify if this really makes a differences in the 1-critical case, otherwise just keep the general case
    if constexpr (Ensure1Criticality) {
      if (_first_dominates(view_a, 0, view_b, 0)) return false;
      return _strictly_contains(view_a, 0, view_b, 0);
    } else {
      for (std::size_t i = 0U; i < b.num_generators(); ++i) {
        // for each generator in b, verify if it is strictly in the cone of at least one generator of a
        bool isContained = false;
        for (std::size_t j = 0U; j < a.num_generators() && !isContained; ++j) {
          // lexicographical order, so if a[j][0] dom b[j][0], than a[j'] can never strictly contain b[i] for all
          // j' > j.
          if (_first_dominates(view_a, j, view_b, i)) return false;
          isContained = _strictly_contains(view_a, j, view_b, i);
        }
        if (!isContained) return false;
      }
      return true;
    }
  }

  /**
   * @brief Returns `true` if and only if the cones generated by @p a are strictly contained in the
   * cones generated by @p b (recall that the cones are positive if `Co` is false and negative if `Co` is true).
   * Both @p a and @p b  have to have the same number of parameters.
   *
   * Note that not all filtration values are comparable. That is, \f$ a \le b \f$ and \f$ b \le a \f$ can both return
   * `false`. If a total order is needed, use @ref is_less_or_equal_than_lexicographically instead.
   */
  friend bool operator<=(const Multi_parameter_filtration &a, const Multi_parameter_filtration &b)
  {
    GUDHI_CHECK(a.num_parameters() == b.num_parameters(),
                std::invalid_argument("Only filtration values with same number of parameters can be compared."));

    if (a.num_generators() == 0 || b.num_generators() == 0) return false;
    if (a.is_nan() || b.is_nan()) return false;
    if (&a == &b) return true;

    const auto &view_a = a.generator_view_;
    const auto &view_b = b.generator_view_;
    // TODO: verify if this really makes a differences in the 1-critical case, otherwise just keep the general case
    if constexpr (Ensure1Criticality) {
      if (_first_strictly_dominates(view_a, 0, view_b, 0)) return false;
      return _contains(view_a, 0, view_b, 0);
    } else {
      // check if this curves is below other's curve
      //  ie for each guy in this, check if there is a guy in other that dominates him
      for (std::size_t i = 0U; i < b.num_generators(); ++i) {
        // for each generator in b, verify if it is in the cone of at least one generator of a
        bool isContained = false;
        for (std::size_t j = 0U; j < a.num_generators() && !isContained; ++j) {
          // lexicographical order, so if a[j][0] strictly dom b[j][0], than a[j'] can never contain b[i] for all
          // j' > j.
          if (_first_strictly_dominates(view_a, j, view_b, i)) return false;
          isContained = _contains(view_a, j, view_b, i);
        }
        if (!isContained) return false;
      }
      return true;
    }
  }

  /**
   * @brief Returns `true` if and only if the cones generated by @p b are contained in or are (partially)
   * equal to the cones generated by @p a (recall that the cones are positive if `Co` is false and negative if `Co` is
   * true).
   * Both @p a and @p b  have to have the same number of parameters.
   *
   * Note that not all filtration values are comparable. That is, \f$ a > b \f$ and \f$ b > a \f$ returning both false
   * does **not** imply \f$ a == b \f$. If a total order is needed, use @ref is_strict_less_than_lexicographically
   * instead.
   */
  friend bool operator>(const Multi_parameter_filtration &a, const Multi_parameter_filtration &b) { return b < a; }

  /**
   * @brief Returns `true` if and only if the cones generated by @p a are contained in or are (partially)
   * equal to the cones generated by @p b (recall that the cones are positive if `Co` is false and negative if `Co` is
   * true).
   * Both @p a and @p b  have to have the same number of parameters.
   *
   * Note that not all filtration values are comparable. That is, \f$ a \ge b \f$ and \f$ b \ge a \f$ can both return
   * `false`. If a total order is needed, use @ref is_less_or_equal_than_lexicographically instead.
   */
  friend bool operator>=(const Multi_parameter_filtration &a, const Multi_parameter_filtration &b) { return b <= a; }

  /**
   * @brief Returns `true` if and only if for each \f$ i,j \f$, \f$ a(i,j) \f$ is equal to \f$ b(i,j) \f$.
   */
  friend bool operator==(const Multi_parameter_filtration &a, const Multi_parameter_filtration &b)
  {
    if (&a == &b) return true;
    // assumes lexicographical order for both
    return a.generators_ == b.generators_;
  }

  /**
   * @brief Returns `true` if and only if \f$ a == b \f$ returns `false`.
   */
  friend bool operator!=(const Multi_parameter_filtration &a, const Multi_parameter_filtration &b) { return !(a == b); }

  // ARITHMETIC OPERATORS

  // opposite
  /**
   * @brief Returns a filtration value such that an entry at index \f$ i,j \f$ is equal to \f$ -f(i,j) \f$.
   *
   * Used conventions:
   * - \f$ -NaN = NaN \f$.
   *
   * @param f Value to opposite.
   * @return The opposite of @p f.
   */
  friend Multi_parameter_filtration operator-(const Multi_parameter_filtration &f)
  {
    using F = Multi_parameter_filtration;

    Underlying_container result(f.generators_);
    std::for_each(result.begin(), result.end(), [](T &v) {
      if (v == F::T_inf)
        v = F::T_m_inf;
      else if (v == F::T_m_inf)
        v = F::T_inf;
      else
        v = -v;
    });
    return Multi_parameter_filtration(std::move(result), f.num_parameters());
  }

  // subtraction
  /**
   * @brief Returns a filtration value such that an entry at index \f$ (g,p) \f$, with
   * \f$ 0 \leq g \leq num_generators \f$ and \f$ 0 \leq p \leq num_parameters \f$ is equal to \f$ f(g,p) - r(p) \f$
   * if \f$ p < length_r \f$ and to \f$ f(g,p) \f$ otherwise.
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
  template <class ValueRange, class = std::enable_if_t<RangeTraits<ValueRange>::has_begin> >
  friend Multi_parameter_filtration operator-(Multi_parameter_filtration f, const ValueRange &r)
  {
    f -= r;
    return f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ (g,p) \f$, with
   * \f$ 0 \leq g \leq num_generators \f$ and \f$ 0 \leq p \leq num_parameters \f$ is equal to \f$ r(p) - f(g,p) \f$
   * if \f$ p < length_r \f$ and to \f$ -f(g,p) \f$ otherwise.
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
  template <class ValueRange,
            class = std::enable_if_t<RangeTraits<ValueRange>::has_begin &&
                                     !std::is_same_v<ValueRange, Multi_parameter_filtration> > >
  friend Multi_parameter_filtration operator-(const ValueRange &r, Multi_parameter_filtration f)
  {
    f._apply_operation(r, [](T &valF, const T &valR) {
      valF = -valF;
      _add(valF, valR);
    });
    return f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ (g,p) \f$ is equal to \f$ f(g,p) - val \f$.
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
  friend Multi_parameter_filtration operator-(Multi_parameter_filtration f, const T &val)
  {
    f -= val;
    return f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ (g,p) \f$ is equal to \f$ val - f(g,p) \f$.
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
  friend Multi_parameter_filtration operator-(const T &val, Multi_parameter_filtration f)
  {
    f._apply_operation(val, [](T &valF, const T &valR) {
      valF = -valF;
      _add(valF, valR);
    });
    return f;
  }

  /**
   * @brief Modifies the first parameter such that an entry at index \f$ (g,p) \f$, with
   * \f$ 0 \leq g \leq num_generators \f$ and \f$ 0 \leq p \leq num_parameters \f$ is equal to \f$ f(g,p) - r(p) \f$
   * if \f$ p < length_r \f$ and to \f$ f(g,p) \f$ otherwise.
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
  template <class ValueRange, class = std::enable_if_t<RangeTraits<ValueRange>::has_begin> >
  friend Multi_parameter_filtration &operator-=(Multi_parameter_filtration &f, const ValueRange &r)
  {
    f._apply_operation(r, [](T &valF, const T &valR) { _subtract(valF, valR); });
    return f;
  }

  /**
   * @brief Modifies the first parameter such that an entry at index \f$ (g,p) \f$ is equal to \f$ f(g,p) - val \f$.
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
  friend Multi_parameter_filtration &operator-=(Multi_parameter_filtration &f, const T &val)
  {
    f._apply_operation(val, [](T &valF, const T &valR) { _subtract(valF, valR); });
    return f;
  }

  // addition
  /**
   * @brief Returns a filtration value such that an entry at index \f$ (g,p) \f$, with
   * \f$ 0 \leq g \leq num_generators \f$ and \f$ 0 \leq p \leq num_parameters \f$ is equal to \f$ f(g,p) + r(p) \f$
   * if \f$ p < length_r \f$ and to \f$ f(g,p) \f$ otherwise.
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
  template <class ValueRange, class = std::enable_if_t<RangeTraits<ValueRange>::has_begin> >
  friend Multi_parameter_filtration operator+(Multi_parameter_filtration f, const ValueRange &r)
  {
    f += r;
    return f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ (g,p) \f$, with
   * \f$ 0 \leq g \leq num_generators \f$ and \f$ 0 \leq p \leq num_parameters \f$ is equal to \f$ r(p) + f(g,p) \f$
   * if \f$ p < length_r \f$ and to \f$ f(g,p) \f$ otherwise.
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
  template <class ValueRange,
            class = std::enable_if_t<RangeTraits<ValueRange>::has_begin &&
                                     !std::is_same_v<ValueRange, Multi_parameter_filtration> > >
  friend Multi_parameter_filtration operator+(const ValueRange &r, Multi_parameter_filtration f)
  {
    f += r;
    return f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ (g,p) \f$ is equal to \f$ f(g,p) + val \f$.
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
  friend Multi_parameter_filtration operator+(Multi_parameter_filtration f, const T &val)
  {
    f += val;
    return f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ (g,p) \f$ is equal to \f$ val + f(g,p) \f$.
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
  friend Multi_parameter_filtration operator+(const T &val, Multi_parameter_filtration f)
  {
    f += val;
    return f;
  }

  /**
   * @brief Modifies the first parameter such that an entry at index \f$ (g,p) \f$, with
   * \f$ 0 \leq g \leq num_generators \f$ and \f$ 0 \leq p \leq num_parameters \f$ is equal to \f$ f(g,p) + r(p) \f$
   * if \f$ p < length_r \f$ and to \f$ f(g,p) \f$ otherwise.
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
  template <class ValueRange, class = std::enable_if_t<RangeTraits<ValueRange>::has_begin> >
  friend Multi_parameter_filtration &operator+=(Multi_parameter_filtration &f, const ValueRange &r)
  {
    f._apply_operation(r, [](T &valF, const T &valR) { _add(valF, valR); });
    return f;
  }

  /**
   * @brief Modifies the first parameter such that an entry at index \f$ (g,p) \f$ is equal to \f$ f(g,p) + val \f$.
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
  friend Multi_parameter_filtration &operator+=(Multi_parameter_filtration &f, const T &val)
  {
    f._apply_operation(val, [](T &valF, const T &valR) { _add(valF, valR); });
    return f;
  }

  // multiplication
  /**
   * @brief Returns a filtration value such that an entry at index \f$ (g,p) \f$, with
   * \f$ 0 \leq g \leq num_generators \f$ and \f$ 0 \leq p \leq num_parameters \f$ is equal to \f$ f(g,p) * r(p) \f$
   * if \f$ p < length_r \f$ and to \f$ f(g,p) \f$ otherwise.
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
  template <class ValueRange, class = std::enable_if_t<RangeTraits<ValueRange>::has_begin> >
  friend Multi_parameter_filtration operator*(Multi_parameter_filtration f, const ValueRange &r)
  {
    f *= r;
    return f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ (g,p) \f$, with
   * \f$ 0 \leq g \leq num_generators \f$ and \f$ 0 \leq p \leq num_parameters \f$ is equal to \f$ r(p) * f(g,p) \f$
   * if \f$ p < length_r \f$ and to \f$ f(g,p) \f$ otherwise.
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
  template <class ValueRange,
            class = std::enable_if_t<RangeTraits<ValueRange>::has_begin &&
                                     !std::is_same_v<ValueRange, Multi_parameter_filtration> > >
  friend Multi_parameter_filtration operator*(const ValueRange &r, Multi_parameter_filtration f)
  {
    f *= r;
    return f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ (g,p) \f$ is equal to \f$ f(g,p) * val \f$.
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
  friend Multi_parameter_filtration operator*(Multi_parameter_filtration f, const T &val)
  {
    f *= val;
    return f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ (g,p) \f$ is equal to \f$ val * f(g,p) \f$.
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
  friend Multi_parameter_filtration operator*(const T &val, Multi_parameter_filtration f)
  {
    f *= val;
    return f;
  }

  /**
   * @brief Modifies the first parameter such that an entry at index \f$ (g,p) \f$, with
   * \f$ 0 \leq g \leq num_generators \f$ and \f$ 0 \leq p \leq num_parameters \f$ is equal to \f$ f(g,p) * r(p) \f$
   * if \f$ p < length_r \f$ and to \f$ f(g,p) \f$ otherwise.
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
  template <class ValueRange, class = std::enable_if_t<RangeTraits<ValueRange>::has_begin> >
  friend Multi_parameter_filtration &operator*=(Multi_parameter_filtration &f, const ValueRange &r)
  {
    f._apply_operation(r, [](T &valF, const T &valR) { _multiply(valF, valR); });
    return f;
  }

  /**
   * @brief Modifies the first parameter such that an entry at index \f$ (g,p) \f$ is equal to \f$ f(g,p) * val \f$.
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
  friend Multi_parameter_filtration &operator*=(Multi_parameter_filtration &f, const T &val)
  {
    f._apply_operation(val, [](T &valF, const T &valR) { _multiply(valF, valR); });
    return f;
  }

  // division
  /**
   * @brief Returns a filtration value such that an entry at index \f$ (g,p) \f$, with
   * \f$ 0 \leq g \leq num_generators \f$ and \f$ 0 \leq p \leq num_parameters \f$ is equal to \f$ f(g,p) / r(p) \f$
   * if \f$ p < length_r \f$ and to \f$ f(g,p) \f$ otherwise.
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
  template <class ValueRange, class = std::enable_if_t<RangeTraits<ValueRange>::has_begin> >
  friend Multi_parameter_filtration operator/(Multi_parameter_filtration f, const ValueRange &r)
  {
    f /= r;
    return f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ (g,p) \f$, with
   * \f$ 0 \leq g \leq num_generators \f$ and \f$ 0 \leq p \leq num_parameters \f$ is equal to \f$ r(p) / f(g,p) \f$
   * if \f$ p < length_r \f$ and to \f$ f(g,p) \f$ otherwise.
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
  template <class ValueRange,
            class = std::enable_if_t<RangeTraits<ValueRange>::has_begin &&
                                     !std::is_same_v<ValueRange, Multi_parameter_filtration> > >
  friend Multi_parameter_filtration operator/(const ValueRange &r, Multi_parameter_filtration f)
  {
    f._apply_operation(r, [](T &valF, const T &valR) {
      T tmp = valF;
      valF = valR;
      _divide(valF, tmp);
    });
    return f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ (g,p) \f$ is equal to \f$ f(g,p) / val \f$.
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
  friend Multi_parameter_filtration operator/(Multi_parameter_filtration f, const T &val)
  {
    f /= val;
    return f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ (g,p) \f$ is equal to \f$ val / f(g,p) \f$.
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
  friend Multi_parameter_filtration operator/(const T &val, Multi_parameter_filtration f)
  {
    f._apply_operation(val, [](T &valF, const T &valR) {
      T tmp = valF;
      valF = valR;
      _divide(valF, tmp);
    });
    return f;
  }

  /**
   * @brief Modifies the first parameter such that an entry at index \f$ (g,p) \f$, with
   * \f$ 0 \leq g \leq num_generators \f$ and \f$ 0 \leq p \leq num_parameters \f$ is equal to \f$ f(g,p) / r(p) \f$
   * if \f$ p < length_r \f$ and to \f$ f(g,p) \f$ otherwise.
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
  template <class ValueRange, class = std::enable_if_t<RangeTraits<ValueRange>::has_begin> >
  friend Multi_parameter_filtration &operator/=(Multi_parameter_filtration &f, const ValueRange &r)
  {
    f._apply_operation(r, [](T &valF, const T &valR) { _divide(valF, valR); });
    return f;
  }

  /**
   * @brief Modifies the first parameter such that an entry at index \f$ (g,p) \f$ is equal to \f$ f(g,p) / val \f$.
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
  friend Multi_parameter_filtration &operator/=(Multi_parameter_filtration &f, const T &val)
  {
    f._apply_operation(val, [](T &valF, const T &valR) { _divide(valF, valR); });
    return f;
  }

  // MODIFIERS

  /**
   * @brief Sets the number of generators. If there were less generators before, new generators with default values are
   * constructed. If there were more generators before, the exceed of generators is destroyed (any generator with index
   * higher or equal to @p g to be more precise). If @p g is zero, the methods does nothing.
   *
   * Fails to compile if `Ensure1Criticality` is true.
   *
   * @warning All new generators will be set to infinity (`Co` is true) or -infinity (`Co` is false). That is, the new
   * filtration value is not minimal anymore. Make sure to fill them with real generators or to remove them before
   * using other methods.
   *
   * @warning Be sure to call @ref simplify if necessary after initializing all the generators. Most methods will have
   * an undefined behaviour if the set of generators is not minimal or sorted.
   *
   * @param g New number of generators.
   */
  void set_num_generators(size_type g)
  {
    static_assert(!Ensure1Criticality, "Number of generators cannot be set for a 1-critical only filtration value.");

    if (g == 0) return;
    generators_.resize(g * num_parameters(), _get_default_value());
    generator_view_.update_data(generators_.data());  // in case it was relocated
    generator_view_.update_extent(0, g);
  }

  /**
   * @brief Adds the given generator to the filtration value such that the set remains minimal and sorted.
   * It is therefore possible that the generator is ignored if it does not generated any new lifetime or that
   * old generators disappear if they are overshadowed by the new one.
   *
   * @tparam GeneratorRange Range of elements convertible to `T`. Must have a begin(), end() method and the iterator
   * type should satisfy the requirements of the standard `LegacyForwardIterator`.
   * @param x New generator to add. Has to have the same number of parameters than @ref num_parameters().
   * @return true If and only if the generator is actually added to the set of generators.
   * @return false Otherwise.
   */
  template <class GeneratorRange = std::initializer_list<T>,
            class = std::enable_if_t<RangeTraits<GeneratorRange>::has_begin> >
  bool add_generator(const GeneratorRange &x)
  {
    return add_generator(x.begin(), x.end());
  }

  /**
   * @brief Adds the given generator to the filtration value such that the set remains minimal and sorted.
   * It is therefore possible that the generator is ignored if it does not generated any new lifetime or that
   * old generators disappear if they are overshadowed by the new one.
   *
   * @tparam Iterator Iterator class satisfying the requirements of the standard `LegacyForwardIterator`.
   * The dereferenced type has to be convertible to `T`.
   * @param genStart Iterator pointing to the begining of the range.
   * @param genEnd Iterator pointing to the end of the range.
   * @return true If and only if the generator is actually added to the set of generators.
   * @return false Otherwise.
   */
  template <class Iterator>
  bool add_generator(Iterator genStart, Iterator genEnd)
  {
    GUDHI_CHECK(std::distance(genStart, genEnd) == static_cast<int>(num_parameters()),
                std::invalid_argument("Wrong range size. Should correspond to the number of parameters."));

    const int newIndex = -1;

    std::size_t end = num_generators();
    std::vector<int> indices(end);
    std::iota(indices.begin(), indices.end(), 0);

    if (_generator_can_be_added(genStart, 0, end, indices)) {
      indices.resize(end);
      indices.push_back(newIndex);
      _build_from(indices, newIndex, genStart, genEnd);
      if constexpr (Ensure1Criticality) {
        if (generator_view_.extent(0) != 1)
          throw std::logic_error("Multiparameter filtration value is not 1-critical anymore.");
      }
      return true;
    }

    return false;
  }

  /**
   * @brief Adds the given generator to the filtration value without any verifications or simplifications at the end
   * of the set.
   *
   * Fails to compile if `Ensure1Criticality` is true.
   *
   * @warning If the resulting set of generators is not minimal or sorted after modification, some methods will have an
   * undefined behaviour. Be sure to call @ref simplify() before using them.
   *
   * @tparam GeneratorRange Range of elements convertible to `T`. Must have a begin(), end() and size() method.
   * @param x New generator to add. Must have the same number of parameters than @ref num_parameters().
   */
  template <class GeneratorRange = std::initializer_list<T>,
            class = std::enable_if_t<RangeTraits<GeneratorRange>::has_begin> >
  void add_guaranteed_generator(const GeneratorRange &x)
  {
    static_assert(!Ensure1Criticality, "Cannot add additional generator to a 1-critical only filtration value.");

    GUDHI_CHECK(x.size() == num_parameters(),
                std::invalid_argument("Wrong range size. Should correspond to the number of parameters."));

    generators_.insert(generators_.end(), x.begin(), x.end());
    generator_view_.update_data(generators_.data());  // in case it was relocated
    generator_view_.update_extent(0, num_generators() + 1);
  }

  /**
   * @brief Simplifies the current set of generators such that it becomes minimal. Also orders it in increasing
   * lexicographical order. Only necessary if generators were added "by hand" without verification either trough the
   * constructor or with @ref add_guaranteed_generator "", etc.
   */
  void simplify()
  {
    if constexpr (Ensure1Criticality) {
      return;
    } else {
      std::size_t end = 0;
      std::vector<int> indices(num_generators());
      std::iota(indices.begin(), indices.end(), 0);

      for (std::size_t curr = 0; curr < num_generators(); ++curr) {
        if (_generator_can_be_added(
                generators_.begin() + generator_view_.mapping()(indices[curr], 0), 0, end, indices)) {
          std::swap(indices[end], indices[curr]);
          ++end;
        }
      }

      indices.resize(end);
      _build_from(indices);
    }
  }

  /**
   * @brief Removes all empty generators from the filtration value. If @p include_infinities is true, it also
   * removes the generators at infinity or minus infinity. As empty generators are not possible (except if number of
   * parameters is 0), the method does nothing except sorting the set of generators if @p include_infinities is false.
   * Exists mostly for interface purposes.
   * If the set of generators is empty after removals, it is set to minus infinity if `Co` is false or to infinity
   * if `Co` is true.
   *
   * @warning If the resulting set of generators is not minimal after the removals/sorting, some methods will have an
   * undefined behaviour. Be sure to call @ref simplify() before using them.
   *
   * @param include_infinities If true, removes also infinity values.
   */
  void remove_empty_generators(bool include_infinities = false)
  {
    // TODO: verify if this really makes a differences in the 1-critical case, otherwise just keep the general case
    if constexpr (Ensure1Criticality) {
      if (!include_infinities) return;
      bool allNaN = true, allInf = true;
      bool toEmpty = true;
      for (size_type p = 0; p < num_parameters() && toEmpty; ++p) {
        if (!_is_nan(generators_[p])) allNaN = false;
        if constexpr (Co) {
          if (generators_[p] != T_m_inf) allInf = false;
        } else {
          if (generators_[p] != T_inf) allInf = false;
        }
        if (!allNaN && !allInf) toEmpty = false;
      }
      if (toEmpty) generators_.clear();
    } else {
      std::vector<int> indices;
      indices.reserve(num_generators());
      for (int i = 0; i < static_cast<int>(num_generators()); ++i) {
        if (!include_infinities || _is_finite(i)) indices.push_back(i);
      }
      _build_from(indices);  // sorts
    }

    if (generators_.empty()) {
      generators_.resize(num_parameters(), _get_default_value());
      generator_view_.update_data(generators_.data());  // in case it was relocated
      generator_view_.update_extent(0, 1);
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
   * @param exclude_infinite_values If true, values at infinity or minus infinity are not affected.
   * @return true If the filtration value was actually modified.
   * @return false Otherwise.
   */
  template <class GeneratorRange = std::initializer_list<value_type>,
            class = std::enable_if_t<RangeTraits<GeneratorRange>::has_begin> >
  bool push_to_least_common_upper_bound(const GeneratorRange &x, bool exclude_infinite_values = false)
  {
    GUDHI_CHECK(x.size() == num_parameters(),
                std::invalid_argument("Wrong range size. Should correspond to the number of parameters."));

    bool xIsInf = true, xIsMinusInf = true, xIsNaN = true;
    bool thisIsInf = true, thisIsMinusInf = true, thisIsNaN = true;

    // if one is not finite, we can avoid the heavy simplification process
    _get_infinity_statuses(generator_view_, x, thisIsInf, thisIsMinusInf, thisIsNaN, xIsInf, xIsMinusInf, xIsNaN);

    if (thisIsInf || thisIsNaN || xIsNaN || xIsMinusInf || (xIsInf && exclude_infinite_values)) return false;

    if (xIsInf || thisIsMinusInf) {
      generators_ = Underlying_container(x.begin(), x.end());
      generator_view_.update_data(generators_.data());  // in case it was relocated
      generator_view_.update_extent(0, 1);
      return true;
    }

    bool modified = false;

    auto push_generator_value = [&](T &val, T valX) {
      if (exclude_infinite_values && (valX == T_inf || valX == T_m_inf)) return;
      modified |= val < valX;
      val = valX > val ? valX : val;
    };

    // TODO: verify if this really makes a differences in the 1-critical case, otherwise just keep the general case
    if constexpr (Ensure1Criticality) {
      auto it = x.begin();
      for (size_type p = 0; p < num_parameters(); ++p) {
        push_generator_value(generators_[p], *it);
        ++it;
      }
    } else {
      for (size_type g = 0; g < num_generators(); ++g) {
        auto it = x.begin();
        for (size_type p = 0; p < num_parameters(); ++p) {
          push_generator_value(generator_view_(g, p), *it);
          ++it;
        }
      }

      if (modified && num_generators() > 1) simplify();
    }

    return modified;
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
            class = std::enable_if_t<RangeTraits<GeneratorRange>::has_begin> >
  bool pull_to_greatest_common_lower_bound(const GeneratorRange &x, bool exclude_infinite_values = false)
  {
    GUDHI_CHECK(x.size() == num_parameters(),
                std::invalid_argument("Wrong range size. Should correspond to the number of parameters."));

    bool xIsInf = true, xIsMinusInf = true, xIsNaN = true;
    bool thisIsInf = true, thisIsMinusInf = true, thisIsNaN = true;

    // if one is not finite, we can avoid the heavy simplification process
    _get_infinity_statuses(generator_view_, x, thisIsInf, thisIsMinusInf, thisIsNaN, xIsInf, xIsMinusInf, xIsNaN);

    if (xIsInf || thisIsNaN || xIsNaN || thisIsMinusInf || (xIsMinusInf && exclude_infinite_values)) return false;

    if (thisIsInf || xIsMinusInf) {
      generators_ = Underlying_container(x.begin(), x.end());
      generator_view_.update_data(generators_.data());  // in case it was relocated
      generator_view_.update_extent(0, 1);
      return true;
    }

    bool modified = false;

    auto pull_generator_value = [&](T &val, T valX) {
      if (exclude_infinite_values && (valX == T_inf || valX == T_m_inf)) return;
      modified |= val > valX;
      val = valX < val ? valX : val;
    };

    // TODO: verify if this really makes a differences in the 1-critical case, otherwise just keep the general case
    if constexpr (Ensure1Criticality) {
      auto it = x.begin();
      for (size_type p = 0; p < num_parameters(); ++p) {
        pull_generator_value(generators_[p], *it);
        ++it;
      }
    } else {
      for (size_type g = 0; g < num_generators(); ++g) {
        auto it = x.begin();
        for (size_type p = 0; p < num_parameters(); ++p) {
          pull_generator_value(generator_view_(g, p), *it);
          ++it;
        }
      }

      if (modified && num_generators() > 1) simplify();
    }

    return modified;
  }

  /**
   * @brief Projects the filtration value into the given grid. If @p coordinate is false, the entries are set to
   * the nearest upper bound value with the same parameter in the grid. Otherwise, the entries are set to the indices
   * of those nearest upper bound values.
   * The grid has to be represented as a vector of ordered ranges of values convertible into `T`. An index
   * \f$ i \f$ of the vector corresponds to the same parameter as the index \f$ i \f$ in a generator of the filtration
   * value. The ranges correspond to the possible values of the parameters, ordered by increasing value, forming
   * therefore all together a 2D grid.
   *
   * @tparam OneDimArray A range of values convertible into `T` ordered by increasing value. Has to implement
   * a begin, end and operator[] method.
   * @param grid Vector of @p OneDimArray with size at least number of filtration parameters.
   * @param coordinate If true, the values are set to the coordinates of the projection in the grid. If false,
   * the values are set to the values at the coordinates of the projection.
   */
  template <typename OneDimArray>
  void project_onto_grid(const std::vector<OneDimArray> &grid, bool coordinate = true)
  {
    GUDHI_CHECK(
        grid.size() >= num_parameters(),
        std::invalid_argument("The grid should not be smaller than the number of parameters in the filtration value."));

    auto project_generator_value = [&](T &val, const OneDimArray &filtration) {
      auto v = static_cast<typename OneDimArray::value_type>(val);
      auto d = std::distance(filtration.begin(), std::lower_bound(filtration.begin(), filtration.end(), v));
      if (d != 0 && std::abs(v - filtration[d]) > std::abs(v - filtration[d - 1])) {
        --d;
      }
      val = coordinate ? static_cast<T>(d) : static_cast<T>(filtration[d]);
    };

    // TODO: verify if this really makes a differences in the 1-critical case, otherwise just keep the general case
    if constexpr (Ensure1Criticality) {
      for (size_type p = 0; p < num_parameters(); ++p) {
        project_generator_value(generators_[p], grid[p]);
      }
    } else {
      for (size_type g = 0; g < num_generators(); ++g) {
        for (size_type p = 0; p < num_parameters(); ++p) {
          project_generator_value(generator_view_(g, p), grid[p]);
        }
      }

      if (!coordinate && num_generators() > 1) simplify();
    }
  }

  // FONCTIONNALITIES

  /**
   * @brief Returns a generator with the minimal values of all parameters in any generator of the given filtration
   * value. That is, the greatest lower bound of all generators.
   */
  friend Multi_parameter_filtration factorize_below(const Multi_parameter_filtration &f)
  {
    if (f.num_generators() <= 1) return f;

    bool nan = true;
    Underlying_container result(f.num_parameters(), T_inf);
    for (size_type p = 0; p < f.num_parameters(); ++p) {
      for (size_type g = 0; g < f.num_generators(); ++g) {
        T val = f(g, p);
        if (!_is_nan(val)) {
          nan = false;
          result[p] = val < result[p] ? val : result[p];
        }
      }
      if (nan)
        result[p] = std::numeric_limits<T>::quiet_NaN();
      else
        nan = true;
    }
    return Multi_parameter_filtration(std::move(result), f.num_parameters());
  }

  /**
   * @brief Returns a generator with the maximal values of all parameters in any generator of the given filtration
   * value. That is, the least upper bound of all generators.
   */
  friend Multi_parameter_filtration factorize_above(const Multi_parameter_filtration &f)
  {
    if (f.num_generators() <= 1) return f;

    bool nan = true;
    Underlying_container result(f.num_parameters(), T_m_inf);
    for (size_type p = 0; p < f.num_parameters(); ++p) {
      for (size_type g = 0; g < f.num_generators(); ++g) {
        T val = f(g, p);
        if (!_is_nan(val)) {
          nan = false;
          result[p] = val > result[p] ? val : result[p];
        }
      }
      if (nan)
        result[p] = std::numeric_limits<T>::quiet_NaN();
      else
        nan = true;
    }
    return Multi_parameter_filtration(std::move(result), f.num_parameters());
  }

  /**
   * @brief Computes the smallest (resp. the greatest if `Co` is true) scalar product of the all generators with the
   * given vector.
   *
   * @tparam U Arithmetic type of the result. Default value: `T`.
   * @param f Filtration value.
   * @param x Vector of coefficients.
   * @return Scalar product of @p f with @p x.
   */
  template <typename U = T>
  friend U compute_linear_projection(const Multi_parameter_filtration &f, const std::vector<U> &x)
  {
    auto project_generator = [&](size_type g) -> U {
      U projection = 0;
      std::size_t size = std::min(x.size(), f.num_parameters());
      for (std::size_t i = 0; i < size; i++) projection += x[i] * static_cast<U>(f(g, i));
      return projection;
    };

    if (f.num_generators() == 1) return project_generator(0);

    if constexpr (Co) {
      U projection = std::numeric_limits<U>::lowest();
      for (size_type g = 0; g < f.num_generators(); ++g) {
        // Order in the max important to spread possible NaNs
        projection = std::max(project_generator(g), projection);
      }
      return projection;
    } else {
      U projection = std::numeric_limits<U>::max();
      for (size_type g = 0; g < f.num_generators(); ++g) {
        // Order in the min important to spread possible NaNs
        projection = std::min(project_generator(g), projection);
      }
      return projection;
    }
  }

  /**
   * @brief Computes the euclidean distance from the first parameter to the second parameter as the minimum of
   * all Euclidean distances between a generator of @p f and a generator of @p other.
   *
   * @param f Source filtration value.
   * @param other Target filtration value.
   * @return Euclidean distance between @p f and @p other.
   */
  template <typename U = T>
  friend U compute_euclidean_distance_to(const Multi_parameter_filtration &f, const Multi_parameter_filtration &other)
  {
    GUDHI_CHECK(f.num_parameters() == other.num_parameters(),
                std::invalid_argument("We cannot compute the distance between two points of different dimensions."));

    // TODO: verify if this really makes a differences in the 1-critical case, otherwise just keep the general case
    if constexpr (Ensure1Criticality) {
      return _compute_frobenius_norm(f.num_parameters(),
                                     [&](size_type p) -> T { return f.generators_[p] - other.generators_[p]; });
    } else {
      U res = std::numeric_limits<U>::max();
      for (size_type g1 = 0; g1 < f.num_generators(); ++g1) {
        for (size_type g2 = 0; g2 < other.num_generators(); ++g2) {
          // Euclidean distance as a Frobenius norm with matrix 1 x p and values 'f(g1, p) - other(g2, p)'
          // Order in the min important to spread possible NaNs
          res = std::min(
              _compute_frobenius_norm(f.num_parameters(), [&](size_type p) -> T { return f(g1, p) - other(g2, p); }),
              res);
        }
      }
      return res;
    }
  }

  /**
   * @brief Computes the norm of the given filtration value.
   *
   * The filtration value is seen as a \f$ num_generators x num_parameters \f$ matrix and a standard Frobenius norm
   * is computed from it: the square root of the sum of the squares of all elements in the matrix.
   *
   * @param f Filtration value.
   * @return The norm of @p f.
   */
  template <typename U = T>
  friend U compute_norm(const Multi_parameter_filtration &f)
  {
    // Frobenius norm with matrix g x p based on Euclidean norm
    return _compute_frobenius_norm(f.num_entries(), [&](size_type i) -> T { return f.generators_[i]; });
  }

  /**
   * @brief Computes the coordinates in the given grid, corresponding to the nearest upper bounds of the entries
   * in the given filtration value.
   * The grid has to be represented as a vector of vectors of ordered values convertible into `OutValue`. An index
   * \f$ i \f$ of the vector corresponds to the same parameter as the index \f$ i \f$ in a generator of the filtration
   * value. The ranges correspond to the possible values of the parameters, ordered by increasing value, forming
   * therefore all together a 2D grid.
   *
   * @tparam OutValue Signed arithmetic type. Default value: std::int32_t.
   * @tparam U Type which is convertible into `OutValue`.
   * @param f Filtration value to project.
   * @param grid Vector of vectors to project into.
   * @return Filtration value \f$ out \f$ whose entry correspond to the indices of the projected values. That is,
   * the projection of \f$ f(g,p) \f$ is \f$ grid[p][out(g,p)] \f$.
   */
  template <typename OutValue = std::int32_t, typename U = T>
  friend Multi_parameter_filtration<OutValue, Co, Ensure1Criticality> compute_coordinates_in_grid(
      Multi_parameter_filtration f,
      const std::vector<std::vector<U> > &grid)
  {
    // TODO: by replicating the code of "project_onto_grid", this could be done with just one copy
    // instead of two. But it is not clear if it is really worth it, i.e., how much the change in type is really
    // necessary in the use cases. To see later.
    f.project_onto_grid(grid);
    if constexpr (std::is_same_v<OutValue, T>) {
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
   * @tparam U Signed arithmetic type.
   * @param f Filtration value storing coordinates compatible with `grid`.
   * @param grid Vector of vector.
   * @return Filtration value \f$ out \f$ whose entry correspond to \f$ out(g,p) = grid[p][f(g,p)] \f$.
   */
  template <typename U>
  friend Multi_parameter_filtration<U, Co, Ensure1Criticality> evaluate_coordinates_in_grid(
      const Multi_parameter_filtration &f,
      const std::vector<std::vector<U> > &grid)
  {
    GUDHI_CHECK(grid.size() >= f.num_parameters(),
                std::invalid_argument(
                    "The size of the grid should correspond to the number of parameters in the filtration value."));

    U grid_inf = Multi_parameter_filtration<U, Co, Ensure1Criticality>::T_inf;
    std::vector<U> outVec(f.num_entries());

    // TODO: verify if this really makes a differences in the 1-critical case, otherwise just keep the general case
    if constexpr (Ensure1Criticality) {
      for (size_type p = 0; p < f.num_parameters(); ++p) {
        const std::vector<U> &filtration = grid[p];
        const T &c = f.generators_[p];
        outVec[p] = (c == T_inf ? grid_inf : filtration[c]);
      }
    } else {
      for (size_type g = 0; g < f.num_generators(); ++g) {
        for (size_type p = 0; p < f.num_parameters(); ++p) {
          const std::vector<U> &filtration = grid[p];
          const T &c = f(g, p);
          outVec[f.generator_view_.mapping()(g, p)] = (c == T_inf ? grid_inf : filtration[c]);
        }
      }
    }

    Multi_parameter_filtration<U, Co, Ensure1Criticality> out(std::move(outVec), f.num_parameters());
    if constexpr (!Ensure1Criticality)
      if (out.num_generators() > 1) out.simplify();
    return out;
  }

  // UTILITIES

  /**
   * @brief Outstream operator.
   */
  friend std::ostream &operator<<(std::ostream &stream, const Multi_parameter_filtration &f)
  {
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
  friend std::istream &operator>>(std::istream &stream, Multi_parameter_filtration &f)
  {
    size_type num_gen;
    size_type num_param;
    char delimiter;
    stream >> delimiter;  // (
    stream >> delimiter;  // k
    stream >> delimiter;  // =
    stream >> num_gen;
    if (!stream.good()) throw std::invalid_argument("Invalid incoming stream format for Multi_parameter_filtration.");
    stream >> delimiter;  // )
    stream >> delimiter;  // (
    stream >> delimiter;  // p
    stream >> delimiter;  // =
    stream >> num_param;
    if (!stream.good()) throw std::invalid_argument("Invalid incoming stream format for Multi_parameter_filtration.");
    f.generators_.resize(num_gen * num_param);
    f.generator_view_ = Viewer(f.generators_.data(), num_gen, num_param);
    stream >> delimiter;  // )
    stream >> delimiter;  // [
    if (delimiter != '[') throw std::invalid_argument("Invalid incoming stream format for Multi_parameter_filtration.");
    if (num_gen == 0) return stream;
    for (size_type i = 0; i < num_gen; ++i) {
      stream >> delimiter;  // [
      for (size_type j = 0; j < num_param; ++j) {
        f(i, j) = _get_value<T>(stream);
        if (!stream.good())
          throw std::invalid_argument("Invalid incoming stream format for Multi_parameter_filtration.");
        stream >> delimiter;  // , or last ]
      }
      stream >> delimiter;  // ; or last ]
    }
    if (delimiter != ']') throw std::invalid_argument("Invalid incoming stream format for Multi_parameter_filtration.");

    return stream;
  }

  /**
   * @brief Returns true if and only if the given filtration value is at plus infinity.
   */
  friend bool is_positive_infinity(const Multi_parameter_filtration &f)
  {
    return f.is_plus_inf();
  }

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
  friend bool unify_lifetimes(Multi_parameter_filtration &f1, const Multi_parameter_filtration &f2)
  {
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
        auto start = f2.begin();
        start += g * f2.num_parameters();
        auto end = start + f2.num_parameters();
        modified |= f1.add_generator(start, end);
      }
      return modified;
    }
  }

  /**
   * @brief Stores in the first argument the origins of the cones in the intersection of the positive
   * (negative if `Co` is true) cones generated by the two arguments.
   *
   * @param f1 First set of cones which will be modified.
   * @param f2 Second set of cones. Should have the same number of parameters than the first one.
   * @return true If the first argument was actually modified.
   * @return false Otherwise.
   */
  friend bool intersect_lifetimes(Multi_parameter_filtration &f1, const Multi_parameter_filtration &f2)
  {
    GUDHI_CHECK(f1.num_parameters() == f2.num_parameters(),
                std::invalid_argument("Cannot intersect two filtration values with different number of parameters."));

    if constexpr (Ensure1Criticality) {
      if constexpr (Co) {
        return f1.pull_to_greatest_common_lower_bound(f2);
      } else {
        return f1.push_to_least_common_upper_bound(f2);
      }
    } else {
      bool f2IsInf = true, f2IsMinusInf = true, f2IsNaN = true;
      bool f1IsInf = true, f1IsMinusInf = true, f1IsNaN = true;

      // if one is not finite, we can avoid the heavy simplification process
      _get_infinity_statuses(f1.generator_view_, f2, f1IsInf, f1IsMinusInf, f1IsNaN, f2IsInf, f2IsMinusInf, f2IsNaN);

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

      const size_type num_param = f1.num_parameters();
      Multi_parameter_filtration res(num_param, Co ? T_m_inf : T_inf);
      std::vector<T> newGen(num_param);
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
  friend char *serialize_value_to_char_buffer(const Multi_parameter_filtration &value, char *start)
  {
    const size_type length = value.generators_.size();
    const size_type num_param = value.num_parameters();
    const std::size_t arg_size = sizeof(T) * length;
    const std::size_t type_size = sizeof(size_type);
    memcpy(start, &num_param, type_size);
    memcpy(start + type_size, &length, type_size);
    memcpy(start + (type_size * 2), value.generators_.data(), arg_size);
    return start + arg_size + (type_size * 2);
  }

  /**
   * @brief Deserialize the value from a buffer at given pointer and stores it in given value.
   *
   * @param value Value to fill with the deserialized filtration value.
   * @param start Pointer to the start of the space in the buffer where the serialization is stored.
   * @return End position of the serialization in the buffer.
   */
  friend const char *deserialize_value_from_char_buffer(Multi_parameter_filtration &value, const char *start)
  {
    const std::size_t type_size = sizeof(size_type);
    size_type length;
    size_type num_param;
    memcpy(&num_param, start, type_size);
    memcpy(&length, start + type_size, type_size);
    std::size_t arg_size = sizeof(T) * length;
    value.generators_.resize(length);
    memcpy(value.generators_.data(), start + (type_size * 2), arg_size);
    value.generator_view_ =
        Viewer(value.generators_.data(), num_param == 0 ? 0 : value.generators_.size() / num_param, num_param);
    return start + arg_size + (type_size * 2);
  }

  /**
   * @brief Returns the serialization size of the given filtration value.
   */
  friend std::size_t get_serialization_size_of(const Multi_parameter_filtration &value)
  {
    return (sizeof(size_type) * 2) + (sizeof(T) * value.num_entries());
  }

  /**
   * @brief Plus infinity value of an entry of the filtration value.
   */
  constexpr static const T T_inf = MF_T_inf<T>;

  /**
   * @brief Minus infinity value of an entry of the filtration value.
   */
  constexpr static const T T_m_inf = MF_T_m_inf<T>;

 private:
  Underlying_container generators_; /**< Container of the filtration value elements. */
  Viewer generator_view_;           /**< Matrix view of the container. Has to be created after generators_. */

  /**
   * @brief Default value of an element in the filtration value.
   */
  constexpr static T _get_default_value() { return Co ? T_inf : T_m_inf; }

  /**
   * @brief Verifies if @p b is strictly contained in the positive cone originating in `a`.
   */
  static bool _strictly_contains(const Viewer &a, size_type g_a, const Viewer &b, size_type g_b)
  {
    bool isSame = true;
    for (auto i = 0U; i < a.extent(1); ++i) {
      T a_i, b_i;
      if constexpr (Co) {
        a_i = b(g_b, i);
        b_i = a(g_a, i);
      } else {
        a_i = a(g_a, i);
        b_i = b(g_b, i);
      }
      if (a_i > b_i || _is_nan(a_i) || _is_nan(b_i)) return false;
      if (isSame && a_i != b_i) isSame = false;
    }
    return !isSame;
  }

  /**
   * @brief Verifies if @p b is contained in the positive cone originating in `a`.
   */
  static bool _contains(const Viewer &a, size_type g_a, const Viewer &b, size_type g_b)
  {
    for (std::size_t i = 0U; i < a.extent(1); ++i) {
      T a_i, b_i;
      if constexpr (Co) {
        a_i = b(g_b, i);
        b_i = a(g_a, i);
      } else {
        a_i = a(g_a, i);
        b_i = b(g_b, i);
      }
      if (a_i > b_i || (!_is_nan(a_i) && _is_nan(b_i)) || (_is_nan(a_i) && !_is_nan(b_i))) return false;
    }
    return true;
  }

  /**
   * @brief Verifies if the first element of @p b strictly dominates the first element of `a`.
   */
  static bool _first_strictly_dominates(const Viewer &a, size_type g_a, const Viewer &b, size_type g_b)
  {
    if constexpr (Co) {
      return a(g_a, 0) < b(g_b, 0);
    } else {
      return a(g_a, 0) > b(g_b, 0);
    }
  }

  /**
   * @brief Verifies if the first element of @p b dominates the first element of `a`.
   */
  static bool _first_dominates(const Viewer &a, size_type g_a, const Viewer &b, size_type g_b)
  {
    if constexpr (Co) {
      return a(g_a, 0) <= b(g_b, 0);
    } else {
      return a(g_a, 0) >= b(g_b, 0);
    }
  }

  /**
   * @brief Applies operation on the elements of the filtration value.
   */
  template <class ValueRange, class F, class = std::enable_if_t<RangeTraits<ValueRange>::has_begin> >
  void _apply_operation(const ValueRange &range, F &&operate)
  {
    auto &view = generator_view_;
    for (unsigned int g = 0; g < num_generators(); ++g) {
      auto it = range.begin();
      for (unsigned int p = 0; p < num_parameters() && it != range.end(); ++p) {
        std::forward<F>(operate)(view(g, p), *it);
        ++it;
      }
    }
  }

  /**
   * @brief Applies operation on the elements of the filtration value.
   */
  template <class F>
  void _apply_operation(const T &val, F &&operate)
  {
    auto &gens = generators_;
    for (unsigned int i = 0; i < gens.size(); ++i) {
      std::forward<F>(operate)(gens[i], val);
    }
  }

  template <class GeneratorRange>
  static void _get_infinity_statuses(const Viewer &a,
                                     const GeneratorRange &b,
                                     bool &aIsInf,
                                     bool &aIsMinusInf,
                                     bool &aIsNaN,
                                     bool &bIsInf,
                                     bool &bIsMinusInf,
                                     bool &bIsNaN)
  {
    auto itB = b.begin();
    for (std::size_t i = 0; i < a.extent(1); ++i) {
      if (a(0, i) != T_inf) aIsInf = false;
      if (a(0, i) != T_m_inf) aIsMinusInf = false;
      if (!_is_nan(a(0, i))) aIsNaN = false;
      if (*itB != T_inf) bIsInf = false;
      if (*itB != T_m_inf) bIsMinusInf = false;
      if (!_is_nan(*itB)) bIsNaN = false;
      if (!aIsInf && !aIsMinusInf && !aIsNaN && !bIsInf && !bIsMinusInf && !bIsNaN) return;
      ++itB;
    }
  }

  enum class Rel : std::uint8_t { EQUAL, DOMINATES, IS_DOMINATED, NONE };

  template <class Iterator>
  static Rel _get_domination_relation(const Viewer &a, size_type g_a, Iterator itB)
  {
    bool equal = true;
    bool allGreater = true;
    bool allSmaller = true;
    bool allNaNA = true;
    bool allNaNB = true;
    for (unsigned int i = 0; i < a.extent(1); ++i) {
      if (a(g_a, i) < *itB) {
        if (!allSmaller) return Rel::NONE;
        equal = false;
        allGreater = false;
      } else if (a(g_a, i) > *itB) {
        if (!allGreater) return Rel::NONE;
        equal = false;
        allSmaller = false;
      }
      if (!_is_nan(a(g_a, i))) allNaNA = false;
      if (!_is_nan(*itB)) allNaNB = false;
      ++itB;
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
  bool _generator_can_be_added(Iterator x, std::size_t curr, std::size_t &end, std::vector<int> &indices)
  {
    while (curr != end) {
      Rel res = _get_domination_relation(generator_view_, indices[curr], x);
      if (res == Rel::IS_DOMINATED || res == Rel::EQUAL) return false;  // x dominates or is equal
      if (res == Rel::DOMINATES) {                                      // x is dominated
        --end;
        std::swap(indices[curr], indices[end]);
      } else {  // no relation
        ++curr;
      }
    }
    return true;
  }

  /**
   * @brief Rebuild the generators from the given set.
   */
  template <class Iterator>
  void _build_from(std::vector<int> &indices, const int newIndex, Iterator xStart, Iterator xEnd)
  {
    auto comp = [&](int g1, int g2) -> bool {
      if (g1 == g2) {
        return false;
      }

      if (g1 == newIndex) {
        auto it = xStart;
        for (size_type i = 0; i < num_parameters(); ++i) {
          T v = generator_view_(g2, i);
          if (*it != v) {
            if (_is_nan(v)) return true;
            return *it < v;
          }
          ++it;
        }
        return false;
      }

      if (g2 == newIndex) {
        auto it = xStart;
        for (size_type i = 0; i < num_parameters(); ++i) {
          T v = generator_view_(g1, i);
          if (v != *it) {
            if (_is_nan(*it)) return true;
            return v < *it;
          }
          ++it;
        }
        return false;
      }

      for (size_type i = 0; i < num_parameters(); ++i) {
        T v1 = generator_view_(g1, i);
        T v2 = generator_view_(g2, i);
        if (v1 != v2) {
          if (_is_nan(v2)) return true;
          return v1 < v2;
        }
      }
      return false;
    };
    std::sort(indices.begin(), indices.end(), comp);

    Underlying_container new_container;
    new_container.reserve((indices.size() + 1) * num_parameters());
    for (int i : indices) {
      if (i == newIndex) {
        new_container.insert(new_container.end(), xStart, xEnd);
      } else {
        T *ptr = &generator_view_(i, 0);
        for (size_type p = 0; p < num_parameters(); ++p) {
          new_container.push_back(*ptr);
          ++ptr;
        }
      }
    }
    generators_.swap(new_container);
    generator_view_ = Viewer(generators_.data(), generators_.size() / num_parameters(), num_parameters());
  }

  /**
   * @brief Rebuild the generators from the given set.
   */
  void _build_from(std::vector<int> &indices)
  {
    auto comp = [&](int g1, int g2) -> bool {
      for (std::size_t i = 0; i < num_parameters(); ++i) {
        T v1 = generator_view_(g1, i);
        T v2 = generator_view_(g2, i);
        if (v1 != v2) {
          if (_is_nan(v2)) return true;
          return v1 < v2;
        }
      }
      return false;
    };
    std::sort(indices.begin(), indices.end(), comp);

    Underlying_container new_container;
    new_container.reserve(indices.size() * num_parameters());
    for (int i : indices) {
      T *ptr = &generator_view_(i, 0);
      for (size_type p = 0; p < num_parameters(); ++p) {
        new_container.push_back(*ptr);
        ++ptr;
      }
    }
    generators_.swap(new_container);
    generator_view_ = Viewer(generators_.data(), generators_.size() / num_parameters(), num_parameters());
  }

  bool _is_finite(size_type g)
  {
    bool isInf = true, isMinusInf = true, isNan = true;
    for (size_type p = 0; p < num_parameters(); ++p) {
      T v = generator_view_(g, p);
      if (v != T_inf) isInf = false;
      if (v != T_m_inf) isMinusInf = false;
      if (!_is_nan(v)) isNan = false;
      if (!isInf && !isMinusInf && !isNan) return true;
    }
    return false;
  }

  template <class F, typename U = T>
  static U _compute_frobenius_norm(size_type number_of_elements, F &&norm)
  {
    if (number_of_elements == 1) return std::forward<F>(norm)(0);

    U out = 0;
    for (size_type p = 0; p < number_of_elements; ++p) {
      T v = std::forward<F>(norm)(p);
      out += v * v;
    }
    if constexpr (std::is_integral_v<U>) {
      // to avoid Windows issue that don't know how to cast integers for cmath methods
      return std::sqrt(static_cast<double>(out));
    } else {
      return std::sqrt(out);
    }
  }
};

}  // namespace Gudhi::multi_filtration

namespace std {

template <typename T, bool Co, bool Ensure1Criticality>
class numeric_limits<Gudhi::multi_filtration::Multi_parameter_filtration<T, Co, Ensure1Criticality> >
{
 public:
  using Filtration_value = Gudhi::multi_filtration::Multi_parameter_filtration<T, Co, Ensure1Criticality>;

  static constexpr bool has_infinity = true;
  static constexpr bool has_quiet_NaN = std::numeric_limits<T>::has_quiet_NaN;

  static constexpr Filtration_value infinity(std::size_t p = 1) noexcept { return Filtration_value::inf(p); };

  // non-standard
  static constexpr Filtration_value minus_infinity(std::size_t p = 1) noexcept
  {
    return Filtration_value::minus_inf(p);
  };

  static constexpr Filtration_value max() noexcept(false)
  {
    throw std::logic_error(
        "The max value cannot be represented with no finite numbers of parameters."
        "Use `max(number_of_parameters)` instead");
  };

  static constexpr Filtration_value max(std::size_t p) noexcept
  {
    return Filtration_value(p, std::numeric_limits<T>::max());
  };

  static constexpr Filtration_value lowest(std::size_t p = 1) noexcept { return Filtration_value::minus_inf(p); };

  static constexpr Filtration_value quiet_NaN(std::size_t p = 1) noexcept(false)
  {
    if constexpr (std::numeric_limits<T>::has_quiet_NaN) {
      return Filtration_value::nan(p);
    } else {
      throw std::logic_error("Does not have a NaN value.");
    }
  };
};

}  // namespace std

#endif  // MF_MULTI_PARAMETER_FILTRATION_H_
