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
 * @private
 * @file Multi_parameter_generator.h
 * @author Hannah Schreiber, David Loiseaux
 * @brief Contains the @ref Gudhi::multi_filtration::Multi_parameter_generator class.
 */

#ifndef MF_MULTI_PARAMETER_GENERATOR_H_
#define MF_MULTI_PARAMETER_GENERATOR_H_

#include <algorithm>    //std::lower_bound
#include <cmath>        //std::isnan, std::min
#include <cstddef>      //std::size_t
#include <cstring>      //memcpy
#include <iterator>     //std::distance
#include <ostream>      //std::ostream
#include <limits>       //std::numerical_limits
#include <stdexcept>    //std::logic_error
#include <type_traits>  //std::is_arithmetic
#include <utility>      //std::swap, std::move
#include <vector>
#include <initializer_list>

#include <gudhi/Debug_utils.h>
#include <gudhi/simple_mdspan.h>
#include <gudhi/Multi_filtration/multi_filtration_utils.h>

namespace Gudhi::multi_filtration {

// declaration needed pre C++20 for friends with templates defined inside a class
template <typename U>
U compute_euclidean_distance_to();
template <typename U>
U compute_norm();

/**
 * private
 * @class Multi_parameter_generator Multi_parameter_generator.h gudhi/Multi_parameter_generator.h
 * @ingroup multi_filtration
 *
 * @brief Class encoding a single generator of @ref Dynamic_multi_parameter_filtration "".
 *
 * @tparam T Arithmetic type of an entry for one parameter of a filtration value. Has to be **signed** and
 * to implement `std::isnan(T)`, `std::numeric_limits<T>::has_quiet_NaN`, `std::numeric_limits<T>::quiet_NaN()`,
 * `std::numeric_limits<T>::has_infinity`, `std::numeric_limits<T>::infinity()` and `std::numeric_limits<T>::max()`.
 * If `std::numeric_limits<T>::has_infinity` returns `false`, a call to `std::numeric_limits<T>::infinity()`
 * can simply throw. Examples are the native types `double`, `float` and `int`.
 */
template <typename T>
class Multi_parameter_generator
{
 public:
  using Underlying_container = std::vector<T>; /**< Underlying container for values. */

  // CONSTRUCTORS

  /**
   * @brief Default constructor. Builds generator with given number of parameters.
   * All values are set at -inf.
   *
   * @param number_of_parameters If negative, takes the default value instead. Default value: 1.
   */
  Multi_parameter_generator(int number_of_parameters = 1)
      : generator_(number_of_parameters < 0 ? 1 : number_of_parameters, _get_default_value())
  {}

  /**
   * @brief Builds generator with given number of parameters. All values are initialized at the given value.
   *
   * @param number_of_parameters If negative, is set to 1 instead.
   * @param value Initialization value for every value in the generator.
   */
  Multi_parameter_generator(int number_of_parameters, T value)
      : generator_(number_of_parameters < 0 ? 1 : number_of_parameters, value)
  {}

  /**
   * @brief Builds generator that is initialized with the given range. The number of
   * parameters are therefore deduced from the length of the range.
   *
   * @tparam ValueRange Range of types convertible to `T`. Should have a begin() and end() method.
   * @param range Values of the generator.
   */
  template <class ValueRange = std::initializer_list<T>, class = std::enable_if_t<RangeTraits<ValueRange>::has_begin> >
  Multi_parameter_generator(const ValueRange &range) : generator_(range.begin(), range.end())
  {}

  /**
   * @brief Builds generator that is initialized with the given range. The range is
   * determined from the two given iterators. The number of parameters are therefore deduced from the distance
   * between the two.
   *
   * @tparam Iterator Iterator type that has to satisfy the requirements of standard LegacyInputIterator and
   * dereferenced elements have to be convertible to `T`.
   * @param it_begin Iterator pointing to the start of the range.
   * @param it_end Iterator pointing to the end of the range.
   */
  template <class Iterator>
  Multi_parameter_generator(Iterator it_begin, Iterator it_end) : generator_(it_begin, it_end)
  {}

  /**
   * @brief Builds generator with given number of parameters and values from the given range.
   * The range is represented by @ref Multi_parameter_generator::Underlying_container "" and **moved** into the
   * underlying container of the class.
   *
   * @param generators Values to move.
   */
  Multi_parameter_generator(Underlying_container &&generators) : generator_(std::move(generators)) {}

  /**
   * @brief Copy constructor.
   */
  Multi_parameter_generator(const Multi_parameter_generator &other) = default;

  /**
   * @brief Copy constructor.
   *
   * @tparam U Type convertible into `T`.
   */
  template <typename U>
  Multi_parameter_generator(const Multi_parameter_generator<U> &other) : generator_(other.begin(), other.end())
  {}

  /**
   * @brief Move constructor.
   */
  Multi_parameter_generator(Multi_parameter_generator &&other) noexcept = default;

  ~Multi_parameter_generator() = default;

  /**
   * @brief Assign operator.
   */
  Multi_parameter_generator &operator=(const Multi_parameter_generator &other) = default;

  /**
   * @brief Move assign operator.
   */
  Multi_parameter_generator &operator=(Multi_parameter_generator &&other) noexcept = default;

  /**
   * @brief Assign operator.
   *
   * @tparam U Type convertible into `T`.
   */
  template <typename U>
  Multi_parameter_generator &operator=(const Multi_parameter_generator<U> &other)
  {
    generator_ = Underlying_container(other.begin(), other.end());
    return *this;
  }

  /**
   * @brief Swap operator.
   */
  friend void swap(Multi_parameter_generator &g1, Multi_parameter_generator &g2) noexcept
  {
    g1.generator_.swap(g2.generator_);
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
   * @brief Returns reference to the value of parameter `p`.
   */
  reference operator[](size_type p)
  {
    GUDHI_CHECK(p < generator_.size(), std::out_of_range("Out of bound parameter."));
    return generator_[p];
  }

  /**
   * @brief Returns const reference to the value of parameter `p`.
   */
  const_reference operator[](size_type p) const
  {
    GUDHI_CHECK(p < generator_.size(), std::out_of_range("Out of bound parameter."));
    return generator_[p];
  }

  /**
   * @brief Returns an iterator pointing the begining of the underlying container.
   */
  iterator begin() noexcept { return generator_.begin(); }

  /**
   * @brief Returns an iterator pointing the begining of the underlying container.
   */
  const_iterator begin() const noexcept { return generator_.begin(); }

  /**
   * @brief Returns an iterator pointing the begining of the underlying container.
   */
  const_iterator cbegin() const noexcept { return generator_.cbegin(); }

  /**
   * @brief Returns an iterator pointing the end of the underlying container.
   */
  iterator end() noexcept { return generator_.end(); }

  /**
   * @brief Returns an iterator pointing the end of the underlying container.
   */
  const_iterator end() const noexcept { return generator_.end(); }

  /**
   * @brief Returns an iterator pointing the end of the underlying container.
   */
  const_iterator cend() const noexcept { return generator_.cend(); }

  /**
   * @brief Returns a reverse iterator pointing to the first element from the back of the underlying container.
   */
  reverse_iterator rbegin() noexcept { return generator_.rbegin(); }

  /**
   * @brief Returns a reverse iterator pointing to the first element from the back of the underlying container.
   */
  const_reverse_iterator rbegin() const noexcept { return generator_.rbegin(); }

  /**
   * @brief Returns a reverse iterator pointing to the first element from the back of the underlying container.
   */
  const_reverse_iterator crbegin() const noexcept { return generator_.crbegin(); }

  /**
   * @brief Returns a reverse iterator pointing to the end of the reversed underlying container.
   */
  reverse_iterator rend() noexcept { return generator_.rend(); }

  /**
   * @brief Returns a reverse iterator pointing to the end of the reversed underlying container.
   */
  const_reverse_iterator rend() const noexcept { return generator_.rend(); }

  /**
   * @brief Returns a reverse iterator pointing to the end of the reversed underlying container.
   */
  const_reverse_iterator crend() const noexcept { return generator_.crend(); }

  /**
   * @brief Returns the size of the underlying container. Corresponds exactly to @ref num_parameters(), but enables to
   * use the class as a classic range with a `begin`, `end` and `size` method.
   */
  size_type size() const noexcept { return generator_.size(); }

  // CONVERTERS

  // like numpy
  /**
   * @brief Returns a copy with entries casted into the type given as template parameter.
   *
   * @tparam U New type for the entries.
   * @return Copy with new entry type.
   */
  template <typename U>
  Multi_parameter_generator<U> as_type() const
  {
    std::vector<U> out(generator_.begin(), generator_.end());
    return Multi_parameter_generator<U>(std::move(out));
  }

  // ACCESS

  /**
   * @brief Returns the number of parameters in the filtration value.
   */
  size_type num_parameters() const { return generator_.size(); }

  /**
   * @brief Returns a generator for which @ref is_plus_inf() returns `true`.
   */
  static Multi_parameter_generator inf() { return Multi_parameter_generator(1, T_inf); }

  /**
   * @brief Returns a generator for which @ref is_minus_inf() returns `true`.
   */
  static Multi_parameter_generator minus_inf() { return Multi_parameter_generator(1, T_m_inf); }

  /**
   * @brief Returns a generator for which @ref is_nan() returns `true`.
   */
  static Multi_parameter_generator nan()
  {
    if constexpr (std::numeric_limits<T>::has_quiet_NaN) {
      return Multi_parameter_generator(1, std::numeric_limits<T>::quiet_NaN());
    } else {
      return Multi_parameter_generator(0);
    }
  }

  // DESCRIPTORS

  /**
   * @brief Returns `true` if and only if the generator is considered as plus infinity.
   */
  [[nodiscard]] bool is_plus_inf() const
  {
    for (const T &v : generator_) {
      if (v != T_inf) return false;
    }
    return !generator_.empty();
  }

  /**
   * @brief Returns `true` if and only if the generator is considered as minus infinity.
   */
  [[nodiscard]] bool is_minus_inf() const
  {
    for (const T &v : generator_) {
      if (v != T_m_inf) return false;
    }
    return !generator_.empty();
  }

  /**
   * @brief Returns `true` if and only if the generator is considered as NaN.
   */
  [[nodiscard]] bool is_nan() const
  {
    if constexpr (std::numeric_limits<T>::has_quiet_NaN) {
      for (const T &v : generator_) {
        if (!std::isnan(v)) return false;
      }
      return true;
    } else {
      return generator_.empty();
    }
  }

  /**
   * @brief Returns `true` if and only if the generator is non-empty and is not considered as plus infinity,
   * minus infinity or NaN.
   */
  [[nodiscard]] bool is_finite() const
  {
    bool isInf = true, isMinusInf = true, isNan = true;
    for (const auto &v : generator_) {
      if (v != T_inf) isInf = false;
      if (v != T_m_inf) isMinusInf = false;
      if (!_is_nan(v)) isNan = false;
      if (!isInf && !isMinusInf && !isNan) return true;
    }
    return false;
  }

  // COMPARAISON OPERATORS

  /**
   * @brief Returns `true` if and only if the positive cone generated by @p b is strictly contained in the
   * cone generated by @p a. Both @p a and @p b  have to have the same number of parameters.
   *
   * Note that not all filtration values are comparable. That is, \f$ a < b \f$ and \f$ b < a \f$ returning both false
   * does **not** imply \f$ a == b \f$.
   */
  friend bool operator<(const Multi_parameter_generator &a, const Multi_parameter_generator &b)
  {
    if (&a == &b) return false;
    if (a.is_nan() || b.is_nan()) return false;
    if (a.is_minus_inf() || b.is_plus_inf()) return true;
    if (b.is_minus_inf() || a.is_plus_inf()) return false;

    GUDHI_CHECK(a.num_parameters() == b.num_parameters(),
                std::invalid_argument("Only generators with same number of parameters can be compared."));

    bool isSame = true;
    auto n = a.size();
    for (size_type i = 0U; i < n; ++i) {
      if (a[i] > b[i]) return false;
      if (isSame && a[i] != b[i]) isSame = false;
    }
    return !isSame;
  }

  /**
   * @brief Returns `true` if and only if the positive cone generated by @p a is strictly contained in the
   * cone generated by @p b. Both @p a and @p b  have to have the same number of parameters.
   *
   * Note that not all filtration values are comparable. That is, \f$ a \le b \f$ and \f$ b \le a \f$ can both return
   * `false`.
   */
  friend bool operator<=(const Multi_parameter_generator &a, const Multi_parameter_generator &b)
  {
    if (a.is_nan() || b.is_nan()) return false;
    if (&a == &b) return true;
    const bool aIsMinusInf = a.is_minus_inf();
    const bool aIsPlusInf = a.is_plus_inf();
    const bool bIsMinusInf = b.is_minus_inf();
    const bool bIsPlusInf = b.is_plus_inf();
    if ((aIsMinusInf && bIsMinusInf) || (aIsPlusInf && bIsPlusInf)) return true;
    if (aIsMinusInf || bIsPlusInf) return true;
    if (bIsMinusInf || aIsPlusInf) return false;

    GUDHI_CHECK(a.num_parameters() == b.num_parameters(),
                std::invalid_argument("Only generators with same number of parameters can be compared."));

    auto n = a.size();
    for (std::size_t i = 0U; i < n; ++i) {
      if (a[i] > b[i]) return false;
    }
    return true;
  }

  /**
   * @brief Returns `true` if and only if the positive cone generated by @p b is contained in or are (partially)
   * equal to the positive cone generated by @p a.
   * Both @p a and @p b  have to have the same number of parameters.
   *
   * Note that not all filtration values are comparable. That is, \f$ a > b \f$ and \f$ b > a \f$ returning both false
   * does **not** imply \f$ a == b \f$.
   */
  friend bool operator>(const Multi_parameter_generator &a, const Multi_parameter_generator &b) { return b < a; }

  /**
   * @brief Returns `true` if and only if the positive cone generated by @p a is contained in or are (partially)
   * equal to the positive cone generated by @p b.
   * Both @p a and @p b  have to have the same number of parameters.
   *
   * Note that not all filtration values are comparable. That is, \f$ a \ge b \f$ and \f$ b \ge a \f$ can both return
   * `false`.
   */
  friend bool operator>=(const Multi_parameter_generator &a, const Multi_parameter_generator &b) { return b <= a; }

  /**
   * @brief Returns `true` if and only if for each \f$ i,j \f$, \f$ a(i,j) \f$ is equal to \f$ b(i,j) \f$.
   */
  friend bool operator==(const Multi_parameter_generator &a, const Multi_parameter_generator &b)
  {
    if (a.is_nan() || b.is_nan()) return false;
    if (&a == &b) return true;
    if (a.is_plus_inf()) {
      return b.is_plus_inf();
    }
    if (a.is_minus_inf()) {
      return b.is_minus_inf();
    }
    return a.generator_ == b.generator_;
  }

  /**
   * @brief Returns `true` if and only if \f$ a == b \f$ returns `false`.
   */
  friend bool operator!=(const Multi_parameter_generator &a, const Multi_parameter_generator &b) { return !(a == b); }

  // ARITHMETIC OPERATORS

  // opposite
  /**
   * @brief Returns a generator such that an entry at index \f$ i \f$ is equal to \f$ -g(i) \f$.
   *
   * Used conventions:
   * - \f$ -NaN = NaN \f$.
   *
   * @param g Value to opposite.
   * @return The opposite of @p f.
   */
  friend Multi_parameter_generator operator-(const Multi_parameter_generator &g)
  {
    using F = Multi_parameter_generator;

    std::vector<T> result(g.generator_);
    std::for_each(result.begin(), result.end(), [](T &v) {
      if (v == F::T_inf)
        v = F::T_m_inf;
      else if (v == F::T_m_inf)
        v = F::T_inf;
      else
        v = -v;
    });
    return Multi_parameter_generator(std::move(result));
  }

  // subtraction
  /**
   * @brief Returns a generator such that an entry at index \f$ p \f$, with
   * \f$ 0 \leq p \leq num_parameters \f$ is equal to \f$ g(p) - r(p) \f$
   * if \f$ p < length_r \f$ and to \f$ g(p) \f$ otherwise.
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
   * @tparam ValueRange Range with a begin(), end() and size() method.
   * @param g First element of the subtraction.
   * @param r Second element of the subtraction.
   */
  template <class ValueRange, class = std::enable_if_t<RangeTraits<ValueRange>::has_begin> >
  friend Multi_parameter_generator operator-(Multi_parameter_generator g, const ValueRange &r)
  {
    g -= r;
    return g;
  }

  /**
   * @brief Returns a generator such that an entry at index \f$ p \f$, with
   * \f$ 0 \leq p \leq num_parameters \f$ is equal to \f$ r(p) - g(p) \f$
   * if \f$ p < length_r \f$ and to \f$ g(p) \f$ otherwise.
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
   * @tparam ValueRange Range with a begin(), end() and size() method.
   * @param r First element of the subtraction.
   * @param g Second element of the subtraction.
   */
  template <class ValueRange,
            class = std::enable_if_t<RangeTraits<ValueRange>::has_begin &&
                                     !std::is_same_v<ValueRange, Multi_parameter_generator> > >
  friend Multi_parameter_generator operator-(const ValueRange &r, const Multi_parameter_generator &g)
  {
    if (g.num_parameters() == 0) return g;
    if (g.is_nan()) return nan();

    if constexpr (RangeTraits<ValueRange>::is_dynamic_multi_filtration) {
      if (r[0].is_finite()) {
        return g -= r[0];
      }
      if constexpr (!std::numeric_limits<T>::has_quiet_NaN)
        if (r[0].size() == 0) return nan();
      return g -= r[0][0];
    } else {
      if (g.is_finite()) {
        Multi_parameter_generator res = g;
        res._apply_operation(r, [](T &valF, const T &valR) -> bool {
          valF = -valF;
          return _add(valF, valR);
        });
        return res;
      }
      Multi_parameter_generator res(r.begin(), r.end());
      res._apply_operation(g[0], [](T &valRes, const T &valF) -> bool { return _subtract(valRes, valF); });
      return res;
    }
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ p \f$ is equal to \f$ g(p) - val \f$.
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
   * @param g First element of the subtraction.
   * @param val Second element of the subtraction.
   */
  friend Multi_parameter_generator operator-(Multi_parameter_generator g, const T &val)
  {
    g -= val;
    return g;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ p \f$ is equal to \f$ val - g(p) \f$.
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
   * @param g Second element of the subtraction.
   */
  friend Multi_parameter_generator operator-(const T &val, Multi_parameter_generator g)
  {
    g._apply_operation(val, [](T &valF, const T &valR) -> bool {
      valF = -valF;
      return _add(valF, valR);
    });
    return g;
  }

  /**
   * @brief Modifies the first argument such that an entry at index \f$ p \f$, with
   * \f$ 0 \leq p \leq num_parameters \f$ is equal to \f$ g(p) - r(p) \f$
   * if \f$ p < length_r \f$ and to \f$ g(p) \f$ otherwise.
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
   * @tparam ValueRange Range with a begin(), end() and size() method.
   * @param g First element of the subtraction.
   * @param r Second element of the subtraction.
   */
  template <class ValueRange, class = std::enable_if_t<RangeTraits<ValueRange>::has_begin> >
  friend Multi_parameter_generator &operator-=(Multi_parameter_generator &g, const ValueRange &r)
  {
    if (g.num_parameters() == 0 || g.is_nan()) return g;

    if constexpr (RangeTraits<ValueRange>::is_dynamic_multi_filtration) {
      if (r[0].is_finite()) {
        return g -= r[0];
      }
      if constexpr (!std::numeric_limits<T>::has_quiet_NaN)
        if (r[0].size() == 0) {
          g = nan();
          return g;
        }
      return g -= r[0][0];
    } else {
      if (!g.is_finite()) {
        g.generator_.resize(r.size(), g[0]);
      }
      g._apply_operation(r, [](T &valF, const T &valR) -> bool { return _subtract(valF, valR); });
    }

    return g;
  }

  /**
   * @brief Modifies the first parameter such that an entry at index \f$ p \f$ is equal to \f$ g(p) - val \f$.
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
   * @param g First element of the subtraction.
   * @param val Second element of the subtraction.
   */
  friend Multi_parameter_generator &operator-=(Multi_parameter_generator &g, const T &val)
  {
    g._apply_operation(val, [](T &valF, const T &valR) -> bool { return _subtract(valF, valR); });
    return g;
  }

  // addition
  /**
   * @brief Returns a generator such that an entry at index \f$ p \f$, with
   * \f$ 0 \leq p \leq num_parameters \f$ is equal to \f$ g(p) + r(p) \f$
   * if \f$ p < length_r \f$ and to \f$ g(p) \f$ otherwise.
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
   * @tparam ValueRange Range with a begin(), end() and size() method.
   * @param g First element of the addition.
   * @param r Second element of the addition.
   */
  template <class ValueRange, class = std::enable_if_t<RangeTraits<ValueRange>::has_begin> >
  friend Multi_parameter_generator operator+(Multi_parameter_generator g, const ValueRange &r)
  {
    g += r;
    return g;
  }

  /**
   * @brief Returns a generator such that an entry at index \f$ p \f$, with
   * \f$ 0 \leq p \leq num_parameters \f$ is equal to \f$ r(p) + g(p) \f$
   * if \f$ p < length_r \f$ and to \f$ g(p) \f$ otherwise.
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
   * @tparam ValueRange Range with a begin(), end() and size() method.
   * @param r First element of the addition.
   * @param g Second element of the addition.
   */
  template <class ValueRange,
            class = std::enable_if_t<RangeTraits<ValueRange>::has_begin &&
                                     !std::is_same_v<ValueRange, Multi_parameter_generator> > >
  friend Multi_parameter_generator operator+(const ValueRange &r, Multi_parameter_generator g)
  {
    g += r;
    return g;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ p \f$ is equal to \f$ g(p) + val \f$.
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
   * @param g First element of the addition.
   * @param val Second element of the addition.
   */
  friend Multi_parameter_generator operator+(Multi_parameter_generator g, const T &val)
  {
    g += val;
    return g;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ p \f$ is equal to \f$ g(p) + val \f$.
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
   * @param g Second element of the addition.
   */
  friend Multi_parameter_generator operator+(const T &val, Multi_parameter_generator g)
  {
    g += val;
    return g;
  }

  /**
   * @brief Modifies the first argument such that an entry at index \f$ p \f$, with
   * \f$ 0 \leq p \leq num_parameters \f$ is equal to \f$ g(p) + r(p) \f$
   * if \f$ p < length_r \f$ and to \f$ g(p) \f$ otherwise.
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
   * @tparam ValueRange Range with a begin(), end() and size() method.
   * @param g First element of the addition.
   * @param r Second element of the addition.
   */
  template <class ValueRange, class = std::enable_if_t<RangeTraits<ValueRange>::has_begin> >
  friend Multi_parameter_generator &operator+=(Multi_parameter_generator &g, const ValueRange &r)
  {
    if (g.num_parameters() == 0 || g.is_nan()) return g;

    if constexpr (RangeTraits<ValueRange>::is_dynamic_multi_filtration) {
      if (r[0].is_finite()) {
        return g += r[0];
      }
      if constexpr (!std::numeric_limits<T>::has_quiet_NaN)
        if (r[0].size() == 0) {
          g = nan();
          return g;
        }
      return g += r[0][0];
    } else {
      if (!g.is_finite()) {
        g.generator_.resize(r.size(), g[0]);
      }
      g._apply_operation(r, [](T &valF, const T &valR) -> bool { return _add(valF, valR); });
    }

    return g;
  }

  /**
   * @brief Modifies the first parameter such that an entry at index \f$ p \f$ is equal to \f$ g(p) + val \f$.
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
   * @param g First element of the addition.
   * @param val Second element of the addition.
   */
  friend Multi_parameter_generator &operator+=(Multi_parameter_generator &g, const T &val)
  {
    g._apply_operation(val, [](T &valF, const T &valR) -> bool { return _add(valF, valR); });
    return g;
  }

  // multiplication
  /**
   * @brief Returns a generator such that an entry at index \f$ p \f$, with
   * \f$ 0 \leq p \leq num_parameters \f$ is equal to \f$ g(p) * r(p) \f$
   * if \f$ p < length_r \f$ and to \f$ g(p) \f$ otherwise.
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
   * @tparam ValueRange Range with a begin(), end() and size() method.
   * @param g First element of the multiplication.
   * @param r Second element of the multiplication.
   */
  template <class ValueRange, class = std::enable_if_t<RangeTraits<ValueRange>::has_begin> >
  friend Multi_parameter_generator operator*(Multi_parameter_generator g, const ValueRange &r)
  {
    g *= r;
    return g;
  }

  /**
   * @brief Returns a generator such that an entry at index \f$ p \f$, with
   * \f$ 0 \leq p \leq num_parameters \f$ is equal to \f$ g(p) * r(p) \f$
   * if \f$ p < length_r \f$ and to \f$ g(p) \f$ otherwise.
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
   * @tparam ValueRange Range with a begin(), end() and size() method.
   * @param r First element of the multiplication.
   * @param g Second element of the multiplication.
   */
  template <class ValueRange,
            class = std::enable_if_t<RangeTraits<ValueRange>::has_begin &&
                                     !std::is_same_v<ValueRange, Multi_parameter_generator> > >
  friend Multi_parameter_generator operator*(const ValueRange &r, Multi_parameter_generator g)
  {
    g *= r;
    return g;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ p \f$ is equal to \f$ g(p) * val \f$.
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
   * @param g First element of the multiplication.
   * @param val Second element of the multiplication.
   */
  friend Multi_parameter_generator operator*(Multi_parameter_generator g, const T &val)
  {
    g *= val;
    return g;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ p \f$ is equal to \f$ g(p) * val \f$.
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
   * @param g Second element of the multiplication.
   */
  friend Multi_parameter_generator operator*(const T &val, Multi_parameter_generator g)
  {
    g *= val;
    return g;
  }

  /**
   * @brief Modifies the first argument such that an entry at index \f$ p \f$, with
   * \f$ 0 \leq p \leq num_parameters \f$ is equal to \f$ g(p) * r(p) \f$
   * if \f$ p < length_r \f$ and to \f$ g(p) \f$ otherwise.
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
   * @tparam ValueRange Range with a begin(), end() and size() method.
   * @param g First element of the multiplication.
   * @param r Second element of the multiplication.
   */
  template <class ValueRange, class = std::enable_if_t<RangeTraits<ValueRange>::has_begin> >
  friend Multi_parameter_generator &operator*=(Multi_parameter_generator &g, const ValueRange &r)
  {
    if (g.num_parameters() == 0 || g.is_nan()) return g;

    if constexpr (RangeTraits<ValueRange>::is_dynamic_multi_filtration) {
      if (r[0].is_finite()) {
        return g *= r[0];
      }
      if constexpr (!std::numeric_limits<T>::has_quiet_NaN)
        if (r[0].size() == 0) {
          g = nan();
          return g;
        }
      return g *= r[0][0];
    } else {
      if (!g.is_finite()) {
        g.generator_.resize(r.size(), g[0]);
      }
      g._apply_operation(r, [](T &valF, const T &valR) -> bool { return _multiply(valF, valR); });
    }

    return g;
  }

  /**
   * @brief Modifies the first parameter such that an entry at index \f$ p \f$ is equal to \f$ g(p) * val \f$.
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
   * @param g First element of the multiplication.
   * @param val Second element of the multiplication.
   */
  friend Multi_parameter_generator &operator*=(Multi_parameter_generator &g, const T &val)
  {
    g._apply_operation(val, [](T &valF, const T &valR) -> bool { return _multiply(valF, valR); });
    return g;
  }

  // division
  /**
   * @brief Returns a generator such that an entry at index \f$ p \f$, with
   * \f$ 0 \leq p \leq num_parameters \f$ is equal to \f$ g(p) / r(p) \f$
   * if \f$ p < length_r \f$ and to \f$ g(p) \f$ otherwise.
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
   * @tparam ValueRange Range with a begin(), end() and size() method.
   * @param g First element of the division.
   * @param r Second element of the division.
   */
  template <class ValueRange, class = std::enable_if_t<RangeTraits<ValueRange>::has_begin> >
  friend Multi_parameter_generator operator/(Multi_parameter_generator g, const ValueRange &r)
  {
    g /= r;
    return g;
  }

  /**
   * @brief Returns a generator such that an entry at index \f$ p \f$, with
   * \f$ 0 \leq p \leq num_parameters \f$ is equal to \f$ r(p) / g(p) \f$
   * if \f$ p < length_r \f$ and to \f$ g(p) \f$ otherwise.
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
   * @tparam ValueRange Range with a begin(), end() and size() method.
   * @param r First element of the division.
   * @param g Second element of the division.
   */
  template <class ValueRange,
            class = std::enable_if_t<RangeTraits<ValueRange>::has_begin &&
                                     !std::is_same_v<ValueRange, Multi_parameter_generator> > >
  friend Multi_parameter_generator operator/(const ValueRange &r, const Multi_parameter_generator &g)
  {
    if (g.num_parameters() == 0) return g;
    if (g.is_nan()) return nan();

    if constexpr (RangeTraits<ValueRange>::is_dynamic_multi_filtration) {
      if (r[0].is_finite()) {
        return r[0] / g;
      }
      if constexpr (!std::numeric_limits<T>::has_quiet_NaN)
        if (r[0].size() == 0) return nan();
      return r[0][0] / g;
    } else {
      if (g.is_finite()) {
        Multi_parameter_generator res = g;
        res._apply_operation(r, [](T &valF, const T &valR) -> bool {
          T tmp = valF;
          valF = valR;
          return _divide(valF, tmp);
        });
        return res;
      }
      Multi_parameter_generator res(r.begin(), r.end());
      res._apply_operation(g[0], [](T &valRes, const T &valF) -> bool { return _divide(valRes, valF); });
      return res;
    }
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ p \f$ is equal to \f$ g(p) / val \f$.
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
   * @param g First element of the division.
   * @param val Second element of the division.
   */
  friend Multi_parameter_generator operator/(Multi_parameter_generator g, const T &val)
  {
    g /= val;
    return g;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ p \f$ is equal to \f$ val / g(p) \f$.
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
   * @param g Second element of the division.
   */
  friend Multi_parameter_generator operator/(const T &val, Multi_parameter_generator g)
  {
    g._apply_operation(val, [](T &valF, const T &valR) -> bool {
      T tmp = valF;
      valF = valR;
      return _divide(valF, tmp);
    });
    return g;
  }

  /**
   * @brief Modifies the first argument such that an entry at index \f$ p \f$, with
   * \f$ 0 \leq p \leq num_parameters \f$ is equal to \f$ g(p) / r(p) \f$
   * if \f$ p < length_r \f$ and to \f$ g(p) \f$ otherwise.
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
   * @tparam ValueRange Range with a begin(), end() and size() method.
   * @param g First element of the division.
   * @param r Second element of the division.
   */
  template <class ValueRange, class = std::enable_if_t<RangeTraits<ValueRange>::has_begin> >
  friend Multi_parameter_generator &operator/=(Multi_parameter_generator &g, const ValueRange &r)
  {
    if (g.num_parameters() == 0 || g.is_nan()) return g;

    if constexpr (RangeTraits<ValueRange>::is_dynamic_multi_filtration) {
      if (r[0].is_finite()) {
        return g /= r[0];
      }
      if constexpr (!std::numeric_limits<T>::has_quiet_NaN)
        if (r[0].size() == 0) {
          g = nan();
          return g;
        }
      return g /= r[0][0];
    } else {
      if (!g.is_finite()) {
        g.generator_.resize(r.size(), g[0]);
      }
      g._apply_operation(r, [](T &valF, const T &valR) -> bool { return _divide(valF, valR); });
    }

    return g;
  }

  /**
   * @brief Modifies the first parameter such that an entry at index \f$ p \f$ is equal to \f$ g(p) / val \f$.
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
   * @param g First element of the division.
   * @param val Second element of the division.
   */
  friend Multi_parameter_generator &operator/=(Multi_parameter_generator &g, const T &val)
  {
    g._apply_operation(val, [](T &valF, const T &valR) -> bool { return _divide(valF, valR); });
    return g;
  }

  // MODIFIERS

  /**
   * @brief If the filtration value is at +/- infinity or NaN, the underlying container will probably only contain
   * one element, representing the value. This method forces the underlying container to contain explicitly
   * `number_of_parameters` elements. If the container is empty, it will fill all new elements with the default
   * value -inf. If the container had more than one element, it does nothing and fails in Debug Mode if the number was
   * different from `number_of_parameters`.
   */
  void force_size_to_number_of_parameters(int number_of_parameters)
  {
    if (number_of_parameters < 1) return;

    if (generator_.size() > 1) {
      GUDHI_CHECK(static_cast<std::size_t>(number_of_parameters) == generator_.size(),
                  std::invalid_argument("Cannot force size to another number of parameters than set."));
      return;
    }

    auto val = generator_.empty() ? _get_default_value() : generator_[0];
    generator_.resize(number_of_parameters, val);
  }

  /**
   * @brief Sets the generator to the least common upper bound between it and the given value.
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
    if (is_nan() || is_plus_inf()) return false;
    if (x.size() == 0) {
      *this = nan();
      return true;
    }

    if constexpr (RangeTraits<GeneratorRange>::is_dynamic_multi_filtration) {
      return push_to_least_common_upper_bound(x[0], exclude_infinite_values);
    } else {
      if (is_minus_inf()) {
        generator_ = Underlying_container(x.begin(), x.end());
        if (is_minus_inf()) {
          *this = minus_inf();
          return false;
        }
        if (is_nan()) *this = nan();
        if (is_plus_inf()) *this = inf();
        return true;
      }

      bool modified = false;

      auto it = x.begin();
      for (size_type p = 0; p < generator_.size(); ++p) {
        T valX = *it;
        if (!exclude_infinite_values || (valX != T_inf && valX != T_m_inf)) {
          modified |= generator_[p] < valX;
          generator_[p] = valX > generator_[p] ? valX : generator_[p];
        }
        if (it + 1 != x.end()) ++it;
      }

      return modified;
    }
  }

  /**
   * @brief Sets each the generator to the greatest common lower bound between it and the given value.
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
    if (is_nan() || is_minus_inf()) return false;
    if (x.size() == 0) {
      *this = nan();
      return true;
    }

    if constexpr (RangeTraits<GeneratorRange>::is_dynamic_multi_filtration) {
      return pull_to_greatest_common_lower_bound(x[0], exclude_infinite_values);
    } else {
      if (is_plus_inf()) {
        generator_ = Underlying_container(x.begin(), x.end());
        if (is_plus_inf()) {
          *this = inf();
          return false;
        }
        if (is_nan()) *this = nan();
        if (is_minus_inf()) *this = minus_inf();
        return true;
      }

      bool modified = false;

      auto it = x.begin();
      for (size_type p = 0; p < generator_.size() && it != x.end(); ++p) {
        T valX = *it;
        if (!exclude_infinite_values || (valX != T_inf && valX != T_m_inf)) {
          modified |= generator_[p] > valX;
          generator_[p] = valX < generator_[p] ? valX : generator_[p];
        }
        if (it + 1 != x.end()) ++it;
      }

      return modified;
    }
  }

  /**
   * @brief Projects the generator into the given grid. If @p coordinate is false, the entries are set to
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
        grid.size() >= generator_.size(),
        std::invalid_argument("The grid should not be smaller than the number of parameters in the filtration value."));

    if (is_nan() || generator_.empty()) return;

    if (!is_finite()) generator_.resize(grid.size(), generator_[0]);

    for (size_type p = 0; p < generator_.size(); ++p) {
      const auto &filtration = grid[p];
      auto v = static_cast<typename OneDimArray::value_type>(generator_[p]);
      auto d = std::distance(filtration.begin(), std::lower_bound(filtration.begin(), filtration.end(), v));
      if (d != 0 && std::abs(v - filtration[d]) > std::abs(v - filtration[d - 1])) {
        --d;
      }
      generator_[p] = coordinate ? static_cast<T>(d) : static_cast<T>(filtration[d]);
    }
  }

  // FONCTIONNALITIES

  /**
   * @brief Computes the scalar product of the generator with the given vector.
   *
   * @tparam U Arithmetic type of the result. Default value: `T`.
   * @param g Filtration value.
   * @param x Vector of coefficients.
   * @return Scalar product of @p f with @p x.
   */
  template <typename U = T>
  friend U compute_linear_projection(const Multi_parameter_generator &g, const std::vector<U> &x)
  {
    U projection = 0;
    std::size_t size = std::min(x.size(), g.num_parameters());
    for (std::size_t i = 0; i < size; i++) projection += x[i] * static_cast<U>(g[i]);
    return projection;
  }

  /**
   * @brief Computes the euclidean distance from the first parameter to the second parameter.
   *
   * @param g Source filtration value.
   * @param other Target filtration value.
   * @return Euclidean distance between @p f and @p other.
   */
  template <typename U = T>
  friend U compute_euclidean_distance_to(const Multi_parameter_generator &g, const Multi_parameter_generator &other)
  {
    if (!g.is_finite() || !other.is_finite()) {
      if constexpr (std::numeric_limits<T>::has_quiet_NaN)
        return std::numeric_limits<T>::quiet_NaN();
      else
        return T_inf;
    }

    GUDHI_CHECK(g.num_parameters() == other.num_parameters(),
                std::invalid_argument("We cannot compute the distance between two points of different dimensions."));

    if (g.num_parameters() == 1) return g[0] - other[0];

    U out = 0;
    for (size_type p = 0; p < g.num_parameters(); ++p) {
      T v = g[p] - other[p];
      out += v * v;
    }
    if constexpr (std::is_integral_v<U>) {
      // to avoid Windows issue that don't know how to cast integers for cmath methods
      return std::sqrt(static_cast<double>(out));
    } else {
      return std::sqrt(out);
    }
  }

  /**
   * @brief Computes the sum of the squares of all values in the given generator.
   *
   * @param g Filtration value.
   * @return The norm of @p f.
   */
  template <typename U = T>
  friend U compute_squares(const Multi_parameter_generator &g)
  {
    U out = 0;
    for (size_type p = 0; p < g.num_parameters(); ++p) {
      out += g[p] * g[p];
    }
    return out;
  }

  /**
   * @brief Computes the values in the given grid corresponding to the coordinates given by the given generator.
   * That is, if \f$ out \f$ is the result, \f$ out(g,p) = grid[p][f(g,p)] \f$. Assumes therefore, that the
   * values stored in the filtration value corresponds to indices existing in the given grid.
   *
   * @tparam U Signed arithmetic type.
   * @param g Filtration value storing coordinates compatible with `grid`.
   * @param grid Vector of vector.
   * @return Filtration value \f$ out \f$ whose entry correspond to \f$ out(g,p) = grid[p][f(g,p)] \f$.
   */
  template <typename U>
  friend Multi_parameter_generator<U> evaluate_coordinates_in_grid(const Multi_parameter_generator &g,
                                                                   const std::vector<std::vector<U> > &grid)
  {
    GUDHI_CHECK(grid.size() >= g.num_parameters(),
                std::invalid_argument(
                    "The size of the grid should correspond to the number of parameters in the filtration value."));

    U grid_inf = Multi_parameter_generator<U>::T_inf;
    std::vector<U> outVec(g.num_parameters());

    for (size_type p = 0; p < g.num_parameters(); ++p) {
      outVec[p] = (g[p] == T_inf ? grid_inf : grid[p][g[p]]);
    }

    return Multi_parameter_generator<U>(std::move(outVec));
  }

  // UTILITIES

  /**
   * @brief Outstream operator.
   */
  friend std::ostream &operator<<(std::ostream &stream, const Multi_parameter_generator &g)
  {
    if (g.is_plus_inf()) {
      stream << "[inf, ..., inf]";
      return stream;
    }
    if (g.is_minus_inf()) {
      stream << "[-inf, ..., -inf]";
      return stream;
    }
    if (g.is_nan()) {
      stream << "[NaN]";
      return stream;
    }

    const size_type num_param = g.num_parameters();

    stream << "[";
    for (size_type p = 0; p < num_param; ++p) {
      stream << g[p];
      if (p < num_param - 1) stream << ", ";
    }
    stream << "]";

    return stream;
  }

  /**
   * @brief Instream operator.
   */
  friend std::istream &operator>>(std::istream &stream, Multi_parameter_generator &g)
  {
    char firstCharacter;
    stream >> firstCharacter;
    if (firstCharacter != '[')
      throw std::invalid_argument("Invalid incoming stream format for Multi_parameter_generator.");
    g.generator_.clear();
    auto pos = stream.tellg();
    stream >> firstCharacter;
    if (firstCharacter == ']') return stream;
    if (firstCharacter == 'i') {
      while (firstCharacter != ']') {
        stream >> firstCharacter;
      }
      g = Multi_parameter_generator::inf();
      return stream;
    }
    if (firstCharacter == '-') {
      stream >> firstCharacter;
      if (firstCharacter == 'i') {
        while (firstCharacter != ']') {
          stream >> firstCharacter;
        }
        g = Multi_parameter_generator::minus_inf();
        return stream;
      }  // else could be a negative number
    }
    if (firstCharacter == 'N') {
      while (firstCharacter != ']') {
        stream >> firstCharacter;
      }
      g = Multi_parameter_generator::nan();
      return stream;
    }

    stream.seekg(pos, std::ios_base::beg);
    char delimiter = '\0';
    while (delimiter != ']') {
      g.generator_.push_back(_get_value<T>(stream));
      if (!stream.good()) throw std::invalid_argument("Invalid incoming stream format for Multi_parameter_generator.");
      stream >> delimiter;
    }

    return stream;
  }

  /**
   * @brief Serialize given value into the buffer at given pointer.
   *
   * @param value Value to serialize.
   * @param start Pointer to the start of the space in the buffer where to store the serialization.
   * @return End position of the serialization in the buffer.
   */
  friend char *serialize_value_to_char_buffer(const Multi_parameter_generator &value, char *start)
  {
    const size_type length = value.generator_.size();
    const std::size_t arg_size = sizeof(T) * length;
    const std::size_t type_size = sizeof(size_type);
    memcpy(start, &length, type_size);
    memcpy(start + type_size, value.generator_.data(), arg_size);
    return start + arg_size + type_size;
  }

  /**
   * @brief Deserialize the value from a buffer at given pointer and stores it in given value.
   *
   * @param value Value to fill with the deserialized filtration value.
   * @param start Pointer to the start of the space in the buffer where the serialization is stored.
   * @return End position of the serialization in the buffer.
   */
  friend const char *deserialize_value_from_char_buffer(Multi_parameter_generator &value, const char *start)
  {
    const std::size_t type_size = sizeof(size_type);
    size_type length;
    memcpy(&length, start, type_size);
    std::size_t arg_size = sizeof(T) * length;
    value.generator_.resize(length);
    memcpy(value.generator_.data(), start + type_size, arg_size);
    return start + arg_size + type_size;
  }

  /**
   * @brief Returns the serialization size of the given filtration value.
   */
  friend std::size_t get_serialization_size_of(const Multi_parameter_generator &value)
  {
    return sizeof(size_type) + (sizeof(T) * value.generator_.size());
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
  Underlying_container generator_; /**< Container of the filtration value elements. */

  /**
   * @brief Default value of an element in the filtration value.
   */
  constexpr static T _get_default_value() { return T_m_inf; }

  /**
   * @brief Applies operation on the elements of the filtration value.
   */
  template <class ValueRange, class F>
  void _apply_operation(const ValueRange &range, F &&operate)
  {
    if (range.size() < generator_.size()) {
      auto it = range.begin();
      for (size_type p = 0; p < generator_.size() && it != range.end(); ++p) {
        std::forward<F>(operate)(generator_[p], *it);
        ++it;
      }
    } else {
      bool allPlusInf = true;
      bool allMinusInf = true;
      bool allNaN = true;

      auto it = range.begin();
      for (size_type p = 0; p < generator_.size() && it != range.end(); ++p) {
        bool isNotNaN = std::forward<F>(operate)(generator_[p], *it);
        if (generator_[p] != T_inf) allPlusInf = false;
        if (generator_[p] != T_m_inf) allMinusInf = false;
        if (isNotNaN) allNaN = false;
        ++it;
      }

      if (allPlusInf)
        *this = inf();
      else if (allMinusInf)
        *this = minus_inf();
      else if (allNaN)
        *this = nan();
    }
  }

  /**
   * @brief Applies operation on the elements of the filtration value.
   */
  template <class F>
  void _apply_operation(const T &val, F &&operate)
  {
    bool allNaN = true;

    for (auto &p : generator_) {
      if (std::forward<F>(operate)(p, val)) allNaN = false;
    }

    if (allNaN) *this = nan();
  }
};

}  // namespace Gudhi::multi_filtration

namespace std {

template <typename T>
class numeric_limits<Gudhi::multi_filtration::Multi_parameter_generator<T> >
{
 public:
  using Generator = Gudhi::multi_filtration::Multi_parameter_generator<T>;

  static constexpr bool has_infinity = true;
  static constexpr bool has_quiet_NaN = true;

  static constexpr Generator infinity() noexcept { return Generator::inf(); };

  // non-standard
  static constexpr Generator minus_infinity() noexcept { return Generator::minus_inf(); };

  static constexpr Generator max() noexcept(false)
  {
    throw std::logic_error(
        "The max value cannot be represented with no finite numbers of parameters."
        "Use `max(number_of_parameters)` instead");
  };

  static constexpr Generator max(std::size_t p) noexcept { return Generator(p, std::numeric_limits<T>::max()); };

  static constexpr Generator quiet_NaN() noexcept { return Generator::nan(); };
};

}  // namespace std

#endif  // MF_MULTI_PARAMETER_GENERATOR_H_
