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

#include <algorithm>        //std::lower_bound
#include <cmath>            //std::isnan, std::min
#include <cstddef>          //std::size_t
#include <cstdint>          //std::int32_t
#include <cstring>          //memcpy
#include <ostream>          //std::ostream
#include <limits>           //std::numerical_limits
#include <type_traits>      //std::is_arithmetic
#include <utility>          //std::swap, std::move
#include <numeric>          //std::iota
#include <vector>
#include <initializer_list>

#include <gudhi/Debug_utils.h>
#include <gudhi/Simple_mdspan.h>

namespace Gudhi::multi_filtration {

// declaration needed pre C++20 for friends with templates defined inside a class
template <typename U>
U compute_euclidean_distance_to();
template <typename U>
U compute_norm();

template <typename T, bool co = false>
class Multi_parameter_filtration
{
 public:
  using Underlying_container = std::vector<T>;

  // CONSTRUCTORS

  Multi_parameter_filtration(int number_of_parameters)
      : generators_(number_of_parameters < 0 ? 2 : number_of_parameters, _get_default_value()),
        generator_view_(generators_.data(), generators_.empty() ? 0 : 1, generators_.size())
  {}

  Multi_parameter_filtration(int number_of_parameters, T value)
      : generators_(number_of_parameters < 0 ? 2 : number_of_parameters, value),
        generator_view_(generators_.data(), generators_.empty() ? 0 : 1, generators_.size())
  {}

  template <class ValueRange = std::initializer_list<T> >
  Multi_parameter_filtration(const ValueRange &range)
      : generators_(range.begin(), range.end()),
        generator_view_(generators_.data(), generators_.empty() ? 0 : 1, generators_.size())
  {}

  template <class Iterator>
  Multi_parameter_filtration(Iterator it_begin, Iterator it_end)
      : generators_(it_begin, it_end),
        generator_view_(generators_.data(), generators_.empty() ? 0 : 1, generators_.size())
  {}

  Multi_parameter_filtration(const Underlying_container &generators, int number_of_parameters)
      : generators_(generators),
        generator_view_(generators_.data(),
                        number_of_parameters == 0 ? 0 : generators_.size() / number_of_parameters,
                        number_of_parameters)
  {
    GUDHI_CHECK(number_of_parameters != 0 || generators_.empty(),
                "Number of parameters cannot be 0 if the container is not empty.");
  }

  Multi_parameter_filtration(Underlying_container &&generators, int number_of_parameters)
      : generators_(std::move(generators)),
        generator_view_(generators_.data(),
                        number_of_parameters == 0 ? 0 : generators_.size() / number_of_parameters,
                        number_of_parameters)
  {
    GUDHI_CHECK(number_of_parameters != 0 || generators_.empty(),
                "Number of parameters cannot be 0 if the container is not empty.");
  }

  // VECTOR-LIKE

  using value_type = T;
  using size_type = typename Underlying_container::size_type;
  using difference_type = typename Underlying_container::difference_type;
  using reference = value_type &;
  using const_reference = const value_type &;
  using pointer = typename Underlying_container::pointer;
  using const_pointer = typename Underlying_container::const_pointer;

  reference operator()(size_type g, size_type n) { return generator_view_(g, n); }

  const_reference operator()(size_type g, size_type n) const { return generator_view_(g, n); }

  template <class IndexRange = std::initializer_list<size_type> >
  reference operator[](const IndexRange &indices)
  {
    GUDHI_CHECK(indices.size() == 2,
                "Exactly 2 indices allowed only: first the generator number, second the parameter number.");
    auto it = indices.begin();
    return generator_view_(*it, *(++it));
  }

  template <class IndexRange = std::initializer_list<size_type> >
  const_reference operator[](const IndexRange &indices) const
  {
    GUDHI_CHECK(indices.size() == 2,
                "Exactly 2 indices allowed only: first the generator number, second the parameter number.");
    auto it = indices.begin();
    return generator_view_(*it, *(++it));
  }

  /**
   * @brief Reserves space for the given number of generators in the underlying container.
   *
   * @param n Number of generators.
   */
  void reserve(size_type number_of_generators) { generators_.reserve(num_parameters() * number_of_generators); }

  // CONVERTERS

  // like numpy
  /**
   * @brief Returns a copy with entries casted into the type given as template parameter.
   *
   * @tparam U New type for the entries.
   * @return Copy with new entry type.
   */
  template <typename U>
  Multi_parameter_filtration<U, co> as_type() const
  {
    std::vector<U> out(generators_.begin(), generators_.end());
    return Multi_parameter_filtration<U, co>(std::move(out), num_parameters());
  }

  // ACCESS

  /**
   * @brief Returns the number of parameters of the finite filtration value. If the value is "inf", "-inf" or "NaN",
   * returns 1.
   *
   * @return Number of parameters.
   */
  size_type num_parameters() const { return generator_view_.extent(1); }

  size_type num_generators() const { return generator_view_.extent(0); }

  size_type num_entries() const { return generators_.size(); }

  /**
   * @brief Returns a filtration value for which @ref is_plus_inf() returns `true`.
   *
   * @return Infinity.
   */
  static Multi_parameter_filtration inf(int number_of_parameters)
  {
    return Multi_parameter_filtration(number_of_parameters, T_inf);
  }

  /**
   * @brief Returns a filtration value for which @ref is_minus_inf() returns `true`.
   *
   * @return Minus infinity.
   */
  static Multi_parameter_filtration minus_inf(int number_of_parameters)
  {
    return Multi_parameter_filtration(number_of_parameters, -T_inf);
  }

  /**
   * @brief Returns a filtration value for which @ref is_nan() returns `true`. If `T` does not support NaN values,
   * all values in the array will be 0 and will not be recognized as NaN.
   *
   * @return NaN.
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
   * @brief Returns `true` if and only if the filtration value is considered as infinity.
   */
  bool is_plus_inf() const
  {
    for (const T &v : generators_) {
      if (v != T_inf) return false;
    }
    return true;
  }

  /**
   * @brief Returns `true` if and only if the filtration value is considered as minus infinity.
   */
  bool is_minus_inf() const
  {
    for (const T &v : generators_) {
      if (v != -T_inf) return false;
    }
    return true;
  }

  /**
   * @brief Returns `true` if and only if the filtration value is considered as NaN.
   */
  bool is_nan() const
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
   * @brief Returns `true` if and only if the filtration value is non-empty and is not considered as infinity,
   * minus infinity or NaN.
   */
  bool is_finite() const
  {
    bool isInf = true, isMinusInf = true, isNan = true;
    for (const auto &v : generators_) {
      if (v != T_inf) isInf = false;
      if (v != -T_inf) isMinusInf = false;
      if (!_is_nan(v)) isNan = false;
      if (!isInf && !isMinusInf && !isNan) return true;
    }
    return false;
  }

  // COMPARAISON OPERATORS

  /**
   * @brief Returns `true` if and only if for each \f$ i \f$, \f$ a[i] \f$ is strictly smaller than \f$ b[i] \f$.
   * If @p a and @p b are both not infinite or NaN, they have to have the same number of parameters.
   *
   * Note that not all filtration values are comparable. That is, \f$ a < b \f$ and \f$ b < a \f$ returning both false
   * does **not** imply \f$ a == b \f$.
   */
  friend bool operator<(const Multi_parameter_filtration &a, const Multi_parameter_filtration &b)
  {
    GUDHI_CHECK(a.num_parameters() == b.num_parameters(),
                "Only filtration values with same number of parameters can be compared.");

    auto view_a = a.generator_view_;
    auto view_b = b.generator_view_;
    for (std::size_t i = 0u; i < b.num_generators(); ++i) {
      // for each generator in b, verify if it is strictly in the cone of at least one generator of a
      bool isContained = false;
      for (std::size_t j = 0u; j < a.num_generators() && !isContained; ++j) {
        // lexicographical order, so if a[j][0] dom b[j][0], than a[j'] can never strictly contain b[i] for all j' > j.
        if (_first_dominates(view_a, j, view_b, i)) return false;
        isContained = _strictly_contains(view_a, j, view_b, i);
      }
      if (!isContained) return false;
    }
    return true;
  }

  /**
   * @brief Returns `true` if and only if for each \f$ i \f$, \f$ a[i] \f$ is smaller or equal than \f$ b[i] \f$.
   * If @p a and @p b are both not infinite or NaN, they have to have the same number of parameters.
   *
   * Note that not all filtration values are comparable. That is, \f$ a \le b \f$ and \f$ b \le a \f$ can both return
   * `false`.
   */
  friend bool operator<=(const Multi_parameter_filtration &a, const Multi_parameter_filtration &b)
  {
    GUDHI_CHECK(a.num_parameters() == b.num_parameters(),
                "Only filtration values with same number of parameters can be compared.");

    auto view_a = a.generator_view_;
    auto view_b = b.generator_view_;
    // check if this curves is below other's curve
    //  ie for each guy in this, check if there is a guy in other that dominates him
    for (std::size_t i = 0u; i < b.num_generators(); ++i) {
      // for each generator in b, verify if it is in the cone of at least one generator of a
      bool isContained = false;
      for (std::size_t j = 0u; j < a.num_generators() && !isContained; ++j) {
        // lexicographical order, so if a[j][0] strictly dom b[j][0], than a[j'] can never contain b[i] for all j' > j.
        if (_first_strictly_dominates(view_a, j, view_b, i)) return false;
        isContained = _contains(view_a, j, view_b, i);
      }
      if (!isContained) return false;
    }
    return true;
  }

  /**
   * @brief Returns `true` if and only if for each \f$ i \f$, \f$ a[i] \f$ is strictly greater than \f$ b[i] \f$.
   * If @p a and @p b are both not infinite or NaN, they have to have the same number of parameters.
   *
   * Note that not all filtration values are comparable. That is, \f$ a > b \f$ and \f$ b > a \f$ returning both false
   * does **not** imply \f$ a == b \f$.
   */
  friend bool operator>(const Multi_parameter_filtration &a, const Multi_parameter_filtration &b) { return b < a; }

  /**
   * @brief Returns `true` if and only if for each \f$ i \f$, \f$ a[i] \f$ is greater or equal than \f$ b[i] \f$.
   * If @p a and @p b are both not infinite or NaN, they have to have the same number of parameters.
   *
   * Note that not all filtration values are comparable. That is, \f$ a \ge b \f$ and \f$ b \ge a \f$ can both return
   * `false`.
   */
  friend bool operator>=(const Multi_parameter_filtration &a, const Multi_parameter_filtration &b) { return b <= a; }

  /**
   * @brief Returns `true` if and only if for each \f$ i \f$, \f$ a[i] \f$ is equal to \f$ b[i] \f$.
   */
  friend bool operator==(const Multi_parameter_filtration &a, const Multi_parameter_filtration &b)
  {
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
   * @brief Returns a filtration value such that an entry at index \f$ i \f$ is equal to \f$ -f[i] \f$.
   *
   * Used conventions:
   * - \f$ -NaN = NaN \f$.
   *
   * @param f Value to opposite.
   * @return The opposite of @p f.
   */
  friend Multi_parameter_filtration operator-(const Multi_parameter_filtration &f)
  {
    std::vector<T> result(f.generators_);
    std::for_each(result.begin(), result.end(), [](T &v) { v = -v; });
    return Multi_parameter_filtration(std::move(result), f.num_parameters());
  }

  // subtraction
  /**
   * @brief Returns a filtration value such that an entry at index \f$ i \f$ is equal to \f$ a[i] - b[i] \f$.
   * If @p a and @p b are both not infinite or NaN, they have to have the same number of parameters.
   *
   * Used conventions:
   * - \f$ inf - inf = NaN \f$,
   * - \f$ -inf - (-inf) = NaN \f$,
   * - \f$ NaN - b = NaN \f$,
   * - \f$ a - NaN = NaN \f$.
   *
   * If `std::numeric_limits<T>::has_quiet_NaN` is false, then the returned filtration value will be @ref nan()
   * if any operation results in NaN, not only if all operations result in NaN.
   *
   * @param a First element of the subtraction.
   * @param b Second element of the subtraction.
   * @return Entry-wise \f$ a - b \f$.
   */
  template <class ValueRange, class = std::enable_if_t<!std::is_arithmetic_v<ValueRange> > >
  friend Multi_parameter_filtration operator-(Multi_parameter_filtration f, const ValueRange &r)
  {
    f -= r;
    return f;
  }

  template <class ValueRange, class = std::enable_if_t<!std::is_arithmetic_v<ValueRange> > >
  friend Multi_parameter_filtration operator-(const ValueRange &r, Multi_parameter_filtration f)
  {
    return f._apply_operation(r, [](T &valF, const T &valR) {
      valF = -valF;
      _add(valF, valR);
    });
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ i \f$ is equal to \f$ f[i] - val \f$.
   *
   * Used conventions:
   * - \f$ inf - inf = NaN \f$,
   * - \f$ -inf - (-inf) = NaN \f$,
   * - \f$ NaN - b = NaN \f$,
   * - \f$ a - NaN = NaN \f$.
   *
   * If `std::numeric_limits<T>::has_quiet_NaN` is false, then the returned filtration value will be @ref nan()
   * if any operation results in NaN, not only if all operations result in NaN.
   *
   * @param f First element of the subtraction.
   * @param val Second element of the subtraction.
   * @return Entry-wise \f$ f - val \f$.
   */
  friend Multi_parameter_filtration operator-(Multi_parameter_filtration f, const T &val)
  {
    f -= val;
    return f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ i \f$ is equal to \f$ val - f[i] \f$.
   *
   * Used conventions:
   * - \f$ inf - inf = NaN \f$,
   * - \f$ -inf - (-inf) = NaN \f$,
   * - \f$ NaN - b = NaN \f$,
   * - \f$ a - NaN = NaN \f$.
   *
   * If `std::numeric_limits<T>::has_quiet_NaN` is false, then the returned filtration value will be @ref nan()
   * if any operation results in NaN, not only if all operations result in NaN.
   *
   * @param val First element of the subtraction.
   * @param f Second element of the subtraction.
   * @return Entry-wise \f$ val - f \f$.
   */
  friend Multi_parameter_filtration operator-(const T &val, Multi_parameter_filtration f)
  {
    return f._apply_operation(val, [](T &valF, const T &valR) {
      valF = -valF;
      _add(valF, valR);
    });
  }

  /**
   * @brief Modifies the first parameters such that an entry at index \f$ i \f$ is equal to
   * \f$ result[i] - to\_subtract[i] \f$.
   * If @p result and @p to_subtract are both not infinite or NaN, they have to have the same number of parameters.
   *
   * Used conventions:
   * - \f$ inf - inf = NaN \f$,
   * - \f$ -inf - (-inf) = NaN \f$,
   * - \f$ NaN - b = NaN \f$,
   * - \f$ a - NaN = NaN \f$.
   *
   * If `std::numeric_limits<T>::has_quiet_NaN` is false, then the returned filtration value will be @ref nan()
   * if any operation results in NaN, not only if all operations result in NaN.
   *
   * @param result First element of the subtraction.
   * @param to_subtract Second element of the subtraction.
   * @return Entry-wise \f$ result - to\_subtract \f$.
   */
  template <class ValueRange, class = std::enable_if_t<!std::is_arithmetic_v<ValueRange> > >
  friend Multi_parameter_filtration &operator-=(Multi_parameter_filtration &f, const ValueRange &r)
  {
    return f._apply_operation(r, [](T &valF, const T &valR) { _subtract(valF, valR); });
  }

  /**
   * @brief Modifies the first parameters such that an entry at index \f$ i \f$ is equal to
   * \f$ result[i] - to\_subtract \f$.
   *
   * Used conventions:
   * - \f$ inf - inf = NaN \f$,
   * - \f$ -inf - (-inf) = NaN \f$,
   * - \f$ NaN - b = NaN \f$,
   * - \f$ a - NaN = NaN \f$.
   *
   * If `std::numeric_limits<T>::has_quiet_NaN` is false, then the returned filtration value will be @ref nan()
   * if any operation results in NaN, not only if all operations result in NaN.
   *
   * @param result First element of the subtraction.
   * @param to_subtract Second element of the subtraction.
   * @return Entry-wise \f$ result - to\_subtract \f$.
   */
  friend Multi_parameter_filtration &operator-=(Multi_parameter_filtration &f, const T &val)
  {
    return f._apply_operation(val, [](T &valF, const T &valR) { _subtract(valF, valR); });
  }

  // addition
  /**
   * @brief Returns a filtration value such that an entry at index \f$ i \f$ is equal to \f$ a[i] + b[i] \f$.
   * If @p a and @p b are both not infinite or NaN, they have to have the same number of parameters.
   *
   * Used conventions:
   * - \f$ NaN + b = NaN \f$,
   * - \f$ a + NaN = NaN \f$.
   *
   * If `std::numeric_limits<T>::has_quiet_NaN` is false, then the returned filtration value will be @ref nan()
   * if any operation results in NaN, not only if all operations result in NaN.
   *
   * @param a First element of the addition.
   * @param b Second element of the addition.
   * @return Entry-wise \f$ a + b \f$.
   */
  template <class ValueRange, class = std::enable_if_t<!std::is_arithmetic_v<ValueRange> > >
  friend Multi_parameter_filtration operator+(Multi_parameter_filtration f, const ValueRange &r)
  {
    f += r;
    return f;
  }

  template <class ValueRange, class = std::enable_if_t<!std::is_arithmetic_v<ValueRange> > >
  friend Multi_parameter_filtration operator+(const ValueRange &r, Multi_parameter_filtration f)
  {
    f += r;
    return f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ i \f$ is equal to \f$ f[i] + val \f$.
   *
   * Used conventions:
   * - \f$ NaN + b = NaN \f$,
   * - \f$ a + NaN = NaN \f$.
   *
   * If `std::numeric_limits<T>::has_quiet_NaN` is false, then the returned filtration value will be @ref nan()
   * if any operation results in NaN, not only if all operations result in NaN.
   *
   * @param f First element of the addition.
   * @param val Second element of the addition.
   * @return Entry-wise \f$ f + val \f$.
   */
  friend Multi_parameter_filtration operator+(Multi_parameter_filtration f, const T &val)
  {
    f += val;
    return f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ i \f$ is equal to \f$ val + f[i] \f$.
   *
   * Used conventions:
   * - \f$ NaN + b = NaN \f$,
   * - \f$ a + NaN = NaN \f$.
   *
   * If `std::numeric_limits<T>::has_quiet_NaN` is false, then the returned filtration value will be @ref nan()
   * if any operation results in NaN, not only if all operations result in NaN.
   *
   * @param val First element of the addition.
   * @param f Second element of the addition.
   * @return Entry-wise \f$ val + f \f$.
   */
  friend Multi_parameter_filtration operator+(const T &val, Multi_parameter_filtration f)
  {
    f += val;
    return f;
  }

  /**
   * @brief Modifies the first parameters such that an entry at index \f$ i \f$ is equal to
   * \f$ result[i] + to\_add[i] \f$.
   * If @p result and @p to_add are both not infinite or NaN, they have to have the same number of parameters.
   *
   * Used conventions:
   * - \f$ NaN + b = NaN \f$,
   * - \f$ a + NaN = NaN \f$.
   *
   * If `std::numeric_limits<T>::has_quiet_NaN` is false, then the returned filtration value will be @ref nan()
   * if any operation results in NaN, not only if all operations result in NaN.
   *
   * @param result First element of the addition.
   * @param to_add Second element of the addition.
   * @return Entry-wise \f$ result + to\_add \f$.
   */
  template <class ValueRange, class = std::enable_if_t<!std::is_arithmetic_v<ValueRange> > >
  friend Multi_parameter_filtration &operator+=(Multi_parameter_filtration &f, const ValueRange &r)
  {
    return f._apply_operation(r, [](T &valF, const T &valR) { _add(valF, valR); });
  }

  /**
   * @brief Modifies the first parameters such that an entry at index \f$ i \f$ is equal to
   * \f$ result[i] + to\_add \f$.
   *
   * Used conventions:
   * - \f$ NaN + b = NaN \f$,
   * - \f$ a + NaN = NaN \f$.
   *
   * If `std::numeric_limits<T>::has_quiet_NaN` is false, then the returned filtration value will be @ref nan()
   * if any operation results in NaN, not only if all operations result in NaN.
   *
   * @param result First element of the addition.
   * @param to_add Second element of the addition.
   * @return Entry-wise \f$ result + to\_add \f$.
   */
  friend Multi_parameter_filtration &operator+=(Multi_parameter_filtration &f, const T &val)
  {
    return f._apply_operation(val, [](T &valF, const T &valR) { _add(valF, valR); });
  }

  // multiplication
  /**
   * @brief Returns a filtration value such that an entry at index \f$ i \f$ is equal to \f$ a[i] * b[i] \f$.
   * If @p a and @p b are both not infinite or NaN, they have to have the same number of parameters.
   *
   * Used conventions:
   * - \f$ inf * 0 = NaN \f$,
   * - \f$ 0 * inf = NaN \f$,
   * - \f$ -inf * 0 = NaN \f$,
   * - \f$ 0 * -inf = NaN \f$,
   * - \f$ NaN * b = NaN \f$,
   * - \f$ a * NaN = NaN \f$.
   *
   * If `std::numeric_limits<T>::has_quiet_NaN` is false, then the returned filtration value will be @ref nan()
   * if any operation results in NaN, not only if all operations result in NaN.
   *
   * @param a First element of the multiplication.
   * @param b Second element of the multiplication.
   * @return Entry-wise \f$ a * b \f$.
   */
  template <class ValueRange, class = std::enable_if_t<!std::is_arithmetic_v<ValueRange> > >
  friend Multi_parameter_filtration operator*(Multi_parameter_filtration f, const ValueRange &r)
  {
    f *= r;
    return f;
  }

  template <class ValueRange, class = std::enable_if_t<!std::is_arithmetic_v<ValueRange> > >
  friend Multi_parameter_filtration operator*(const ValueRange &r, Multi_parameter_filtration f)
  {
    f *= r;
    return f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ i \f$ is equal to \f$ f[i] * val \f$.
   *
   * Used conventions:
   * - \f$ inf * 0 = NaN \f$,
   * - \f$ 0 * inf = NaN \f$,
   * - \f$ -inf * 0 = NaN \f$,
   * - \f$ 0 * -inf = NaN \f$,
   * - \f$ NaN * b = NaN \f$,
   * - \f$ a * NaN = NaN \f$.
   *
   * If `std::numeric_limits<T>::has_quiet_NaN` is false, then the returned filtration value will be @ref nan()
   * if any operation results in NaN, not only if all operations result in NaN.
   *
   * @param f First element of the multiplication.
   * @param val Second element of the multiplication.
   * @return Entry-wise \f$ f * val \f$.
   */
  friend Multi_parameter_filtration operator*(Multi_parameter_filtration f, const T &val)
  {
    f *= val;
    return f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ i \f$ is equal to \f$ val * f[i] \f$.
   *
   * Used conventions:
   * - \f$ inf * 0 = NaN \f$,
   * - \f$ 0 * inf = NaN \f$,
   * - \f$ -inf * 0 = NaN \f$,
   * - \f$ 0 * -inf = NaN \f$,
   * - \f$ NaN * b = NaN \f$,
   * - \f$ a * NaN = NaN \f$.
   *
   * If `std::numeric_limits<T>::has_quiet_NaN` is false, then the returned filtration value will be @ref nan()
   * if any operation results in NaN, not only if all operations result in NaN.
   *
   * @param val First element of the multiplication.
   * @param f Second element of the multiplication.
   * @return Entry-wise \f$ val * f \f$.
   */
  friend Multi_parameter_filtration operator*(const T &val, Multi_parameter_filtration f)
  {
    f *= val;
    return f;
  }

  /**
   * @brief Modifies the first parameters such that an entry at index \f$ i \f$ is equal to
   * \f$ result[i] * to\_mul[i] \f$.
   * If @p result and @p to_mul are both not infinite or NaN, they have to have the same number of parameters.
   *
   * Used conventions:
   * - \f$ inf * 0 = NaN \f$,
   * - \f$ 0 * inf = NaN \f$,
   * - \f$ -inf * 0 = NaN \f$,
   * - \f$ 0 * -inf = NaN \f$,
   * - \f$ NaN * b = NaN \f$,
   * - \f$ a * NaN = NaN \f$.
   *
   * If `std::numeric_limits<T>::has_quiet_NaN` is false, then the returned filtration value will be @ref nan()
   * if any operation results in NaN, not only if all operations result in NaN.
   *
   * @param result First element of the multiplication.
   * @param to_mul Second element of the multiplication.
   * @return Entry-wise \f$ result * to\_mul \f$.
   */
  template <class ValueRange, class = std::enable_if_t<!std::is_arithmetic_v<ValueRange> > >
  friend Multi_parameter_filtration &operator*=(Multi_parameter_filtration &f, const ValueRange &r)
  {
    return f._apply_operation(r, [](T &valF, const T &valR) { _multiply(valF, valR); });
  }

  /**
   * @brief Modifies the first parameters such that an entry at index \f$ i \f$ is equal to
   * \f$ result[i] * to\_mul \f$.
   *
   * Used conventions:
   * - \f$ inf * 0 = NaN \f$,
   * - \f$ 0 * inf = NaN \f$,
   * - \f$ -inf * 0 = NaN \f$,
   * - \f$ 0 * -inf = NaN \f$,
   * - \f$ NaN * b = NaN \f$,
   * - \f$ a * NaN = NaN \f$.
   *
   * If `std::numeric_limits<T>::has_quiet_NaN` is false, then the returned filtration value will be @ref nan()
   * if any operation results in NaN, not only if all operations result in NaN.
   *
   * @param result First element of the multiplication.
   * @param to_mul Second element of the multiplication.
   * @return Entry-wise \f$ result * to\_mul \f$.
   */
  friend Multi_parameter_filtration &operator*=(Multi_parameter_filtration &f, const T &val)
  {
    return f._apply_operation(val, [](T &valF, const T &valR) { _multiply(valF, valR); });
  }

  // division
  /**
   * @brief Returns a filtration value such that an entry at index \f$ i \f$ is equal to \f$ a[i] / b[i] \f$.
   * If @p a and @p b are both not infinite or NaN, they have to have the same number of parameters.
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
   * If `std::numeric_limits<T>::has_quiet_NaN` is false, then the returned filtration value will be @ref nan()
   * if any operation results in NaN, not only if all operations result in NaN.
   *
   * @param a First element of the division.
   * @param b Second element of the division.
   * @return Entry-wise \f$ a / b \f$.
   */
  template <class ValueRange, class = std::enable_if_t<!std::is_arithmetic_v<ValueRange> > >
  friend Multi_parameter_filtration operator/(Multi_parameter_filtration f, const ValueRange &r)
  {
    f /= r;
    return f;
  }

  template <class ValueRange, class = std::enable_if_t<!std::is_arithmetic_v<ValueRange> > >
  friend Multi_parameter_filtration operator/(const ValueRange &r, Multi_parameter_filtration f)
  {
    return f._apply_operation(r, [](T &valF, const T &valR) {
      T tmp = valF;
      valF = valR;
      _divide(valF, tmp);
    });
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ i \f$ is equal to \f$ f[i] / val \f$.
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
   * If `std::numeric_limits<T>::has_quiet_NaN` is false, then the returned filtration value will be @ref nan()
   * if any operation results in NaN, not only if all operations result in NaN.
   *
   * @param f First element of the division.
   * @param val Second element of the division.
   * @return Entry-wise \f$ f / val \f$.
   */
  friend Multi_parameter_filtration operator/(Multi_parameter_filtration f, const T &val)
  {
    f /= val;
    return f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ i \f$ is equal to \f$ val / f[i] \f$.
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
   * If `std::numeric_limits<T>::has_quiet_NaN` is false, then the returned filtration value will be @ref nan()
   * if any operation results in NaN, not only if all operations result in NaN.
   *
   * @param val First element of the division.
   * @param f Second element of the division.
   * @return Entry-wise \f$ val / f \f$.
   */
  friend Multi_parameter_filtration operator/(const T &val, Multi_parameter_filtration f)
  {
    return f._apply_operation(val, [](T &valF, const T &valR) {
      T tmp = valF;
      valF = valR;
      _divide(valF, tmp);
    });
  }

  /**
   * @brief Modifies the first parameters such that an entry at index \f$ i \f$ is equal to
   * \f$ result[i] / to\_div[i] \f$.
   * If @p result and @p to_div are both not infinite or NaN, they have to have the same number of parameters.
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
   * If `std::numeric_limits<T>::has_quiet_NaN` is false, then the returned filtration value will be @ref nan()
   * if any operation results in NaN, not only if all operations result in NaN.
   *
   * @param result First element of the division.
   * @param to_div Second element of the division.
   * @return Entry-wise \f$ result / to\_div \f$.
   */
  template <class ValueRange, class = std::enable_if_t<!std::is_arithmetic_v<ValueRange> > >
  friend Multi_parameter_filtration &operator/=(Multi_parameter_filtration &f, const ValueRange &r)
  {
    return f._apply_operation(r, [](T &valF, const T &valR) { _divide(valF, valR); });
  }

  /**
   * @brief Modifies the first parameters such that an entry at index \f$ i \f$ is equal to
   * \f$ result[i] / to\_div \f$.
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
   * If `std::numeric_limits<T>::has_quiet_NaN` is false, then the returned filtration value will be @ref nan()
   * if any operation results in NaN, not only if all operations result in NaN.
   *
   * @param result First element of the division.
   * @param to_div Second element of the division.
   * @return Entry-wise \f$ result / to\_div \f$.
   */
  friend Multi_parameter_filtration &operator/=(Multi_parameter_filtration &f, const T &val)
  {
    return f._apply_operation(val, [](T &valF, const T &valR) { _divide(valF, valR); });
  }

  // MODIFIERS

  /**
   * @brief Sets the number of generators. If there were less generators before, new empty generators are constructed.
   * If there were more generators before, the exceed of generators is destroyed (any generator with index higher or
   * equal than @p n to be more precise). If @p n is zero, the methods does nothing. A filtration value should never
   * be empty.
   *
   * @warning All empty generators have 0 parameters. This can be problematic for some methods if there are also
   * non empty generators in the container. Make sure to fill them with real generators or to remove them before
   * using those methods.
   *
   * @warning Be sure to call @ref simplify if necessary after setting all the generators. Most methods will have an
   * undefined behaviour if the set of generators is not minimal or sorted.
   *
   * @param n New number of generators.
   */
  void set_num_generators(size_type g)
  {
    if (g == 0) return;
    generators_.resize(g * num_parameters());
    generator_view_.update_extent(0, g);
  }

  /**
   * @brief Adds the given generator to the filtration value such that the sets remains minimal.
   * It is therefore possible that the generator is ignored if it does not generated any new lifetime or that
   * old generators disappear if they are overshadowed by the new one.
   * @pre If all are finite, the new generator has to have the same number of parameters than the others.
   *
   * @param x New generator to add.
   * @return true If and only if the generator is actually added to the set of generators.
   * @return false Otherwise.
   */
  template <class GeneratorRange = std::initializer_list<T> >
  bool add_generator(const GeneratorRange &x)
  {
    GUDHI_CHECK(x.size() == num_parameters(), "Wrong range size. Should correspond to the number of parameters.");

    const int newIndex = -1;

    std::size_t end = num_generators();
    std::vector<int> indices(end);
    std::iota(indices.begin(), indices.end(), 0);

    if (_generator_can_be_added(x.begin(), 0, end, indices)) {
      indices.resize(end);
      indices.push_back(newIndex);
      _build_from(indices, newIndex, x);
      return true;
    }

    return false;
  }

  /**
   * @brief Adds the given generator to the filtration value without any verifications or simplifications at the end
   * of the set.
   *
   * @warning If the resulting set of generators is not minimal or sorted after modification, some methods will have an
   * undefined behaviour. Be sure to call @ref simplify() before using them.
   *
   * @param x
   */
  template <class GeneratorRange = std::initializer_list<T> >
  void add_guaranteed_generator(const GeneratorRange &x)
  {
    generators_.insert(generators_.end(), x.begin(), x.end());
    generator_view_.update_extent(0, num_generators() + 1);
  }

  /**
   * @brief Simplifies the current set of generators such that it becomes minimal. Also orders it in increasing
   * lexicographical order. Only necessary if generators were added "by hand" without verification either trough the
   * constructor or with @ref add_guaranteed_generator "", etc.
   */
  void simplify()
  {
    std::size_t end = 0;
    std::vector<int> indices(num_generators());
    std::iota(indices.begin(), indices.end(), 0);

    for (std::size_t curr = 0; curr < num_generators(); ++curr) {
      if (_generator_can_be_added(generators_.begin() + generator_view_.get_index(indices[curr], 0), 0, end, indices)) {
        std::swap(indices[end], indices[curr]);
        ++end;
      }
    }

    indices.resize(end);
    _build_from(indices);
  }

  /**
   * @brief Removes all empty generators from the filtration value. If @p include_infinities is true, it also
   * removes the generators at infinity or minus infinity. As empty generators are not possible (except if number of
   * parameters is 0), the method does nothing except sorting the set of generators if @p include_infinities is false.
   * Exists mostly for interface purposes.
   * If the set of generators is empty after removals, it is set to minus infinity if `co` is false or to infinity
   * if `co` is true.
   *
   * @warning If the resulting set of generators is not minimal after the removals/sorting, some methods will have an
   * undefined behaviour. Be sure to call @ref simplify() before using them.
   *
   * @param include_infinities If true, removes also infinity values.
   */
  void remove_empty_generators(bool include_infinities = false)
  {
    std::vector<int> indices;
    indices.reserve(num_generators());
    for (unsigned int i = 0; i < num_generators(); ++i) {
      if (!include_infinities || _is_finite(i)) indices.push_back(i);
    }
    _build_from(indices);  // sorts

    if (generators_.empty()) {
      generators_.resize(num_parameters(), _get_default_value());
      generator_view_.update_extent(0, 1);
    }
  }

  /**
   * @brief Sets the filtration value to the least common upper bound between the current value and the given value.
   *
   * More formally, it pushes the current filtration value to the cone \f$ \{ y \in \mathbb R^n : y \ge x \} \f$
   * originating in the given filtration value \f$ x \f$. The resulting value corresponds to the intersection of both
   * cones: \f$ \mathrm{this} = \min \{ y \in \mathbb R^n : y \ge this \} \cap \{ y \in \mathbb R^n : y \ge x \} \f$.
   *
   * @param x The target filtration value towards which to push.
   * @return True if and only if the value of this actually changed.
   */
  template <class GeneratorRange = std::initializer_list<value_type> >
  bool push_to_least_common_upper_bound(const GeneratorRange &x, bool exclude_infinite_values = false)
  {
    GUDHI_CHECK(x.size() == num_parameters(), "Wrong range size. Should correspond to the number of parameters.");

    bool xIsInf = true, xIsMinusInf = true, xIsNaN = true;
    bool thisIsInf = true, thisIsMinusInf = true, thisIsNaN = true;

    // if one is not finite, we can avoid the heavy simplification process
    _get_infinity_statuses(generator_view_, x, thisIsInf, thisIsMinusInf, thisIsNaN, xIsInf, xIsMinusInf, xIsNaN);

    if (thisIsInf || thisIsNaN || xIsNaN || xIsMinusInf || (xIsInf && exclude_infinite_values)) return false;

    if (xIsInf || thisIsMinusInf) {
      generators_ = {x};
      generator_view_.update_extent(0, 1);
      return true;
    }

    bool modified = false;
    for (size_type g = 0; g < num_generators(); ++g) {
      auto it = x.begin();
      for (unsigned int p = 0; p < num_parameters(); ++p) {
        T valX = *it;
        ++it;
        if (exclude_infinite_values && (valX == T_inf || valX == -T_inf)) continue;
        T &val = generator_view_(g, p);
        modified |= val < valX;
        val = valX > val ? valX : val;
      }
    }

    if (modified && num_generators() > 1) simplify();

    return modified;
  }

  /**
   * @brief Sets the filtration value to the greatest common lower bound between the current value and the given value.
   *
   * More formally, it pulls the current filtration value to the cone \f$ \{ y \in \mathbb R^n : y \le x \} \f$
   * originating in the given filtration value \f$ x \f$. The resulting value corresponds to the intersection of both
   * cones: \f$ \mathrm{this} = \min \{ y \in \mathbb R^n : y \le this \} \cap \{ y \in \mathbb R^n : y \le x \} \f$.
   *
   * @param x The target filtration value towards which to pull.
   * @return True if and only if the value of this actually changed.
   */
  template <class GeneratorRange = std::initializer_list<value_type> >
  bool pull_to_greatest_common_lower_bound(const GeneratorRange &x, bool exclude_infinite_values = false)
  {
    GUDHI_CHECK(x.size() == num_parameters(), "Wrong range size. Should correspond to the number of parameters.");

    bool xIsInf = true, xIsMinusInf = true, xIsNaN = true;
    bool thisIsInf = true, thisIsMinusInf = true, thisIsNaN = true;

    // if one is not finite, we can avoid the heavy simplification process
    _get_infinity_statuses(generator_view_, x, thisIsInf, thisIsMinusInf, thisIsNaN, xIsInf, xIsMinusInf, xIsNaN);

    if (xIsInf || thisIsNaN || xIsNaN || thisIsMinusInf || (xIsMinusInf && exclude_infinite_values)) return false;

    if (thisIsInf || xIsMinusInf) {
      generators_ = {x};
      generator_view_.update_extent(0, 1);
      return true;
    }

    bool modified = false;
    for (size_type g = 0; g < num_generators(); ++g) {
      auto it = x.begin();
      for (unsigned int p = 0; p < num_parameters(); ++p) {
        T valX = *it;
        ++it;
        if (exclude_infinite_values && (valX == T_inf || valX == -T_inf)) continue;
        T &val = generator_view_(g, p);
        modified |= val > valX;
        val = valX < val ? valX : val;
      }
    }

    if (modified && num_generators() > 1) simplify();

    return modified;
  }

  /**
   * @brief Projects the filtration value into the given grid. If @p coordinate is false, the entries are set to
   * the nearest upper bound value with the same parameter in the grid. Otherwise, the entries are set to the indices
   * of those nearest upper bound values.
   * The grid has to be represented as a vector of ordered ranges of values convertible into `T`. An index
   * \f$ i \f$ of the vector corresponds to the same parameter as the index \f$ i \f$ in the filtration value.
   * The ranges correspond to the possible values of the parameters, ordered by increasing value, forming therefore
   * all together a 2D grid.
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
    GUDHI_CHECK(grid.size() >= num_parameters(),
                "The grid should not be smaller than the number of parameters in the filtration value.");

    for (size_type g = 0; g < num_generators(); ++g) {
      for (size_type p = 0; p < num_parameters(); ++p) {
        const auto &filtration = grid[p];
        auto d = std::distance(filtration.begin(),
                               std::lower_bound(filtration.begin(),
                                                filtration.end(),
                                                static_cast<typename OneDimArray::value_type>(generator_view_(g, p))));
        generator_view_(g, p) = coordinate ? static_cast<T>(d) : static_cast<T>(filtration[d]);
      }
    }

    if (!coordinate) simplify();
  }

  // FONCTIONNALITIES

  /**
   * @brief Returns a generator with the minimal values of all parameters in any generator of the given filtration
   * value. That is, the greatest lower bound of all generators.
   */
  friend Multi_parameter_filtration factorize_below(const Multi_parameter_filtration &f)
  {
    if (f.num_generators() == 0) return f;

    std::vector<T> result(f.num_parameters(), T_inf);
    for (size_type g = 0; g < f.num_generators(); ++g) {
      for (size_type p = 0; p < f.num_parameters(); ++p) {
        result[p] = std::min(result[p], f(g, p));
      }
    }
    return Multi_parameter_filtration(std::move(result), f.num_parameters());
  }

  /**
   * @brief Returns a generator with the maximal values of all parameters in any generator of the given filtration
   * value. That is, the least upper bound of all generators.
   */
  friend Multi_parameter_filtration factorize_above(const Multi_parameter_filtration &f)
  {
    if (f.num_generators() == 0) return f;

    std::vector<T> result(f.num_parameters(), -T_inf);
    for (size_type g = 0; g < f.num_generators(); ++g) {
      for (size_type p = 0; p < f.num_parameters(); ++p) {
        result[p] = std::max(result[p], f(g, p));
      }
    }
    return Multi_parameter_filtration(std::move(result), f.num_parameters());
  }

  /**
   * @brief Computes the scalar product of the given filtration value with the given vector.
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

    if constexpr (co) {
      U projection = std::numeric_limits<U>::lowest();
      for (size_type g = 0; g < f.num_generators(); ++g) {
        projection = std::max(projection, project_generator(g));
      }
      return projection;
    } else {
      U projection = std::numeric_limits<U>::max();
      for (size_type g = 0; g < f.num_generators(); ++g) {
        projection = std::min(projection, project_generator(g));
      }
      return projection;
    }
  }

  /**
   * @brief Computes the euclidean distance from the first parameter to the second parameter.
   *
   * @param f Start filtration value.
   * @param other End filtration value.
   * @return Euclidean distance between @p f and @p other.
   */
  template <typename U = T>
  friend U compute_euclidean_distance_to(const Multi_parameter_filtration &f, const Multi_parameter_filtration &other)
  {
    GUDHI_CHECK(f.num_parameters() == other.num_parameters(),
                "We cannot compute the distance between two points of different dimensions.");

    U res = std::numeric_limits<U>::max();
    for (size_type g1 = 0; g1 < f.num_generators(); ++g1) {
      for (size_type g2 = 0; g2 < other.num_generators(); ++g2) {
        // Euclidean distance as a Frobenius norm with matrix 1 x p and values 'f(g1, p) - other(g2, p)'
        res = std::min(res, _compute_frobenius_norm(f.num_parameters(), [&](size_type p) -> U {
                         return f(g1, p) - other(g2, p);
                       }));
      }
    }
    return res;
  }

  /**
   * @brief Computes the norm of the given filtration value.
   *
   * @param f Filtration value.
   * @return The norm of @p f.
   */
  template <typename U = T>
  friend U compute_norm(const Multi_parameter_filtration &f)
  {
    // Frobenius norm with matrix g x p based on Euclidean norm
    return _compute_frobenius_norm(f.num_generators(), [&](size_type g) -> U {
      // Euclidean norm as Frobenius norm with matrix of rank 1 x p and values 'f(g, p)'
      return _compute_frobenius_norm(f.num_parameters(), [&](size_type p) -> U { return f(g, p); });
    });
  }

  /**
   * @brief Computes the coordinates in the given grid, corresponding to the nearest upper bounds of the entries
   * in the given filtration value.
   * The grid has to be represented as a vector of vectors of ordered values convertible into `out_type`. An index
   * \f$ i \f$ of the vector corresponds to the same parameter as the index \f$ i \f$ in the filtration value.
   * The inner vectors correspond to the possible values of the parameters, ordered by increasing value,
   * forming therefore all together a 2D grid.
   *
   * @tparam out_type Signed arithmetic type. Default value: std::int32_t.
   * @tparam U Type which is convertible into `out_type`.
   * @param f Filtration value to project.
   * @param grid Vector of vectors to project into.
   * @return Filtration value \f$ out \f$ whose entry correspond to the indices of the projected values. That is,
   * the projection of \f$ f[i] \f$ is \f$ grid[i][out[i]] \f$.
   */
  template <typename OutValue = std::int32_t, typename U = T>
  friend Multi_parameter_filtration<OutValue, co> compute_coordinates_in_grid(Multi_parameter_filtration f,
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
   * value. That is, if \f$ out \f$ is the result, \f$ out[i] = grid[i][f[i]] \f$. Assumes therefore, that the
   * values stored in the filtration value corresponds to indices existing in the given grid.
   *
   * @tparam U Signed arithmetic type.
   * @param f Filtration value storing coordinates compatible with `grid`.
   * @param grid Vector of vector.
   * @return Filtration value \f$ out \f$ whose entry correspond to \f$ out[i] = grid[i][f[i]] \f$.
   */
  template <typename U>
  friend Multi_parameter_filtration<U, co> evaluate_coordinates_in_grid(const Multi_parameter_filtration &f,
                                                                        const std::vector<std::vector<U> > &grid)
  {
    GUDHI_CHECK(grid.size() > f.num_parameters(),
                "The size of the grid should correspond to the number of parameters in the filtration value.");

    U grid_inf = Multi_parameter_filtration<U, co>::T_inf;
    std::vector<U> outVec(f.num_entries());

    for (size_type g = 0; g < f.num_generators(); ++g) {
      for (size_type p = 0; p < f.num_parameters(); ++p) {
        const std::vector<U> &filtration = grid[p];
        const T &c = f(g, p);
        outVec[f.generator_view_.get_index(g, p)] = (c == f.T_inf ? grid_inf : filtration[c]);
      }
    }

    Multi_parameter_filtration<U, co> out(std::move(outVec), f.num_parameters());
    out.simplify();
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

    stream << "( k = " << num_gen << " ) [ ";
    for (size_type g = 0; g < num_gen; ++g) {
      stream << "[";
      for (size_type p = 0; p < num_param - 1; ++p) {
        stream << f(g, p);
        if (p < num_param - 1) stream << ", ";
      }
      stream << "]";
      if (g < num_gen - 1) stream << "; ";
    }
    stream << " ]";

    return stream;
  }

  friend bool unify_lifetimes(Multi_parameter_filtration &f1, const Multi_parameter_filtration &f2)
  {
    return f1.pull_to_greatest_common_lower_bound(f2);
  }

  friend bool intersect_lifetimes(Multi_parameter_filtration &f1, const Multi_parameter_filtration &f2)
  {
    return f1.push_to_least_common_upper_bound(f2);
  }

  friend char *serialize_trivial(const Multi_parameter_filtration &value, char *start)
  {
    const size_type length = value.generators_.size();
    const size_type num_param = value.num_parameters();
    const std::size_t arg_size = sizeof(T) * length;
    const std::size_t type_size = sizeof(size_type);
    memcpy(start, &num_param, type_size);
    memcpy(start, &length, type_size);
    memcpy(start + type_size, value.generators_.data(), arg_size);
    return start + arg_size + type_size;
  }

  friend const char *deserialize_trivial(Multi_parameter_filtration &value, const char *start)
  {
    const std::size_t type_size = sizeof(size_type);
    size_type length;
    size_type num_param;
    memcpy(&num_param, start, type_size);
    memcpy(&length, start, type_size);
    std::size_t arg_size = sizeof(T) * length;
    value.generators_.resize(length);
    memcpy(value.generators_.data(), start + type_size, arg_size);
    value.generator_view_ = Simple_mdspan<T>(
        value.generators_.data(), num_param == 0 ? 0 : value.generators_.size() / num_param, num_param);
    return start + arg_size + type_size;
  }

  friend std::size_t get_serialization_size_of(const Multi_parameter_filtration &value)
  {
    return sizeof(size_type) * 2 + sizeof(T) * value.num_entries();
  }

 private:
  Underlying_container generators_;
  Gudhi::Simple_mdspan<T> generator_view_;  // has to be created after generators_
  /**
   * @brief Infinity value of an entry of the filtration value.
   */
  constexpr static const T T_inf =
      std::numeric_limits<T>::has_infinity ? std::numeric_limits<T>::infinity() : std::numeric_limits<T>::max();

  constexpr static T _get_default_value() { return co ? T_inf : -T_inf; }

  constexpr static bool _is_nan(T val)
  {
    if constexpr (std::is_integral_v<T>) {
      // to avoid Windows issue which don't know how to cast integers for cmath methods
      return false;
    } else {
      return std::isnan(val);
    }
  }

  /**
   * @brief Verifies if @p b is strictly contained in the positive cone originating in `a`.
   */
  static bool _strictly_contains(const Gudhi::Simple_mdspan<T> &a,
                                 size_type g_a,
                                 const Gudhi::Simple_mdspan<T> &b,
                                 size_type g_b)
  {
    bool isSame = true;
    for (auto i = 0u; i < a.extent(1); ++i) {
      T a_i, b_i;
      if constexpr (co) {
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
  static bool _contains(const Gudhi::Simple_mdspan<T> &a,
                        size_type g_a,
                        const Gudhi::Simple_mdspan<T> &b,
                        size_type g_b)
  {
    for (std::size_t i = 0u; i < a.extent(1); ++i) {
      T a_i, b_i;
      if constexpr (co) {
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

  static bool _first_strictly_dominates(const Gudhi::Simple_mdspan<T> &a,
                                        size_type g_a,
                                        const Gudhi::Simple_mdspan<T> &b,
                                        size_type g_b)
  {
    if constexpr (co) {
      return a(g_a, 0) < b(g_b, 0);
    } else {
      return a(g_a, 0) > b(g_b, 0);
    }
  }

  static bool _first_dominates(const Gudhi::Simple_mdspan<T> &a,
                               size_type g_a,
                               const Gudhi::Simple_mdspan<T> &b,
                               size_type g_b)
  {
    if constexpr (co) {
      return a(g_a, 0) <= b(g_b, 0);
    } else {
      return a(g_a, 0) >= b(g_b, 0);
    }
  }

  template <class ValueRange, class F, class = std::enable_if_t<!std::is_arithmetic_v<ValueRange> > >
  Multi_parameter_filtration &_apply_operation(const ValueRange &range, F &&operate)
  {
    auto &view = generator_view_;
    for (unsigned int g = 0; g < num_generators(); ++g) {
      auto it = range.begin();
      for (unsigned int p = 0; p < num_parameters() && it != range.end(); ++p) {
        operate(view(g, p), *it);
        ++it;
      }
    }
    return *this;
  }

  template <class F>
  Multi_parameter_filtration &_apply_operation(const T &val, F &&operate)
  {
    auto &gens = generators_;
    for (unsigned int i = 0; i < gens.size(); ++i) {
      operate(gens[i], val);
    }
    return *this;
  }

  constexpr static bool _subtract(T &v1, T v2) { return add_(v1, -v2); }

  constexpr static bool _add(T &v1, T v2)
  {
    if (_is_nan(v1) || _is_nan(v2) || (v1 == T_inf && v2 == -T_inf) || (v1 == -T_inf && v2 == T_inf)) {
      v1 = std::numeric_limits<T>::quiet_NaN();
      return false;
    }
    if (v1 == T_inf || v1 == -T_inf) {
      return true;
    }
    if (v2 == T_inf || v2 == -T_inf) {
      v1 = v2;
      return true;
    }

    v1 += v2;
    return true;
  }

  constexpr static bool _multiply(T &v1, T v2)
  {
    bool v1_is_infinite = v1 == T_inf || v1 == -T_inf;
    bool v2_is_infinite = v2 == T_inf || v2 == -T_inf;

    if (_is_nan(v1) || _is_nan(v2) || (v1_is_infinite && v2 == 0) || (v1 == 0 && v2_is_infinite)) {
      v1 = std::numeric_limits<T>::quiet_NaN();
      return false;
    }

    if ((v1 == T_inf && v2 > 0) || (v1 == -T_inf && v2 < 0) || (v1 < 0 && v2 == -T_inf) || (v1 > 0 && v2 == T_inf)) {
      v1 = T_inf;
      return true;
    }

    if ((v1 == T_inf && v2 < 0) || (v1 == -T_inf && v2 > 0) || (v1 > 0 && v2 == -T_inf) || (v1 < 0 && v2 == T_inf)) {
      v1 = -T_inf;
      return true;
    }

    v1 *= v2;
    return true;
  }

  constexpr static bool _divide(T &v1, T v2)
  {
    bool v1_is_infinite = v1 == T_inf || v1 == -T_inf;
    bool v2_is_infinite = v2 == T_inf || v2 == -T_inf;

    if (_is_nan(v1) || _is_nan(v2) || v2 == 0 || (v1_is_infinite && v2_is_infinite)) {
      v1 = std::numeric_limits<T>::quiet_NaN();
      return false;
    }

    if (v1 == 0 || (v1_is_infinite && v2 > 0)) return true;

    if (v1_is_infinite && v2 < 0) {
      v1 = -v1;
      return true;
    }

    if (v2_is_infinite) {
      v1 = 0;
      return true;
    }

    v1 /= v2;
    return true;
  }

  template <class GeneratorRange>
  static void _get_infinity_statuses(const Gudhi::Simple_mdspan<T> &a,
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
      if (a(0, i) != -T_inf) aIsMinusInf = false;
      if (!_is_nan(a(0, i))) aIsNaN = false;
      if (*itB != T_inf) bIsInf = false;
      if (*itB != -T_inf) bIsMinusInf = false;
      if (!_is_nan(*itB)) bIsNaN = false;
      if (!aIsInf && !aIsMinusInf && !aIsNaN && !bIsInf && !bIsMinusInf && !bIsNaN) return;
      ++itB;
    }
  }

  enum class Rel { EQUAL, DOMINATES, IS_DOMINATED, NONE };

  template <class Iterator>
  static Rel _get_domination_relation(const Gudhi::Simple_mdspan<T> &a, size_type g_a, Iterator itB)
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

    if constexpr (co) {
      if (allSmaller) return Rel::DOMINATES;
      return Rel::IS_DOMINATED;
    } else {
      if (allGreater) return Rel::DOMINATES;
      return Rel::IS_DOMINATED;
    }
  }

  // assumes between 'curr' and 'end' everything is simplified:
  // no nan values and if there is an inf/-inf, then 'end - curr == 1'
  // modifies generators_ only if true is returned.
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

  template <class GeneratorRange = std::initializer_list<T> >
  void _build_from(std::vector<int> &indices, const int newIndex, const GeneratorRange &x)
  {
    auto comp = [&](int g1, int g2) -> bool {
      if (g1 == g2) {
        return false;
      }

      if (g1 == newIndex) {
        auto it = x.begin();
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
        auto it = x.begin();
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
        new_container.insert(new_container.end(), x.begin(), x.end());
      } else {
        T *ptr = &generator_view_(i, 0);
        for (size_type p = 0; p < num_parameters(); ++p) {
          new_container.push_back(*ptr);
          ++ptr;
        }
      }
    }
    generators_.swap(new_container);
    generator_view_ = Simple_mdspan<T>(generators_.data(), generators_.size() / num_parameters(), num_parameters());
  }

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
    generator_view_ = Simple_mdspan<T>(generators_.data(), generators_.size() / num_parameters(), num_parameters());
  }

  bool _is_finite(size_type g)
  {
    bool isInf = true, isMinusInf = true, isNan = true;
    for (size_type p = 0; p < num_parameters(); ++p) {
      T v = generator_view_(g, p);
      if (v != T_inf) isInf = false;
      if (v != -T_inf) isMinusInf = false;
      if (!_is_nan(v)) isNan = false;
      if (!isInf && !isMinusInf && !isNan) return true;
    }
    return false;
  }

  template <class F, typename U = T, class... Args>
  static U _compute_frobenius_norm(size_type number_of_elements, F &&norm)
  {
    if (number_of_elements == 1) return norm(0);

    U out = 0;
    for (size_type p = 0; p < number_of_elements; ++p) {
      T v = norm(p);
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

template <typename T, bool co>
class numeric_limits<Gudhi::multi_filtration::Multi_parameter_filtration<T, co> >
{
 public:
  static constexpr bool has_infinity = true;
  static constexpr bool has_quiet_NaN = std::numeric_limits<T>::has_quiet_NaN;

  static constexpr Gudhi::multi_filtration::Multi_parameter_filtration<T, co> infinity() noexcept(false)
  {
    throw std::logic_error(
        "The infinite value cannot be represented with no finite numbers of parameters."
        "Use `infinity(number_of_parameters)` instead");
  };

  static constexpr Gudhi::multi_filtration::Multi_parameter_filtration<T, co> infinity(std::size_t p) noexcept
  {
    return Gudhi::multi_filtration::Multi_parameter_filtration<T, co>::inf(p);
  };

  // non-standard
  static constexpr Gudhi::multi_filtration::Multi_parameter_filtration<T, co> minus_infinity() noexcept(false)
  {
    throw std::logic_error(
        "The infinite value cannot be represented with no finite numbers of parameters."
        "Use `minus_infinity(number_of_parameters)` instead");
  };

  static constexpr Gudhi::multi_filtration::Multi_parameter_filtration<T, co> minus_infinity(std::size_t p) noexcept
  {
    return Gudhi::multi_filtration::Multi_parameter_filtration<T, co>::minus_inf(p);
  };

  static constexpr Gudhi::multi_filtration::Multi_parameter_filtration<T, co> max() noexcept(false)
  {
    throw std::logic_error(
        "The max value cannot be represented with no finite numbers of parameters."
        "Use `max(number_of_parameters)` instead");
  };

  static constexpr Gudhi::multi_filtration::Multi_parameter_filtration<T, co> max(std::size_t p) noexcept
  {
    return Gudhi::multi_filtration::Multi_parameter_filtration<T, co>(p, std::numeric_limits<T>::max());
  };

  static constexpr Gudhi::multi_filtration::Multi_parameter_filtration<T, co> quiet_NaN() noexcept(false)
  {
    throw std::logic_error(
        "The NaN value cannot be represented with no finite numbers of parameters."
        "Use `quiet_NaN(number_of_parameters)` instead");
  };

  static constexpr Gudhi::multi_filtration::Multi_parameter_filtration<T, co> quiet_NaN(std::size_t p) noexcept(false)
  {
    if constexpr (std::numeric_limits<T>::has_quiet_NaN) {
      return Gudhi::multi_filtration::Multi_parameter_filtration<T, co>::nan(p);
    } else {
      throw std::logic_error("Does not have a NaN value.");
    }
  };
};

}  // namespace std

#endif  // MF_MULTI_PARAMETER_FILTRATION_H_
