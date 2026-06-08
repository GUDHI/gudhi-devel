/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2025 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @private
 * @file multi_filtration_utils.h
 * @author Hannah Schreiber
 */

#ifndef MF_UTILS_H_
#define MF_UTILS_H_

#include <cstddef>
#include <istream>
#include <stdexcept>
#include <type_traits>
#include <limits>
#include <cmath>

namespace Gudhi {

namespace multi_filtration {

/**
 * @private
 */
template <typename T>
class RangeTraits
{
 private:
  static auto check_begin(...) -> std::false_type;
  template <typename U>
  static auto check_begin(U x) -> decltype(x.begin(), std::true_type{});

  static auto check_dynamic_filtration(...) -> std::false_type;
  template <typename U>
  static auto check_dynamic_filtration(U x) -> decltype(x.operator[](std::size_t{}), std::true_type{});

  static auto check_filtration(...) -> std::false_type;
  template <typename U>
  static auto check_filtration(U x) -> decltype(x.ensures_1_criticality(), std::true_type{});

 public:
  static constexpr bool has_begin = decltype(check_begin(std::declval<T>()))::value;
  static constexpr bool is_multi_filtration = decltype(check_filtration(std::declval<T>()))::value;
  static constexpr bool is_dynamic_multi_filtration =
      is_multi_filtration && decltype(check_dynamic_filtration(std::declval<T>()))::value;
};

/**
 * @private
 */
template <typename T>
constexpr bool _is_nan(T val)
{
  if constexpr (std::is_integral_v<T>) {
    // to avoid Windows issue which don't know how to cast integers for cmath methods
    return false;
  } else {
    return std::isnan(val);
  }
};

/**
 * @private
 * @brief Infinity value of an entry of the filtration value.
 */
template <typename T>
constexpr const T MF_T_inf =
    std::numeric_limits<T>::has_infinity ? std::numeric_limits<T>::infinity() : std::numeric_limits<T>::max();

/**
 * @private
 * @brief Minus infinity value of an entry of the filtration value.
 */
template <typename T>
constexpr const T MF_T_m_inf =
    std::numeric_limits<T>::has_infinity ? -std::numeric_limits<T>::infinity() : std::numeric_limits<T>::lowest();

/**
 * @private
 * @brief Adds v1 and v2, stores the result in v1 and returns true if and only if v1 was modified.
 */
template <typename T>
constexpr bool _add(T &v1, T v2)
{
  if (_is_nan(v1) || _is_nan(v2) || (v1 == MF_T_inf<T> && v2 == MF_T_m_inf<T>) ||
      (v1 == MF_T_m_inf<T> && v2 == MF_T_inf<T>)) {
    v1 = std::numeric_limits<T>::quiet_NaN();
    return false;
  }
  if (v1 == MF_T_inf<T> || v1 == MF_T_m_inf<T>) {
    return true;
  }
  if (v2 == MF_T_inf<T> || v2 == MF_T_m_inf<T>) {
    v1 = v2;
    return true;
  }

  v1 += v2;
  return true;
};

/**
 * @private
 * @brief Subtracts v1 and v2, stores the result in v1 and returns true if and only if v1 was modified.
 */
template <typename T>
constexpr bool _subtract(T &v1, T v2)
{
  return _add(v1, v2 == MF_T_inf<T> ? MF_T_m_inf<T> : (v2 == MF_T_m_inf<T> ? MF_T_inf<T> : -v2));
};

/**
 * @private
 * @brief Multiplies v1 and v2, stores the result in v1 and returns true if and only if v1 was modified.
 */
template <typename T>
constexpr bool _multiply(T &v1, T v2)
{
  bool v1_is_infinite = v1 == MF_T_inf<T> || v1 == MF_T_m_inf<T>;
  bool v2_is_infinite = v2 == MF_T_inf<T> || v2 == MF_T_m_inf<T>;

  if (_is_nan(v1) || _is_nan(v2) || (v1_is_infinite && v2 == 0) || (v1 == 0 && v2_is_infinite)) {
    v1 = std::numeric_limits<T>::quiet_NaN();
    return false;
  }

  if ((v1 == MF_T_inf<T> && v2 > 0) || (v1 == MF_T_m_inf<T> && v2 < 0) || (v1 < 0 && v2 == MF_T_m_inf<T>) ||
      (v1 > 0 && v2 == MF_T_inf<T>)) {
    v1 = MF_T_inf<T>;
    return true;
  }

  if ((v1 == MF_T_inf<T> && v2 < 0) || (v1 == MF_T_m_inf<T> && v2 > 0) || (v1 > 0 && v2 == MF_T_m_inf<T>) ||
      (v1 < 0 && v2 == MF_T_inf<T>)) {
    v1 = MF_T_m_inf<T>;
    return true;
  }

  v1 *= v2;
  return true;
};

/**
 * @private
 * @brief Divides v1 and v2, stores the result in v1 and returns true if and only if v1 was modified.
 */
template <typename T>
constexpr bool _divide(T &v1, T v2)
{
  bool v1_is_infinite = v1 == MF_T_inf<T> || v1 == MF_T_m_inf<T>;
  bool v2_is_infinite = v2 == MF_T_inf<T> || v2 == MF_T_m_inf<T>;

  if (_is_nan(v1) || _is_nan(v2) || v2 == 0 || (v1_is_infinite && v2_is_infinite)) {
    v1 = std::numeric_limits<T>::quiet_NaN();
    return false;
  }

  if (v1 == 0 || (v1_is_infinite && v2 > 0)) return true;

  if (v1_is_infinite && v2 < 0) {
    v1 = v1 == MF_T_inf<T> ? MF_T_m_inf<T> : (v1 == MF_T_m_inf<T> ? MF_T_inf<T> : -v1);
    return true;
  }

  if (v2_is_infinite) {
    v1 = 0;
    return true;
  }

  v1 /= v2;
  return true;
};

template <class T>
T _get_value(std::istream &stream)
{
  if constexpr (std::numeric_limits<T>::has_infinity) {
    auto pos = stream.tellg();
    char first;
    stream >> first;
    if (first == 'i') {
      stream >> first;  // n
      stream >> first;  // f
      return std::numeric_limits<T>::infinity();
    }
    if (first == '-') {
      stream >> first;  // i
      if (first == 'i') {
        stream >> first;  // n
        stream >> first;  // f
        return -std::numeric_limits<T>::infinity();
      }  // else could be a negative number
    }
    if (first == 'n') {
      if constexpr (std::numeric_limits<T>::has_quiet_NaN) {
        stream >> first;  // a
        stream >> first;  // n
        return std::numeric_limits<T>::quiet_NaN();
      } else {
        throw std::invalid_argument("Wrong input stream format for value, no nan values allowed.");
      }
    }
    stream.seekg(pos, std::ios_base::beg);
  }

  T val;
  stream >> val;
  if (stream.fail()) throw std::invalid_argument("Wrong input stream format for value.");
  return val;
};

}  // namespace multi_filtration

}  // namespace Gudhi

#endif  // MF_UTILS_H_
