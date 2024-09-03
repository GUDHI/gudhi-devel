/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       David Loiseaux
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - 2024/08 Hannah Schreiber: Generalization to all signed arithmetic types for T + doc
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef ONE_CRITICAL_FILTRATIONS_H_
#define ONE_CRITICAL_FILTRATIONS_H_

#include <algorithm>  //std::lower_bound
#include <cmath>      //std::isnan
#include <cstddef>    //std::size_t
#include <cstdint>    //std::int32_t
#include <ostream>    //std::ostream
#include <limits>     //std::numerical_limits
#include <stdexcept>  //std::logic_error
#include <vector>

#include <gudhi/Debug_utils.h>

namespace Gudhi::multi_filtration {

/**
 * @class One_critical_filtration one_critical_filtration.h gudhi/one_critical_filtration.h
 * @ingroup multi_filtration
 *
 * @brief Class encoding the apparition time, i.e., filtration value of an object 
 * (e.g., simplex, cell, abstract algebraic generator) in the setting of 1-critical multiparameter filtrations. 
 * The class can be used as a vector whose indices correspond to one parameter each.
 * It also follows numpy-like broadcast semantic.
 *
 * @details Inherits of `std::vector<T>`. Overloads `std::numeric_limits` such that:
 * - `std::numeric_limits<One_critical_filtration<T> >::has_infinity` returns `true`,
 * - `std::numeric_limits<One_critical_filtration<T> >::infinity()` returns @ref One_critical_filtration<T>::inf() "",
 * - `std::numeric_limits<One_critical_filtration<T> >::minus_infinity()` returns
 *   @ref One_critical_filtration<T>::minus_inf() "",
 * - `std::numeric_limits<One_critical_filtration<T> >::max()` throws,
 * - `std::numeric_limits<One_critical_filtration<T> >::max(n)` returns a @ref One_critical_filtration<T> with `n`
 *   parameters evaluated at value `std::numeric_limits<T>::max()`,
 * - `std::numeric_limits<One_critical_filtration<T> >::quiet_NaN()` returns @ref One_critical_filtration<T>::nan() "".
 *
 * One critical simplicial filtrations are filtrations such that the lifetime of each object is a positive cone, e.g.
 *  - \f$ \{ x \in  \mathbb R^2 : x>=(1,2)\} \f$ is valid, while
 *  - \f$ \{ x \in  \mathbb R^2 : x>=(1,2)\} \cap \{x \in \mathbb R^2 : x>=(2,1)\} \f$ is not.
 *
 * If the lifetime corresponds to a union of such positive cones, the filtration is called a multi-critical filtration.
 * For those cases, use @ref Multi_critical_filtration instead.
 *
 * @tparam T Arithmetic type of an entry for one parameter of the filtration value. Has to be **signed** and
 * to implement `std::isnan(T)`, `std::numeric_limits<T>::has_quiet_NaN`, `std::numeric_limits<T>::quiet_NaN()`,
 * `std::numeric_limits<T>::has_infinity`, `std::numeric_limits<T>::infinity()` and `std::numeric_limits<T>::max()`.
 * If `std::numeric_limits<T>::has_infinity` returns `false`, a call to `std::numeric_limits<T>::infinity()`
 * can simply throw. Examples are the native types `double`, `float` and `int`.
 */
template <typename T>
class One_critical_filtration : public std::vector<T> {
 private:
  using Base = std::vector<T>;

 public:
  /**
   * @brief Type of the origin of a "lifetime cone", i.e., of a one-critical filtration value.
   * Common with @ref Multi_critical_filtration "". In the 1-critical case, simply the class it-self.
   *
   * @tparam U Type of a value for one parameter within the filtration value. Default value: `T`.
   */
  template <typename U = T>
  using Generator = One_critical_filtration<U>;

  // CONSTRUCTORS

  /**
   * @brief Default constructor. Constructs a value at minus infinity.
   */
  One_critical_filtration() : Base{-T_inf} {};
  /**
   * @brief Constructs a vector of the size of the given number of parameters with -inf as value for each entry.
   *
   * @warning The vector `{-inf, -inf, ...}` with \f$ n > 1 \f$ entries is not considered as "minus infinity" (the
   * method @ref is_minus_inf() will not return true). The `-inf` are just meant as placeholders, at least one entry
   * should be modified by the user. Otherwise, either use the static method @ref minus_inf() or set @p n to 1 instead.
   *
   * @param n Number of parameters.
   */
  One_critical_filtration(int n) : Base(n, -T_inf) {};
  /**
   * @brief Constructs a vector of the size of the given number of parameters and the given value for each entry.
   *
   * @warning If @p value is `inf`, `-inf`, or `NaN`, the vector `{value, value, ...}` with \f$ n > 1 \f$ entries
   * is not wrong but will not be considered as respectively "infinity", "minus infinity" or "NaN" (the corresponding
   * methods @ref is_inf(), @ref is_minus_inf() and @ref is_nan() will return false). For this purpose, please use
   * the static methods @ref inf(), @ref minus_inf() and @ref nan() instead.
   *
   * @param n Number of parameters.
   * @param value Value which will be used for each entry.
   */
  One_critical_filtration(int n, T value) : Base(n, value) {};
  /**
   * @brief Construct a vector from the given initializer list.
   *
   * @param init Initializer list with values for each parameter.
   */
  One_critical_filtration(std::initializer_list<T> init) : Base(init) {};
  /**
   * @brief Construct a vector from the given vector.
   *
   * @param v Vector with values for each parameter.
   */
  One_critical_filtration(const std::vector<T> &v) : Base(v) {};
  /**
   * @brief Construct a vector from the given vector by moving it to the new vector.
   *
   * @param v Vector with values for each parameter.
   */
  One_critical_filtration(std::vector<T> &&v) : Base(std::move(v)) {};
  /**
   * @brief Construct a vector from the range given by the begin and end iterators.
   *
   * @param it_begin Start of the range.
   * @param it_end End of the range.
   */
  One_critical_filtration(typename std::vector<T>::iterator it_begin, typename std::vector<T>::iterator it_end)
      : Base(it_begin, it_end) {};
  /**
   * @brief Construct a vector from the range given by the begin and end const iterators.
   *
   * @param it_begin Start of the range.
   * @param it_end End of the range.
   */
  One_critical_filtration(typename std::vector<T>::const_iterator it_begin,
                          typename std::vector<T>::const_iterator it_end)
      : Base(it_begin, it_end) {};

  /**
   * @brief Assign operator.
   */
  One_critical_filtration &operator=(const One_critical_filtration &a) {
    Base::operator=(a);
    return *this;
  }

  // HERITAGE

  using std::vector<T>::operator[]; /**< Inheritance of entry access. */
  using value_type = T;             /**< Entry type. */

  // CONVERTERS

  /**
   * @brief Cast into `std::vector<T> &`.
   */
  operator std::vector<T> &() const { return *this; }

  /**
   * @brief Cast into `std::vector<T>`.
   */
  operator std::vector<T>() const { return static_cast<std::vector<T> >(*this); }

  // like numpy
  /**
   * @brief Returns a copy with entries casted into the type given as template parameter.
   *
   * @tparam U New type for the entries.
   * @return Copy with new entry type.
   */
  template <typename U>
  One_critical_filtration<U> as_type() const {
    One_critical_filtration<U> out(0);
    out.reserve(Base::size());
    for (std::size_t i = 0u; i < Base::size(); i++) out.push_back(static_cast<U>(Base::operator[](i)));
    return out;
  }

  // ACCESS

  /**
   * @brief Returns the number of parameters of the finite filtration value. If the value is "inf", "-inf" or "NaN",
   * returns 1.
   *
   * @return Number of parameters.
   */
  std::size_t num_parameters() const { return Base::size(); }

  /**
   * @brief Returns a filtration value for which @ref is_inf() returns `true`.
   *
   * @return Infinity.
   */
  constexpr static One_critical_filtration inf() { return {T_inf}; }

  /**
   * @brief Returns a filtration value for which @ref is_minus_inf() returns `true`.
   *
   * @return Minus infinity.
   */
  constexpr static One_critical_filtration minus_inf() { return {-T_inf}; }

  /**
   * @brief Returns a filtration value for which @ref is_nan() returns `true`.
   *
   * @return NaN.
   */
  constexpr static One_critical_filtration nan() {
    if constexpr (std::numeric_limits<T>::has_quiet_NaN) {
      return {std::numeric_limits<T>::quiet_NaN()};
    } else {
      return One_critical_filtration(0);  // to differentiate it from 0, an empty filtration value can't do much anyway.
    }
  }

  // DESCRIPTORS

  // TODO: Accept {-inf, -inf, ...} / {inf, inf, ...} / {NaN, NaN, ...} as resp. -inf / inf / NaN.

  /**
   * @brief Returns `true` if and only if the filtration value is considered as infinity.
   */
  bool is_inf() const {
    if (Base::size() != 1) return false;
    return (Base::operator[](0) == T_inf);
  }

  /**
   * @brief Returns `true` if and only if the filtration value is considered as minus infinity.
   */
  bool is_minus_inf() const {
    if (Base::size() != 1) return false;
    return (Base::operator[](0) == -T_inf);
  }

  /**
   * @brief Returns `true` if and only if the filtration value is considered as NaN.
   */
  bool is_nan() const {
    if constexpr (std::numeric_limits<T>::has_quiet_NaN) {
      if (Base::size() != 1) return false;
      return is_nan_(Base::operator[](0));
    } else {
      return Base::empty();
    }
  }

  /**
   * @brief Returns `true` if and only if the filtration value is non-empty and is not considered as infinity,
   * minus infinity or NaN.
   */
  bool is_finite() const {
    if (Base::size() > 1) return true;
    if (Base::size() == 0) return false;
    auto first_value = Base::operator[](0);  // TODO : Maybe check all entries ?
    if (is_nan_(first_value) || first_value == -T_inf || first_value == T_inf) return false;
    return true;
  }

  // COMPARAISON OPERATORS

  /**
   * @brief Returns `true` if and only if for each \f$ i \f$, \f$ a[i] \f$ is strictly smaller than \f$ b[i] \f$.
   * If @p a and @p b are both not infinite or NaN, they have to have the same number of parameters.
   *
   * Note that not all filtration values are comparable. That is, \f$ a < b \f$ and \f$ b < a \f$ returning both false
   * does **not** imply \f$ a == b \f$.
   */
  friend bool operator<(const One_critical_filtration &a, const One_critical_filtration &b) {
    if (a.is_inf() || a.is_nan() || b.is_nan() || b.is_minus_inf()) return false;
    if (b.is_inf() || a.is_minus_inf()) return true;
    bool isSame = true;
    auto n = a.size();
    GUDHI_CHECK(a.size() == b.size(), "Two filtration values with different number of parameters are not comparable.");
    for (auto i = 0u; i < n; ++i) {
      if (a[i] > b[i]) return false;
      if (isSame && a[i] != b[i]) isSame = false;
    }
    return !isSame;
  }

  /**
   * @brief Returns `true` if and only if for each \f$ i \f$, \f$ a[i] \f$ is smaller or equal than \f$ b[i] \f$.
   * If @p a and @p b are both not infinite or NaN, they have to have the same number of parameters.
   *
   * Note that not all filtration values are comparable. That is, \f$ a \le b \f$ and \f$ b \le a \f$ can both return
   * `false`.
   */
  friend bool operator<=(const One_critical_filtration &a, const One_critical_filtration &b) {
    if (a.is_nan() || b.is_nan()) return false;
    if (b.is_inf() || a.is_minus_inf()) return true;
    if (a.is_inf() || b.is_minus_inf()) return false;
    auto n = a.size();
    GUDHI_CHECK(a.size() == b.size(), "Two filtration values with different number of parameters are not comparable.");
    for (std::size_t i = 0u; i < n; ++i) {
      if (a[i] > b[i]) return false;
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
  friend bool operator>(const One_critical_filtration &a, const One_critical_filtration &b) { return b < a; }

  /**
   * @brief Returns `true` if and only if for each \f$ i \f$, \f$ a[i] \f$ is greater or equal than \f$ b[i] \f$.
   * If @p a and @p b are both not infinite or NaN, they have to have the same number of parameters.
   *
   * Note that not all filtration values are comparable. That is, \f$ a \ge b \f$ and \f$ b \ge a \f$ can both return
   * `false`.
   */
  friend bool operator>=(const One_critical_filtration &a, const One_critical_filtration &b) { return b <= a; }

  /**
   * @brief Returns `true` if and only if for each \f$ i \f$, \f$ a[i] \f$ is equal to \f$ b[i] \f$.
   */
  friend bool operator==(const One_critical_filtration &a, const One_critical_filtration &b) {
    if (a.num_parameters() != b.num_parameters()) return false;
    for (auto i = 0u; i < a.num_parameters(); i++) {
      if (a[i] != b[i]) return false;
    }
    return true;
  }

  /**
   * @brief Returns `true` if and only if \f$ a == b \f$ returns `false`.
   */
  friend bool operator!=(const One_critical_filtration &a, const One_critical_filtration &b) { return !(a == b); }

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
  friend One_critical_filtration operator-(const One_critical_filtration &f) {
    One_critical_filtration result(0);
    result.reserve(f.size());
    for (auto val : f) {
      result.push_back(-val);
    }
    return result;
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
  friend One_critical_filtration operator-(One_critical_filtration a, const One_critical_filtration &b) {
    a -= b;
    return a;
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
  friend One_critical_filtration operator-(One_critical_filtration f, const T &val) {
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
  friend One_critical_filtration operator-(const T &val, One_critical_filtration f) {
    // TODO: in one go
    f = -f;
    f += val;
    return f;
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
  friend One_critical_filtration &operator-=(One_critical_filtration &result,
                                             const One_critical_filtration &to_subtract) {
    if (result.empty()) return result;

    if (result.is_nan() || to_subtract.is_nan() || (result.is_inf() && to_subtract.is_inf()) ||
        (result.is_minus_inf() && to_subtract.is_minus_inf())) {
      result = nan();
      return result;
    }
    if (result.is_inf() || to_subtract.is_minus_inf()) {
      result = inf();
      return result;
    }
    if (result.is_minus_inf() || to_subtract.is_inf()) {
      result = minus_inf();
      return result;
    }

    GUDHI_CHECK(result.size() == to_subtract.size(),
                "Two filtration values with different number of parameters cannot be subtracted.");

    return apply_operation_with_finite_values_(result, to_subtract, subtract_);
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
  friend One_critical_filtration &operator-=(One_critical_filtration &result, const T &to_subtract) {
    if (result.empty()) return result;

    if (result.is_nan() || is_nan_(to_subtract) || (result.is_inf() && to_subtract == T_inf) ||
        (result.is_minus_inf() && to_subtract == -T_inf)) {
      result = nan();
      return result;
    }
    if (result.is_inf() || to_subtract == -T_inf) {
      result = inf();
      return result;
    }
    if (result.is_minus_inf() || to_subtract == T_inf) {
      result = minus_inf();
      return result;
    }

    return apply_scalar_operation_on_finite_value_(result, to_subtract, subtract_);
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
  friend One_critical_filtration operator+(One_critical_filtration a, const One_critical_filtration &b) {
    a += b;
    return a;
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
  friend One_critical_filtration operator+(One_critical_filtration f, const T &val) {
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
  friend One_critical_filtration operator+(const T &val, One_critical_filtration f) {
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
  friend One_critical_filtration &operator+=(One_critical_filtration &result, const One_critical_filtration &to_add) {
    if (result.empty()) return result;

    if (result.is_nan() || to_add.is_nan() || (result.is_inf() && to_add.is_minus_inf()) ||
        (result.is_minus_inf() && to_add.is_inf())) {
      result = nan();
      return result;
    }
    if (result.is_inf() || to_add.is_inf()) {
      result = inf();
      return result;
    }
    if (result.is_minus_inf() || to_add.is_minus_inf()) {
      result = minus_inf();
      return result;
    }

    GUDHI_CHECK(result.size() == to_add.size(),
                "Two filtration values with different number of parameters cannot be added.");

    return apply_operation_with_finite_values_(result, to_add, add_);
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
  friend One_critical_filtration &operator+=(One_critical_filtration &result, const T &to_add) {
    if (result.empty()) return result;

    if (result.is_nan() || is_nan_(to_add) || (result.is_inf() && to_add == -T_inf) ||
        (result.is_minus_inf() && to_add == T_inf)) {
      result = nan();
      return result;
    }
    if (result.is_inf() || to_add == T_inf) {
      result = inf();
      return result;
    }
    if (result.is_minus_inf() || to_add == -T_inf) {
      result = minus_inf();
      return result;
    }

    return apply_scalar_operation_on_finite_value_(result, to_add, add_);
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
  friend One_critical_filtration operator*(One_critical_filtration a, const One_critical_filtration &b) {
    a *= b;
    return a;
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
  friend One_critical_filtration operator*(One_critical_filtration f, const T &val) {
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
  friend One_critical_filtration operator*(const T &val, One_critical_filtration f) {
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
  friend One_critical_filtration &operator*=(One_critical_filtration &result, const One_critical_filtration &to_mul) {
    if (result.empty()) return result;

    if (result.is_nan() || to_mul.is_nan()) {
      result = nan();
      return result;
    }

    bool res_is_infinite = result.is_inf() || result.is_minus_inf();
    bool to_mul_is_infinite = to_mul.is_inf() || to_mul.is_minus_inf();

    if (res_is_infinite && to_mul_is_infinite) {
      if (to_mul.is_minus_inf()) {
        result[0] = -result[0];
      }
      return result;
    }

    if (res_is_infinite || to_mul_is_infinite) {
      const One_critical_filtration &finite = res_is_infinite ? to_mul : result;
      const T infinite = res_is_infinite ? result[0] : to_mul[0];
      result = finite;
      return apply_scalar_operation_on_finite_value_(result, infinite, multiply_);
    }

    GUDHI_CHECK(result.size() == to_mul.size(),
                "Two filtration values with different number of parameters cannot be multiplied.");

    return apply_operation_with_finite_values_(result, to_mul, multiply_);
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
  friend One_critical_filtration &operator*=(One_critical_filtration &result, const T &to_mul) {
    if (result.empty()) return result;

    if (result.is_nan() || is_nan_(to_mul)) {
      result = nan();
      return result;
    }

    if (result.is_inf() || result.is_minus_inf()) {
      if (to_mul == 0)
        result = nan();
      else if (to_mul < 0)
        result[0] = -result[0];
      return result;
    }

    return apply_scalar_operation_on_finite_value_(result, to_mul, multiply_);
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
  friend One_critical_filtration operator/(One_critical_filtration a, const One_critical_filtration &b) {
    a /= b;
    return a;
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
  friend One_critical_filtration operator/(One_critical_filtration f, const T &val) {
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
  friend One_critical_filtration operator/(const T &val, const One_critical_filtration &f) {
    if (f.empty()) return f;
    if (is_nan_(val) || f.is_nan()) return nan();

    One_critical_filtration result(f.size(), val);
    result /= f;
    return result;
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
  friend One_critical_filtration &operator/=(One_critical_filtration &result, const One_critical_filtration &to_div) {
    if (result.empty()) return result;

    bool res_is_infinite = result.is_inf() || result.is_minus_inf();
    bool to_div_is_infinite = to_div.is_inf() || to_div.is_minus_inf();

    if (result.is_nan() || to_div.is_nan() || (res_is_infinite && to_div_is_infinite)) {
      result = nan();
      return result;
    }

    if (to_div_is_infinite) {
      return apply_scalar_operation_on_finite_value_with_all_nan_possible_(result, to_div[0], divide_);
    }

    GUDHI_CHECK(res_is_infinite || result.size() == to_div.size(),
                "Two filtration values with different number of parameters cannot be divided.");

    if (res_is_infinite) {
      result.resize(to_div.size(), result[0]);
    }

    return apply_operation_with_finite_values_(result, to_div, divide_);
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
  friend One_critical_filtration &operator/=(One_critical_filtration &result, const T &to_div) {
    if (result.empty()) return result;

    bool res_is_infinite = result.is_inf() || result.is_minus_inf();
    bool to_div_is_infinite = to_div == T_inf || to_div == -T_inf;

    if (to_div == 0 || is_nan_(to_div) || result.is_nan() || (res_is_infinite && to_div_is_infinite)) {
      result = nan();
      return result;
    }

    if (res_is_infinite) {
      if (to_div < 0) result[0] = -result[0];
      return result;
    }

    return apply_scalar_operation_on_finite_value_with_all_nan_possible_(result, to_div, divide_);
  }

  // MODIFIERS

  /**
   * @brief Sets the filtration value to the least common upper bound between the current value and the given value.
   *
   * More formally, it pushes the current filtration value to the cone \f$ \{ y \in \mathbb R^n : y \ge x \} \f$
   * originating in the given filtration value \f$ x \f$. The resulting value corresponds to the intersection of both
   * cones: \f$ \mathrm{this} = \min \{ y \in \mathbb R^n : y \ge this \} \cap \{ y \in \mathbb R^n : y \ge x \} \f$.
   *
   * @param x The target filtration value towards which to push.
   */
  void push_to_least_common_upper_bound(const One_critical_filtration &x) {
    if (this->is_inf() || this->is_nan() || x.is_nan() || x.is_minus_inf()) return;
    if (x.is_inf() || this->is_minus_inf()) {
      *this = x;
      return;
    }

    GUDHI_CHECK(this->num_parameters() == x.num_parameters(),
                "A filtration value cannot be pushed to another one with different numbers of parameters.");

    for (std::size_t i = 0; i < x.num_parameters(); i++)
      Base::operator[](i) = Base::operator[](i) > x[i] ? Base::operator[](i) : x[i];
  }

  /**
   * @brief Sets the filtration value to the greatest common lower bound between the current value and the given value.
   *
   * More formally, it pulls the current filtration value to the cone \f$ \{ y \in \mathbb R^n : y \le x \} \f$
   * originating in the given filtration value \f$ x \f$. The resulting value corresponds to the intersection of both
   * cones: \f$ \mathrm{this} = \min \{ y \in \mathbb R^n : y \le this \} \cap \{ y \in \mathbb R^n : y \le x \} \f$.
   *
   * @param x The target filtration value towards which to pull.
   */
  void pull_to_greatest_common_lower_bound(const One_critical_filtration &x) {
    if (x.is_inf() || this->is_nan() || x.is_nan() || this->is_minus_inf()) return;
    if (this->is_inf() || x.is_minus_inf()) {
      *this = x;
      return;
    }

    GUDHI_CHECK(this->num_parameters() == x.num_parameters(),
                "A filtration value cannot be pulled to another one with different numbers of parameters.");

    for (std::size_t i = 0u; i < x.num_parameters(); i++)
      Base::operator[](i) = Base::operator[](i) > x[i] ? x[i] : Base::operator[](i);
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
   * @tparam oned_array A range of values convertible into `T` ordered by increasing value. Has to implement
   * a begin, end and operator[] method.
   * @param grid Vector of @p oned_array with size at least number of filtration parameters.
   * @param coordinate If true, the values are set to the coordinates of the projection in the grid. If false,
   * the values are set to the values at the coordinates of the projection.
   */
  template <typename oned_array>
  void project_onto_grid(const std::vector<oned_array> &grid, bool coordinate = true) {
    GUDHI_CHECK(grid.size() >= Base::size(),
                "The grid should not be smaller than the number of parameters in the filtration value.");
    for (std::size_t parameter = 0u; parameter < Base::size(); ++parameter) {
      const auto &filtration = grid[parameter];
      auto d =
          std::distance(filtration.begin(),
                        std::lower_bound(filtration.begin(), filtration.end(),
                                         static_cast<typename oned_array::value_type>(Base::operator[](parameter))));
      Base::operator[](parameter) = coordinate ? static_cast<T>(d) : static_cast<T>(filtration[d]);
    }
  }

  // FONCTIONNALITIES

  /**
   * @brief Computes the scalar product of the given filtration value with the given vector.
   *
   * @tparam U Arithmetic type of the result. Default value: `T`.
   * @param f Filtration value.
   * @param x Vector of coefficients.
   * @return Scalar product of @p f with @p x.
   */
  template <typename U = T>
  friend U compute_linear_projection(const One_critical_filtration &f, const std::vector<U> &x) {
    U projection = 0;
    std::size_t size = std::min(x.size(), f.size());
    for (std::size_t i = 0u; i < size; i++) projection += x[i] * static_cast<U>(f[i]);
    return projection;
  }

  /**
   * @brief Computes the norm of the given filtration value.
   *
   * @param f Filtration value.
   * @return The norm of @p f.
   */
  friend T compute_norm(const One_critical_filtration &f) {
    T out = 0;
    for (auto &val : f) out += (val * val);
    if constexpr (std::is_integral_v<T>){
       //to avoid Windows issue which don't know how to cast integers for cmath methods
      return std::sqrt(static_cast<double>(out));
    } else {
      return std::sqrt(out);
    }
  }

  /**
   * @brief Computes the euclidean distance from the first parameter to the second parameter.
   *
   * @param f Start filtration value.
   * @param other End filtration value.
   * @return Euclidean distance between @p f and @p other.
   */
  friend T compute_euclidean_distance_to(const One_critical_filtration &f, const One_critical_filtration &other) {
    T out = 0;
    for (std::size_t i = 0u; i < other.size(); i++) {
      out += (f[i] - other[i]) * (f[i] - other[i]);
    }
    if constexpr (std::is_integral_v<T>){
       //to avoid Windows issue which don't know how to cast integers for cmath methods
      return std::sqrt(static_cast<double>(out));
    } else {
      return std::sqrt(out);
    }
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
  template <typename out_type = std::int32_t, typename U = T>
  friend One_critical_filtration<out_type> compute_coordinates_in_grid(One_critical_filtration f,
                                                                       const std::vector<std::vector<U> > &grid) {
    // TODO: by replicating the code of "project_onto_grid", this could be done with just one copy
    // instead of two. But it is not clear if it is really worth it, i.e., how much the change in type is really
    // necessary in the use cases. To see later.
    f.project_onto_grid(grid);
    if constexpr (std::is_same_v<out_type, T>) {
      return f;
    } else {
      return f.as_type<out_type>();
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
  friend One_critical_filtration<U> evaluate_coordinates_in_grid(const One_critical_filtration &f,
                                                                 const std::vector<std::vector<U> > &grid) {
    One_critical_filtration<U> pushed_value(f.size());

    GUDHI_CHECK(grid.size() == f.size(),
                "The size of the grid should correspond to the number of parameters in the filtration value.");

    U grid_inf = One_critical_filtration<U>::T_inf;

    for (std::size_t parameter = 0u; parameter < grid.size(); ++parameter) {
      const auto &filtration = grid[parameter];
      const auto &c = f[parameter];
      pushed_value[parameter] = c == f.T_inf ? grid_inf : filtration[c];
    }
    return pushed_value;
  }

  // UTILITIES

  /**
   * @brief Outstream operator.
   */
  friend std::ostream &operator<<(std::ostream &stream, const One_critical_filtration &f) {
    if (f.is_inf()) {
      stream << "[inf, ..., inf]";
      return stream;
    }
    if (f.is_minus_inf()) {
      stream << "[-inf, ..., -inf]";
      return stream;
    }
    if (f.is_nan()) {
      stream << "[NaN]";
      return stream;
    }
    if (f.empty()) {
      stream << "[]";
      return stream;
    }
    stream << "[";
    for (std::size_t i = 0; i < f.size() - 1; i++) {
      stream << f[i] << ", ";
    }
    if (!f.empty()) stream << f.back();
    stream << "]";
    return stream;
  }

 public:
  /**
   * @brief Infinity value of an entry of the filtration value.
   */
  constexpr static const T T_inf =
      std::numeric_limits<T>::has_infinity ? std::numeric_limits<T>::infinity() : std::numeric_limits<T>::max();

  /**
   * @brief Indicates if the class manages multi-critical filtration values.
   */
  constexpr static bool is_multi_critical = false;

 private:
  static bool is_nan_(T val){
    if constexpr (std::is_integral_v<T>){
       //to avoid Windows issue which don't know how to cast integers for cmath methods
      return false;
    } else {
      return std::isnan(val);
    }
  }

  constexpr static bool subtract_(T &v1, T v2) { return add_(v1, -v2); }

  constexpr static bool add_(T &v1, T v2) {
    if (is_nan_(v1) || is_nan_(v2) || (v1 == T_inf && v2 == -T_inf) || (v1 == -T_inf && v2 == T_inf)) {
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

  constexpr static bool multiply_(T &v1, T v2) {
    bool v1_is_infinite = v1 == T_inf || v1 == -T_inf;
    bool v2_is_infinite = v2 == T_inf || v2 == -T_inf;

    if (is_nan_(v1) || is_nan_(v2) || (v1_is_infinite && v2 == 0) || (v1 == 0 && v2_is_infinite)) {
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

  constexpr static bool divide_(T &v1, T v2) {
    bool v1_is_infinite = v1 == T_inf || v1 == -T_inf;
    bool v2_is_infinite = v2 == T_inf || v2 == -T_inf;

    if (is_nan_(v1) || is_nan_(v2) || v2 == 0 || (v1_is_infinite && v2_is_infinite)) {
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

  constexpr static bool update_sign_(T toComp, int &sign) {
    if (toComp == T_inf) {
      if (sign == 0)
        sign = 1;
      else if (sign == -1)
        return false;
    } else if (toComp == -T_inf) {
      if (sign == 0)
        sign = -1;
      else if (sign == 1)
        return false;
    } else {
      return false;
    }

    return true;
  }

  template <typename F>
  static One_critical_filtration &apply_operation_with_finite_values_(One_critical_filtration &result,
                                                                      const One_critical_filtration &to_operate,
                                                                      F &&operate) {
    bool allSameInf = true;
    bool allNaN = true;
    int sign = 0;
    for (auto i = 0u; i < result.size(); ++i) {
      if (operate(result[i], to_operate[i])) {
        allNaN = false;
      } else {
        if constexpr (!std::numeric_limits<T>::has_quiet_NaN) {
          result = nan();
          return result;
        }
      }
      if (allSameInf) allSameInf = update_sign_(result[i], sign);
    }

    if (allSameInf) result = (sign == 1 ? inf() : minus_inf());
    if (allNaN) result = nan();

    return result;
  }

  template <typename F>
  static One_critical_filtration &apply_scalar_operation_on_finite_value_(One_critical_filtration &result,
                                                                          const T &to_operate, F &&operate) {
    for (auto &val : result) {
      if constexpr (std::numeric_limits<T>::has_quiet_NaN) {
        operate(val, to_operate);
      } else {
        if (!operate(val, to_operate)) {
          result = nan();
          return result;
        }
      }
    }

    return result;
  }

  template <typename F>
  static One_critical_filtration &apply_scalar_operation_on_finite_value_with_all_nan_possible_(
      One_critical_filtration &result, const T &to_operate, F &&operate) {
    bool allNaN = true;

    for (auto &val : result) {
      if (operate(val, to_operate)) {
        allNaN = false;
      } else {
        if constexpr (!std::numeric_limits<T>::has_quiet_NaN) {
          result = nan();
          return result;
        }
      }
    }
    if (allNaN) result = nan();

    return result;
  }
};

}  // namespace Gudhi::multi_filtration

namespace std {

template <typename T>
class numeric_limits<Gudhi::multi_filtration::One_critical_filtration<T> > {
 public:
  static constexpr bool has_infinity = true;

  static constexpr Gudhi::multi_filtration::One_critical_filtration<T> infinity() noexcept {
    return Gudhi::multi_filtration::One_critical_filtration<T>::inf();
  };

  // non-standard
  static constexpr Gudhi::multi_filtration::One_critical_filtration<T> minus_infinity() noexcept {
    return Gudhi::multi_filtration::One_critical_filtration<T>::minus_inf();
  };

  static constexpr Gudhi::multi_filtration::One_critical_filtration<T> max() noexcept(false) {
    throw std::logic_error(
        "The maximal value cannot be represented with no finite numbers of parameters."
        "Use `max(number_of_parameters)` instead");
  };

  // non-standard, so I don't want to define a default value.
  static constexpr Gudhi::multi_filtration::One_critical_filtration<T> max(unsigned int n) noexcept {
    return Gudhi::multi_filtration::One_critical_filtration<T>(n, std::numeric_limits<T>::max());
  };

  static constexpr Gudhi::multi_filtration::One_critical_filtration<T> quiet_NaN() noexcept {
    return Gudhi::multi_filtration::One_critical_filtration<T>::nan();
  };
};

}  // namespace std

#endif  // ONE_CRITICAL_FILTRATIONS_H_
