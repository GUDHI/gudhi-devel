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
#include <gudhi/Simple_mdspan.h>
#include <gudhi/multi_filtration_utils.h>

namespace Gudhi::multi_filtration {

// // declaration needed pre C++20 for friends with templates defined inside a class
// template <typename U>
// U compute_euclidean_distance_to();
// template <typename U>
// U compute_norm();

/**
 * @class Multi_parameter_generator Multi_parameter_generator.h gudhi/Multi_parameter_generator.h
 * @ingroup multi_filtration
 *
 * @brief Class encoding the different generators, i.e., apparition times, of a \f$ k \f$-critical
 * \f$\mathbb R^n\f$-filtration value. E.g., the filtration value of a simplex, or, of the algebraic generator of a
 * module presentation. The encoding is compacted into a single vector, so if a lot of non trivial modifications are
 * done (that not only consists of simply adding new generators at the end of the vector), it is probably preferable
 * to use @ref Dynamic_multi_parameter_filtration instead.
 *
 * @details Overloads `std::numeric_limits` such that:
 * - `std::numeric_limits<Multi_parameter_generator>::has_infinity` returns `true`,
 * - `std::numeric_limits<Multi_parameter_generator>::has_quiet_NaN` returns `std::numeric_limits<T>::has_quiet_NaN`,
 * - `std::numeric_limits<Multi_parameter_generator>::infinity(num_param)` returns
 * @ref Multi_parameter_generator::inf(num_param) "",
 * - `std::numeric_limits<Multi_parameter_generator>::minus_infinity(num_param)` returns
 * @ref Multi_parameter_generator::minus_inf(num_param) "",
 * - `std::numeric_limits<Multi_parameter_generator>::max(num_param)` returns a @ref Multi_parameter_generator with
 * 1 generators of `num_param` parameters evaluated at value `std::numeric_limits<T>::max()`,
 * - `std::numeric_limits<Multi_parameter_generator>::quiet_NaN(num_param)` returns
 * @ref Multi_parameter_generator::nan(num_param) if `std::numeric_limits<Multi_parameter_generator>::has_quiet_NaN`
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
template <typename T, bool Co = false>
class Multi_parameter_generator
{
 public:
  using Underlying_container = std::vector<T>; /**< Underlying container for values. */

  // CONSTRUCTORS

  /**
   * @brief Default constructor. Builds filtration value with one generator and given number of parameters.
   * If Co is false, all values are at -inf, if Co is true, all values are at +inf.
   *
   * @param number_of_parameters If negative, takes the default value instead. Default value: 2.
   */
  Multi_parameter_generator(int number_of_parameters = 1)
      : generator_(number_of_parameters < 0 ? 1 : number_of_parameters, _get_default_value())
  {}

  /**
   * @brief Builds filtration value with one generator and given number of parameters.
   * All values are initialized at the given value.
   *
   * @param number_of_parameters If negative, is set to 2 instead.
   * @param value Initialization value for every value in the generator.
   */
  Multi_parameter_generator(int number_of_parameters, T value)
      : generator_(number_of_parameters < 0 ? 2 : number_of_parameters, value)
  {}

  /**
   * @brief Builds filtration value with one generator that is initialized with the given range. The number of
   * parameters are therefore deduced from the length of the range.
   *
   * @tparam ValueRange Range of types convertible to `T`. Should have a begin() and end() method.
   * @param range Values of the generator.
   */
  template <class ValueRange = std::initializer_list<T>, class = std::enable_if_t<RangeTraits<ValueRange>::has_begin> >
  Multi_parameter_generator(const ValueRange &range) : generator_(range.begin(), range.end())
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
  Multi_parameter_generator(Iterator it_begin, Iterator it_end) : generator_(it_begin, it_end)
  {}

  /**
   * @brief Builds filtration value with given number of parameters and values from the given range. Lets \f$ p \f$
   * be the number of parameters. The \f$ p \f$ first elements of the range have to correspond to the first generator,
   * the \f$ p \f$ next elements to the second generator and so on... So the length of the range has to be a multiple
   * of \f$ p \f$ and the number of generators will be \f$ length / p \f$. The range is represented by
   * @ref Multi_parameter_generator::Underlying_container "" and **moved** into the underlying container of the class.
   *
   * @param generators Values to move.
   * @param number_of_parameters Negative values are associated to 0.
   */
  Multi_parameter_generator(Underlying_container &&generators) : generator_(std::move(generators)) {}

  /**
   * @brief Copy constructor.
   */
  Multi_parameter_generator(const Multi_parameter_generator &other) : generator_(other.generator_) {}

  /**
   * @brief Copy constructor.
   *
   * @tparam U Type convertible into `T`.
   */
  template <typename U>
  Multi_parameter_generator(const Multi_parameter_generator<U, Co> &other) : generator_(other.begin(), other.end())
  {}

  /**
   * @brief Assign operator.
   */
  Multi_parameter_generator &operator=(const Multi_parameter_generator &other)
  {
    generator_ = other.generator_;
    return *this;
  }

  /**
   * @brief Assign operator.
   *
   * @tparam U Type convertible into `T`.
   */
  template <typename U>
  Multi_parameter_generator &operator=(const Multi_parameter_generator<U, Co> &other)
  {
    generator_ = Underlying_container(other.begin(), other.end());
    return *this;
  }

  /**
   * @brief Swap operator.
   */
  friend void swap(Multi_parameter_generator &f1, Multi_parameter_generator &f2) noexcept
  {
    f1.generator_.swap(f2.generator_);
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
   * @brief Let \f$ g \f$ be the first value in `indices` and \f$ p \f$ the second value.
   * Returns reference to value of parameter \f$ p \f$ of generator \f$ g \f$.
   *
   * @tparam IndexRange Range with a begin() and size() method.
   * @param indices Range with at least two elements. The first element should correspond to the generator number and
   * the second element to the parameter number.
   */
  reference operator[](size_type p)
  {
    GUDHI_CHECK(p < generator_.size(), "Out of bound parameter.");
    return generator_[p];
  }

  /**
   * @brief Let \f$ g \f$ be the first value in `indices` and \f$ p \f$ the second value.
   * Returns reference to value of parameter \f$ p \f$ of generator \f$ g \f$.
   *
   * @tparam IndexRange Range with a begin() and size() method.
   * @param indices Range with at least two elements. The first element should correspond to the generator number and
   * the second element to the parameter number.
   */
  const_reference operator[](size_type p) const
  {
    GUDHI_CHECK(p < generator_.size(), "Out of bound parameter.");
    return generator_[p];
  }

  /**
   * @brief Returns an iterator pointing the begining of the underlying container. The @ref num_parameter first elements
   * corresponds to the first generator, the @ref num_parameter next to the second and so on.
   */
  iterator begin() noexcept { return generator_.begin(); }

  /**
   * @brief Returns an iterator pointing the begining of the underlying container. The @ref num_parameter first elements
   * corresponds to the first generator, the @ref num_parameter next to the second and so on.
   */
  const_iterator begin() const noexcept { return generator_.begin(); }

  /**
   * @brief Returns an iterator pointing the begining of the underlying container. The @ref num_parameter first elements
   * corresponds to the first generator, the @ref num_parameter next to the second and so on.
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
   * The @ref num_parameter first elements corresponds to the last generator (in parameter reverse order), the
   * @ref num_parameter next to the second to last and so on.
   */
  reverse_iterator rbegin() noexcept { return generator_.rbegin(); }

  /**
   * @brief Returns a reverse iterator pointing to the first element from the back of the underlying container.
   * The @ref num_parameter first elements corresponds to the last generator (in parameter reverse order), the
   * @ref num_parameter next to the second to last and so on.
   */
  const_reverse_iterator rbegin() const noexcept { return generator_.rbegin(); }

  /**
   * @brief Returns a reverse iterator pointing to the first element from the back of the underlying container.
   * The @ref num_parameter first elements corresponds to the last generator (in parameter reverse order), the
   * @ref num_parameter next to the second to last and so on.
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
   * @brief Returns the size of the underlying container. Corresponds exactly to @ref num_entries(), but enables to use
   * the class as a classic range with a `begin`, `end` and `size` method.
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
  Multi_parameter_generator<U, Co> as_type() const
  {
    std::vector<U> out(generator_.begin(), generator_.end());
    return Multi_parameter_generator<U, Co>(std::move(out));
  }

  // ACCESS

  /**
   * @brief Returns the number of parameters in the filtration value.
   */
  size_type num_parameters() const { return generator_.size(); }

  /**
   * @brief Returns a filtration value with given number of parameters for which @ref is_plus_inf() returns `true`.
   */
  static Multi_parameter_generator inf() { return Multi_parameter_generator(1, T_inf); }

  /**
   * @brief Returns a filtration value with given number of parameters for which @ref is_minus_inf() returns `true`.
   */
  static Multi_parameter_generator minus_inf() { return Multi_parameter_generator(1, -T_inf); }

  /**
   * @brief If `std::numeric_limits<T>::has_quiet_NaN` is true, returns a filtration value with given number of
   * parameters for which @ref is_nan() returns `true`. Otherwise, throws.
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
   * @brief Returns `true` if and only if the filtration value is considered as plus infinity.
   */
  bool is_plus_inf() const
  {
    for (const T &v : generator_) {
      if (v != T_inf) return false;
    }
    return !generator_.empty();
  }

  /**
   * @brief Returns `true` if and only if the filtration value is considered as minus infinity.
   */
  bool is_minus_inf() const
  {
    for (const T &v : generator_) {
      if (v != -T_inf) return false;
    }
    return !generator_.empty();
  }

  /**
   * @brief Returns `true` if and only if the filtration value is considered as NaN.
   */
  bool is_nan() const
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
   * @brief Returns `true` if and only if the filtration value is non-empty and is not considered as plus infinity,
   * minus infinity or NaN.
   */
  bool is_finite() const
  {
    bool isInf = true, isMinusInf = true, isNan = true;
    for (const auto &v : generator_) {
      if (v != T_inf) isInf = false;
      if (v != -T_inf) isMinusInf = false;
      if (!_is_nan(v)) isNan = false;
      if (!isInf && !isMinusInf && !isNan) return true;
    }
    return false;
  }

  // COMPARAISON OPERATORS

  /**
   * @brief Returns `true` if and only if the cones generated by @p b are strictly contained in the
   * cones generated by @p a (recall that the cones are positive if `Co` is false and negative if `Co` is true).
   * Both @p a and @p b  have to have the same number of parameters.
   *
   * Note that not all filtration values are comparable. That is, \f$ a < b \f$ and \f$ b < a \f$ returning both false
   * does **not** imply \f$ a == b \f$.
   */
  friend bool operator<(const Multi_parameter_generator &a, const Multi_parameter_generator &b)
  {
    if (a.is_nan() || b.is_nan()) return false;
    if (a.is_minus_inf() || b.is_plus_inf()) return !Co;
    if (b.is_minus_inf() || a.is_plus_inf()) return Co;

    GUDHI_CHECK(a.num_parameters() == b.num_parameters(),
                "Only generators with same number of parameters can be compared.");

    bool isSame = true;
    auto n = a.size();
    for (size_type i = 0u; i < n; ++i) {
      if (a[i] > b[i]) return false;
      if (isSame && a[i] != b[i]) isSame = false;
    }
    return !isSame;
  }

  /**
   * @brief Returns `true` if and only if the cones generated by @p a are strictly contained in the
   * cones generated by @p b (recall that the cones are positive if `Co` is false and negative if `Co` is true).
   * Both @p a and @p b  have to have the same number of parameters.
   *
   * Note that not all filtration values are comparable. That is, \f$ a \le b \f$ and \f$ b \le a \f$ can both return
   * `false`.
   */
  friend bool operator<=(const Multi_parameter_generator &a, const Multi_parameter_generator &b)
  {
    if (a.is_nan() || b.is_nan()) return false;
    const bool aIsMinusInf = a.is_minus_inf();
    const bool aIsPlusInf = a.is_plus_inf();
    const bool bIsMinusInf = b.is_minus_inf();
    const bool bIsPlusInf = b.is_plus_inf();
    if ((aIsMinusInf && bIsMinusInf) || (aIsPlusInf && bIsPlusInf)) return true;
    if (aIsMinusInf || bIsPlusInf) return !Co;
    if (bIsMinusInf || aIsPlusInf) return Co;

    GUDHI_CHECK(a.num_parameters() == b.num_parameters(),
                "Only filtration values with same number of parameters can be compared.");

    auto n = a.size();
    for (std::size_t i = 0u; i < n; ++i) {
      if (a[i] > b[i]) return false;
    }
    return true;
  }

  /**
   * @brief Returns `true` if and only if the cones generated by @p b are contained in or are (partially)
   * equal to the cones generated by @p a (recall that the cones are positive if `Co` is false and negative if `Co` is
   * true).
   * Both @p a and @p b  have to have the same number of parameters.
   *
   * Note that not all filtration values are comparable. That is, \f$ a > b \f$ and \f$ b > a \f$ returning both false
   * does **not** imply \f$ a == b \f$.
   */
  friend bool operator>(const Multi_parameter_generator &a, const Multi_parameter_generator &b) { return b < a; }

  /**
   * @brief Returns `true` if and only if the cones generated by @p a are contained in or are (partially)
   * equal to the cones generated by @p b (recall that the cones are positive if `Co` is false and negative if `Co` is
   * true).
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
    return a.generator_ == b.generator_;
  }

  /**
   * @brief Returns `true` if and only if \f$ a == b \f$ returns `false`.
   */
  friend bool operator!=(const Multi_parameter_generator &a, const Multi_parameter_generator &b) { return !(a == b); }

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
  friend Multi_parameter_generator operator-(const Multi_parameter_generator &f)
  {
    std::vector<T> result(f.generator_);
    std::for_each(result.begin(), result.end(), [](T &v) { v = -v; });
    return Multi_parameter_generator(std::move(result));
  }

  // subtraction
  /**
   * @brief Returns a filtration value such that an entry at index \f$ g,p \f$, with
   * \f$ 0 \leq g \leq num_generators \f$ and \f$ 0 \leq p \leq num_parameters \f$ is equal to \f$ f(g,p) - r(p) \f$
   * if \f$ p \leq length_r \f$ and to \f$ f(g,p) \f$ otherwise.
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
   * @param f First element of the subtraction.
   * @param r Second element of the subtraction.
   */
  template <class ValueRange, class = std::enable_if_t<RangeTraits<ValueRange>::has_begin> >
  friend Multi_parameter_generator operator-(Multi_parameter_generator f, const ValueRange &r)
  {
    f -= r;
    return f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ g,p \f$, with
   * \f$ 0 \leq g \leq num_generators \f$ and \f$ 0 \leq p \leq num_parameters \f$ is equal to \f$ r(p) - f(g,p) \f$
   * if \f$ p \leq length_r \f$ and to \f$ -f(g,p) \f$ otherwise.
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
   * @param f Second element of the subtraction.
   */
  template <class ValueRange,
            class = std::enable_if_t<RangeTraits<ValueRange>::has_begin &&
                                     !std::is_same_v<ValueRange, Multi_parameter_generator> > >
  friend Multi_parameter_generator operator-(const ValueRange &r, const Multi_parameter_generator &f)
  {
    if (f.num_parameters() == 0) return f;
    if (f.is_nan()) return nan();

    if constexpr (RangeTraits<ValueRange>::is_dynamic_multi_filtration){
      if (r[0].is_finite()){
        return f -= r[0];
      } else {
        if constexpr (!std::numeric_limits<T>::has_quiet_NaN) if (r[0].size() == 0) return nan();
        return f -= r[0][0];
      }
    } else {
      if (f.is_finite()) {
        Multi_parameter_generator res = f;
        res._apply_operation(r, [](T &valF, const T &valR) -> bool {
          valF = -valF;
          return _add(valF, valR);
        });
        return res;
      } else {
        Multi_parameter_generator res(r.begin(), r.end());
        res._apply_operation(f[0], [](T &valRes, const T &valF) -> bool { return _subtract(valRes, valF); });
        return res;
      }
    }
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ g,p \f$ is equal to \f$ f(g,p) - val \f$.
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
  friend Multi_parameter_generator operator-(Multi_parameter_generator f, const T &val)
  {
    f -= val;
    return f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ g,p \f$ is equal to \f$ val - f(g,p) \f$.
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
  friend Multi_parameter_generator operator-(const T &val, Multi_parameter_generator f)
  {
    f._apply_operation(val, [](T &valF, const T &valR) -> bool {
      valF = -valF;
      return _add(valF, valR);
    });
    return f;
  }

  /**
   * @brief Modifies the first parameter such that an entry at index \f$ g,p \f$, with
   * \f$ 0 \leq g \leq num_generators \f$ and \f$ 0 \leq p \leq num_parameters \f$ is equal to \f$ f(g,p) - r(p) \f$
   * if \f$ p \leq length_r \f$ and to \f$ f(g,p) \f$ otherwise.
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
   * @param f First element of the subtraction.
   * @param r Second element of the subtraction.
   */
  template <class ValueRange, class = std::enable_if_t<RangeTraits<ValueRange>::has_begin> >
  friend Multi_parameter_generator &operator-=(Multi_parameter_generator &f, const ValueRange &r)
  {
    if (f.num_parameters() == 0 || f.is_nan()) return f;

    if constexpr (RangeTraits<ValueRange>::is_dynamic_multi_filtration){
      if (r[0].is_finite()){
        return f -= r[0];
      } else {
        if constexpr (!std::numeric_limits<T>::has_quiet_NaN)
          if (r[0].size() == 0) {
            f = nan();
            return f;
          }
        return f -= r[0][0];
      }
    } else {
      if (!f.is_finite()) {
        f.generator_.resize(r.size(), f[0]);
      }
      f._apply_operation(r, [](T &valF, const T &valR) -> bool { return _subtract(valF, valR); });
    }
    
    return f;
  }

  /**
   * @brief Modifies the first parameter such that an entry at index \f$ g,p \f$ is equal to \f$ f(g,p) - val \f$.
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
  friend Multi_parameter_generator &operator-=(Multi_parameter_generator &f, const T &val)
  {
    f._apply_operation(val, [](T &valF, const T &valR) -> bool { return _subtract(valF, valR); });
    return f;
  }

  // addition
  /**
   * @brief Returns a filtration value such that an entry at index \f$ g,p \f$, with
   * \f$ 0 \leq g \leq num_generators \f$ and \f$ 0 \leq p \leq num_parameters \f$ is equal to \f$ f(g,p) + r(p) \f$
   * if \f$ p \leq length_r \f$ and to \f$ f(g,p) \f$ otherwise.
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
   * @param f First element of the addition.
   * @param r Second element of the addition.
   */
  template <class ValueRange, class = std::enable_if_t<RangeTraits<ValueRange>::has_begin> >
  friend Multi_parameter_generator operator+(Multi_parameter_generator f, const ValueRange &r)
  {
    f += r;
    return f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ g,p \f$, with
   * \f$ 0 \leq g \leq num_generators \f$ and \f$ 0 \leq p \leq num_parameters \f$ is equal to \f$ r(p) + f(g,p) \f$
   * if \f$ p \leq length_r \f$ and to \f$ f(g,p) \f$ otherwise.
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
   * @param f Second element of the addition.
   */
  template <class ValueRange,
            class = std::enable_if_t<RangeTraits<ValueRange>::has_begin &&
                                     !std::is_same_v<ValueRange, Multi_parameter_generator> > >
  friend Multi_parameter_generator operator+(const ValueRange &r, Multi_parameter_generator f)
  {
    f += r;
    return f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ g,p \f$ is equal to \f$ f(g,p) + val \f$.
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
  friend Multi_parameter_generator operator+(Multi_parameter_generator f, const T &val)
  {
    f += val;
    return f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ g,p \f$ is equal to \f$ val + f(g,p) \f$.
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
  friend Multi_parameter_generator operator+(const T &val, Multi_parameter_generator f)
  {
    f += val;
    return f;
  }

  /**
   * @brief Modifies the first parameter such that an entry at index \f$ g,p \f$, with
   * \f$ 0 \leq g \leq num_generators \f$ and \f$ 0 \leq p \leq num_parameters \f$ is equal to \f$ f(g,p) + r(p) \f$
   * if \f$ p \leq length_r \f$ and to \f$ f(g,p) \f$ otherwise.
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
   * @param f First element of the addition.
   * @param r Second element of the addition.
   */
  template <class ValueRange, class = std::enable_if_t<RangeTraits<ValueRange>::has_begin> >
  friend Multi_parameter_generator &operator+=(Multi_parameter_generator &f, const ValueRange &r)
  {
    if (f.num_parameters() == 0 || f.is_nan()) return f;

    if constexpr (RangeTraits<ValueRange>::is_dynamic_multi_filtration){
      if (r[0].is_finite()){
        return f += r[0];
      } else {
        if constexpr (!std::numeric_limits<T>::has_quiet_NaN)
          if (r[0].size() == 0) {
            f = nan();
            return f;
          }
        return f += r[0][0];
      }
    } else {
      if (!f.is_finite()) {
        f.generator_.resize(r.size(), f[0]);
      }
      f._apply_operation(r, [](T &valF, const T &valR) -> bool { return _add(valF, valR); });
    }

    return f;
  }

  /**
   * @brief Modifies the first parameter such that an entry at index \f$ g,p \f$ is equal to \f$ f(g,p) + val \f$.
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
  friend Multi_parameter_generator &operator+=(Multi_parameter_generator &f, const T &val)
  {
    f._apply_operation(val, [](T &valF, const T &valR) -> bool { return _add(valF, valR); });
    return f;
  }

  // multiplication
  /**
   * @brief Returns a filtration value such that an entry at index \f$ g,p \f$, with
   * \f$ 0 \leq g \leq num_generators \f$ and \f$ 0 \leq p \leq num_parameters \f$ is equal to \f$ f(g,p) * r(p) \f$
   * if \f$ p \leq length_r \f$ and to \f$ f(g,p) \f$ otherwise.
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
   * @param f First element of the multiplication.
   * @param r Second element of the multiplication.
   */
  template <class ValueRange, class = std::enable_if_t<RangeTraits<ValueRange>::has_begin> >
  friend Multi_parameter_generator operator*(Multi_parameter_generator f, const ValueRange &r)
  {
    f *= r;
    return f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ g,p \f$, with
   * \f$ 0 \leq g \leq num_generators \f$ and \f$ 0 \leq p \leq num_parameters \f$ is equal to \f$ r(p) * f(g,p) \f$
   * if \f$ p \leq length_r \f$ and to \f$ f(g,p) \f$ otherwise.
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
   * @param f Second element of the multiplication.
   */
  template <class ValueRange,
            class = std::enable_if_t<RangeTraits<ValueRange>::has_begin &&
                                     !std::is_same_v<ValueRange, Multi_parameter_generator> > >
  friend Multi_parameter_generator operator*(const ValueRange &r, Multi_parameter_generator f)
  {
    f *= r;
    return f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ g,p \f$ is equal to \f$ f(g,p) * val \f$.
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
  friend Multi_parameter_generator operator*(Multi_parameter_generator f, const T &val)
  {
    f *= val;
    return f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ g,p \f$ is equal to \f$ val * f(g,p) \f$.
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
  friend Multi_parameter_generator operator*(const T &val, Multi_parameter_generator f)
  {
    f *= val;
    return f;
  }

  /**
   * @brief Modifies the first parameter such that an entry at index \f$ g,p \f$, with
   * \f$ 0 \leq g \leq num_generators \f$ and \f$ 0 \leq p \leq num_parameters \f$ is equal to \f$ f(g,p) * r(p) \f$
   * if \f$ p \leq length_r \f$ and to \f$ f(g,p) \f$ otherwise.
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
   * @param f First element of the multiplication.
   * @param r Second element of the multiplication.
   */
  template <class ValueRange, class = std::enable_if_t<RangeTraits<ValueRange>::has_begin> >
  friend Multi_parameter_generator &operator*=(Multi_parameter_generator &f, const ValueRange &r)
  {
    if (f.num_parameters() == 0 || f.is_nan()) return f;

    if constexpr (RangeTraits<ValueRange>::is_dynamic_multi_filtration){
      if (r[0].is_finite()){
        return f *= r[0];
      } else {
        if constexpr (!std::numeric_limits<T>::has_quiet_NaN)
          if (r[0].size() == 0) {
            f = nan();
            return f;
          }
        return f *= r[0][0];
      }
    } else {
      if (!f.is_finite()) {
        f.generator_.resize(r.size(), f[0]);
      }
      f._apply_operation(r, [](T &valF, const T &valR) -> bool { return _multiply(valF, valR); });
    }

    return f;
  }

  /**
   * @brief Modifies the first parameter such that an entry at index \f$ g,p \f$ is equal to \f$ f(g,p) * val \f$.
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
  friend Multi_parameter_generator &operator*=(Multi_parameter_generator &f, const T &val)
  {
    f._apply_operation(val, [](T &valF, const T &valR) -> bool { return _multiply(valF, valR); });
    return f;
  }

  // division
  /**
   * @brief Returns a filtration value such that an entry at index \f$ g,p \f$, with
   * \f$ 0 \leq g \leq num_generators \f$ and \f$ 0 \leq p \leq num_parameters \f$ is equal to \f$ f(g,p) / r(p) \f$
   * if \f$ p \leq length_r \f$ and to \f$ f(g,p) \f$ otherwise.
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
   * @param f First element of the division.
   * @param r Second element of the division.
   */
  template <class ValueRange, class = std::enable_if_t<RangeTraits<ValueRange>::has_begin> >
  friend Multi_parameter_generator operator/(Multi_parameter_generator f, const ValueRange &r)
  {
    f /= r;
    return f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ g,p \f$, with
   * \f$ 0 \leq g \leq num_generators \f$ and \f$ 0 \leq p \leq num_parameters \f$ is equal to \f$ r(p) / f(g,p) \f$
   * if \f$ p \leq length_r \f$ and to \f$ f(g,p) \f$ otherwise.
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
   * @param f Second element of the division.
   */
  template <class ValueRange,
            class = std::enable_if_t<RangeTraits<ValueRange>::has_begin &&
                                     !std::is_same_v<ValueRange, Multi_parameter_generator> > >
  friend Multi_parameter_generator operator/(const ValueRange &r, const Multi_parameter_generator &f)
  {
    if (f.num_parameters() == 0) return f;
    if (f.is_nan()) return nan();

    if constexpr (RangeTraits<ValueRange>::is_dynamic_multi_filtration){
      if (r[0].is_finite()){
        return r[0] / f;
      } else {
        if constexpr (!std::numeric_limits<T>::has_quiet_NaN) if (r[0].size() == 0) return nan();
        return r[0][0] / f;
      }
    } else {
      if (f.is_finite()) {
        Multi_parameter_generator res = f;
        res._apply_operation(r, [](T &valF, const T &valR) -> bool {
          T tmp = valF;
          valF = valR;
          return _divide(valF, tmp);
        });
        return res;
      } else {
        Multi_parameter_generator res(r.begin(), r.end());
        res._apply_operation(f[0], [](T &valRes, const T &valF) -> bool { return _divide(valRes, valF); });
        return res;
      }
    }
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ g,p \f$ is equal to \f$ f(g,p) / val \f$.
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
  friend Multi_parameter_generator operator/(Multi_parameter_generator f, const T &val)
  {
    f /= val;
    return f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ g,p \f$ is equal to \f$ val / f(g,p) \f$.
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
  friend Multi_parameter_generator operator/(const T &val, Multi_parameter_generator f)
  {
    f._apply_operation(val, [](T &valF, const T &valR) -> bool {
      T tmp = valF;
      valF = valR;
      return _divide(valF, tmp);
    });
    return f;
  }

  /**
   * @brief Modifies the first parameter such that an entry at index \f$ g,p \f$, with
   * \f$ 0 \leq g \leq num_generators \f$ and \f$ 0 \leq p \leq num_parameters \f$ is equal to \f$ f(g,p) / r(p) \f$
   * if \f$ p \leq length_r \f$ and to \f$ f(g,p) \f$ otherwise.
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
   * @param f First element of the division.
   * @param r Second element of the division.
   */
  template <class ValueRange, class = std::enable_if_t<RangeTraits<ValueRange>::has_begin> >
  friend Multi_parameter_generator &operator/=(Multi_parameter_generator &f, const ValueRange &r)
  {
    if (f.num_parameters() == 0 || f.is_nan()) return f;

    if constexpr (RangeTraits<ValueRange>::is_dynamic_multi_filtration){
      if (r[0].is_finite()){
        return f /= r[0];
      } else {
        if constexpr (!std::numeric_limits<T>::has_quiet_NaN)
          if (r[0].size() == 0) {
            f = nan();
            return f;
          }
        return f /= r[0][0];
      }
    } else {
      if (!f.is_finite()) {
        f.generator_.resize(r.size(), f[0]);
      }
      f._apply_operation(r, [](T &valF, const T &valR) -> bool { return _divide(valF, valR); });
    }

    return f;
  }

  /**
   * @brief Modifies the first parameter such that an entry at index \f$ g,p \f$ is equal to \f$ f(g,p) / val \f$.
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
  friend Multi_parameter_generator &operator/=(Multi_parameter_generator &f, const T &val)
  {
    f._apply_operation(val, [](T &valF, const T &valR) -> bool { return _divide(valF, valR); });
    return f;
  }

  // MODIFIERS

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
    if (is_nan() || is_plus_inf()) return false;
    if (x.size() == 0){
      *this = nan();
      return true;
    }

    if constexpr (RangeTraits<GeneratorRange>::is_dynamic_multi_filtration){
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
        if (!exclude_infinite_values || (valX != T_inf && valX != -T_inf)) {
          modified |= generator_[p] < valX;
          generator_[p] = valX > generator_[p] ? valX : generator_[p];
        }
        if (it + 1 != x.end()) ++it;
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
            class = std::enable_if_t<RangeTraits<GeneratorRange>::has_begin> >
  bool pull_to_greatest_common_lower_bound(const GeneratorRange &x, bool exclude_infinite_values = false)
  {
    if (is_nan() || is_minus_inf()) return false;
    if (x.size() == 0){
      *this = nan();
      return true;
    }

    if constexpr (RangeTraits<GeneratorRange>::is_dynamic_multi_filtration){
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
        if (!exclude_infinite_values || (valX != T_inf && valX != -T_inf)) {
          modified |= generator_[p] > valX;
          generator_[p] = valX < generator_[p] ? valX : generator_[p];
        }
        if (it + 1 != x.end()) ++it;
      }
  
      return modified;
    }
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
    GUDHI_CHECK(grid.size() >= generator_.size(),
                "The grid should not be smaller than the number of parameters in the filtration value.");

    if (is_nan() || generator_.empty()) return;

    if (!is_finite()) generator_.resize(grid.size(), generator_[0]);

    for (size_type p = 0; p < generator_.size(); ++p) {
      const auto &filtration = grid[p];
      auto d = std::distance(
          filtration.begin(),
          std::lower_bound(
              filtration.begin(), filtration.end(), static_cast<typename OneDimArray::value_type>(generator_[p])));
      generator_[p] = coordinate ? static_cast<T>(d) : static_cast<T>(filtration[d]);
    }
  }

  // FONCTIONNALITIES

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
  friend U compute_linear_projection(const Multi_parameter_generator &f, const std::vector<U> &x)
  {
    U projection = 0;
    std::size_t size = std::min(x.size(), f.num_parameters());
    for (std::size_t i = 0; i < size; i++) projection += x[i] * static_cast<U>(f[i]);
    return projection;
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
  friend U compute_euclidean_distance_to(const Multi_parameter_generator &f, const Multi_parameter_generator &other)
  {
    if (!f.is_finite() || !other.is_finite()) {
      if constexpr (std::numeric_limits<T>::has_quiet_NaN)
        return std::numeric_limits<T>::quiet_NaN();
      else
        return T_inf;
    }

    GUDHI_CHECK(f.num_parameters() == other.num_parameters(),
                "We cannot compute the distance between two points of different dimensions.");

    if (f.num_parameters() == 1) return f[0] - other[0];

    U out = 0;
    for (size_type p = 0; p < f.num_parameters(); ++p) {
      T v = f[p] - other[p];
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
   * @brief Computes the norm of the given filtration value.
   *
   * The filtration value is seen as a \f$ num_generators x num_parameters \f$ matrix and a standard Frobenius norm
   * is computed from it: the square root of the sum of the squares of all elements in the matrix.
   *
   * @param f Filtration value.
   * @return The norm of @p f.
   */
  template <typename U = T>
  friend U compute_squares(const Multi_parameter_generator &f)
  {
    U out = 0;
    for (size_type p = 0; p < f.num_parameters(); ++p) {
      out += f[p] * f[p];
    }
    return out;
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
  friend Multi_parameter_generator<U, Co> evaluate_coordinates_in_grid(const Multi_parameter_generator &f,
                                                                       const std::vector<std::vector<U> > &grid)
  {
    GUDHI_CHECK(grid.size() >= f.num_parameters(),
                "The size of the grid should correspond to the number of parameters in the filtration value.");

    U grid_inf = Multi_parameter_generator<U, Co>::T_inf;
    std::vector<U> outVec(f.num_parameters());

    for (size_type p = 0; p < f.num_parameters(); ++p) {
      outVec[p] = (f[p] == f.T_inf ? grid_inf : grid[p][f[p]]);
    }

    return Multi_parameter_generator<U, Co>(std::move(outVec));
  }

  // UTILITIES

  /**
   * @brief Outstream operator.
   */
  friend std::ostream &operator<<(std::ostream &stream, const Multi_parameter_generator &f)
  {
    if (f.is_plus_inf()) {
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

    const size_type num_param = f.num_parameters();

    stream << "[";
    for (size_type p = 0; p < num_param; ++p) {
      stream << f[p];
      if (p < num_param - 1) stream << ", ";
    }
    stream << "]";

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
    return sizeof(size_type) + sizeof(T) * value.generator_.size();
  }

  /**
   * @brief Infinity value of an entry of the filtration value.
   */
  constexpr static const T T_inf =
      std::numeric_limits<T>::has_infinity ? std::numeric_limits<T>::infinity() : std::numeric_limits<T>::max();

 private:
  Underlying_container generator_; /**< Container of the filtration value elements. */

  /**
   * @brief Default value of an element in the filtration value.
   */
  constexpr static T _get_default_value() { return Co ? T_inf : -T_inf; }

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
   * @brief Applies operation on the elements of the filtration value.
   */
  template <class ValueRange, class F>
  void _apply_operation(const ValueRange &range, F &&operate)
  {
    if (range.size() < generator_.size()){
      auto it = range.begin();
      for (size_type p = 0; p < generator_.size() && it != range.end(); ++p) {
        operate(generator_[p], *it);
        ++it;
      }
    } else {
      bool allPlusInf = true;
      bool allMinusInf = true;
      bool allNaN = true;

      auto it = range.begin();
      for (size_type p = 0; p < generator_.size() && it != range.end(); ++p) {
        bool isNotNaN = operate(generator_[p], *it);
        if (generator_[p] != T_inf) allPlusInf = false;
        if (generator_[p] != -T_inf) allMinusInf = false;
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
      if (operate(p, val)) allNaN = false;
    }

    if (allNaN) *this = nan();
  }

  constexpr static bool _subtract(T &v1, T v2) { return _add(v1, -v2); }

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
};

}  // namespace Gudhi::multi_filtration

namespace std {

template <typename T, bool Co>
class numeric_limits<Gudhi::multi_filtration::Multi_parameter_generator<T, Co> >
{
 public:
  using Generator = Gudhi::multi_filtration::Multi_parameter_generator<T, Co>;

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
