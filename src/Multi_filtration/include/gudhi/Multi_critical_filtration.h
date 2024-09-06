/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       David Loiseaux
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - 2024/08 Hannah Schreiber: Optimization and correction + numeric_limits + doc
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file Multi_critical_filtration.h
 * @author David Loiseaux
 * @brief Contains the @ref Gudhi::multi_filtration::Multi_critical_filtration class.
 */

#ifndef MULTI_CRITICAL_FILTRATIONS_H_
#define MULTI_CRITICAL_FILTRATIONS_H_

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include <gudhi/Debug_utils.h>
#include <gudhi/One_critical_filtration.h>

namespace Gudhi::multi_filtration {

/**
 * @class Multi_critical_filtration multi_critical_filtration.h gudhi/multi_critical_filtration.h
 * @ingroup multi_filtration
 *
 * @brief Class encoding the different generators, i.e., apparition times, of a \f$k\f$-critical
 * \f$\mathbb R^n\f$-filtration value, e.g., the filtration of a simplex, or of the algebraic generator of a module
 * presentation. The class can be used as a vector whose indices correspond each to a generator, i.e., a one-critical
 * filtration value. Then, the indices of each generator correspond to a particular parameter. 
 * E.g., \f$ f[i][p] \f$ will be  \f$ p^{\textit{th}} \f$ parameter of the \f$ i^{\textit{th}} \f$ generator 
 * of this filtration value with @ref Multi_critical_filtration \f$ f \f$.
 *
 * @details Overloads `std::numeric_limits` such that:
 * - `std::numeric_limits<Multi_critical_filtration<T,co> >::has_infinity` returns `true`,
 * - `std::numeric_limits<Multi_critical_filtration<T,co> >::infinity()` returns
 * @ref Multi_critical_filtration<T,co>::inf() "",
 * - `std::numeric_limits<Multi_critical_filtration<T,co> >::minus_infinity()` returns
 *   @ref Multi_critical_filtration<T,co>::minus_inf() "",
 * - `std::numeric_limits<Multi_critical_filtration<T,co> >::max()` throws,
 * - `std::numeric_limits<Multi_critical_filtration<T,co> >::max(g,n)` returns a @ref Multi_critical_filtration<T,co>
 * with `g` generators of `n` parameters evaluated at value `std::numeric_limits<T>::max()`,
 * - `std::numeric_limits<Multi_critical_filtration<T,co> >::quiet_NaN()` returns
 * @ref Multi_critical_filtration<T,co>::nan() "".
 *
 * Multi-critical filtrations are filtrations such that the lifetime of each object is union of positive cones in
 * \f$\mathbb R^n\f$, e.g.,
 *  - \f$ \{ x \in \mathbb R^2 : x \ge (1,2)\} \cap \{ x \in \mathbb R^2 : x \ge (2,1)\} \f$ is finitely critical,
 *    and more particularly 2-critical, while
 *  - \f$ \{ x \in \mathbb R^2 : x \ge \mathrm{epigraph}(y \mapsto e^{-y})\} \f$ is not.
 *
 * The particular case of 1-critical filtrations is handled by @ref One_critical_filtration "".
 *
 * @tparam T Arithmetic type of an entry for one parameter of a filtration value. Has to be **signed** and
 * to implement `std::isnan(T)`, `std::numeric_limits<T>::has_quiet_NaN`, `std::numeric_limits<T>::quiet_NaN()`,
 * `std::numeric_limits<T>::has_infinity`, `std::numeric_limits<T>::infinity()` and `std::numeric_limits<T>::max()`.
 * If `std::numeric_limits<T>::has_infinity` returns `false`, a call to `std::numeric_limits<T>::infinity()`
 * can simply throw. Examples are the native types `double`, `float` and `int`.
 * @tparam co If `true`, reverses the poset order, i.e., the order \f$ \le \f$  in \f$ \mathbb R^n \f$ becomes
 * \f$ \ge \f$.
 */
template <typename T, bool co = false>
class Multi_critical_filtration {
 public:
  /**
   * @brief Type of the origin of a "lifetime cone". Common with @ref One_critical_filtration "".
   */
  using Generator = One_critical_filtration<T>;
  using Generators = std::vector<Generator>;                  /**< Container type for the filtration values. */
  using iterator = typename Generators::iterator;             /**< Iterator type for the generator container. */
  using const_iterator = typename Generators::const_iterator; /**< Const iterator type for the generator container. */

  // CONSTRUCTORS

  /**
   * @brief Default constructor. The constructed value will be either at infinity if `co` is true or at minus infinity
   * if `co` is false.
   */
  Multi_critical_filtration() : multi_filtration_(_get_default_filtration_value()) {};
  /**
   * @brief Constructs a filtration value with one generator and @p n parameters.
   * All parameters will be initialized at -inf if `co` is false and at inf if `co` is true.
   *
   * @warning The generator `{-inf, -inf, ...}`/`{inf, inf, ...}` with \f$ n > 1 \f$ entries is not considered as
   * "(minus) infinity" (the resp. methods @ref is_minus_inf() and @ref is_plus_inf() "", as well as the ones of the
   * generator, will not return true). The `-inf/inf` are just meant as placeholders, at least one entry should be
   * modified by the user.
   * Otherwise, either use the static methods @ref minus_inf() or @ref inf(), or set @p n to 1 instead.
   *
   * @param n Number of parameters.
   */
  Multi_critical_filtration(int n) : multi_filtration_(1, Generator(n, _get_default_value())) {};
  /**
   * @brief Constructs a filtration value with one generator and @p n parameters. All parameters will be initialized
   * with @p value.
   *
   * @warning If @p value is `inf`, `-inf`, or `NaN`, the generator `{value, value, ...}` with \f$ n > 1 \f$ entries
   * is not wrong but will not be considered as respectively "infinity", "minus infinity" or "NaN" (the corresponding
   * methods @ref is_plus_inf(), @ref is_minus_inf() and @ref is_nan() will return false). For this purpose, please use
   * the static methods @ref inf(), @ref minus_inf() and @ref nan() instead.
   *
   * @param n Number of parameters.
   * @param value Value which will be used for each entry.
   */
  Multi_critical_filtration(int n, T value) : multi_filtration_(1, Generator(n, value)) {};
  /**
   * @brief Constructs a filtration value with one generator which will be initialzed by the given initializer list.
   *
   * @param init Initializer list with values for each parameter.
   */
  Multi_critical_filtration(std::initializer_list<T> init) : multi_filtration_(1, Generator{init}) {};
  /**
   * @brief Constructs a filtration value with one generator which will be initialzed by the given vector.
   *
   * @param v Vector with values for each parameter.
   */
  Multi_critical_filtration(const std::vector<T> &v) : multi_filtration_(1, Generator{v}) {};
  /**
   * @brief Constructs a filtration value with one generator to which the given vector is moved to.
   *
   * @param v Vector with values for each parameter.
   */
  Multi_critical_filtration(std::vector<T> &&v) : multi_filtration_(1, Generator{std::move(v)}) {};
  /**
   * @brief Constructs filtration value with as many generators than elements in the given vector and initialize
   * them with them.
   * If the vector is empty, then the filtration value is either initialized at infinity if `co` is true or at
   * minus infinity if `co` is false.
   * @pre All generators in the vector have to have the same number of parameters, i.e., size.
   * Furthermore, the generators have to be a minimal generating set.
   *
   * @warning If the set of generators is not minimal or not sorted, the behaviour of most methods is undefined.
   * It is possible to call @ref simplify() after construction if there is a doubt to ensure this property.
   *
   * @param v Vector of generators.
   */
  Multi_critical_filtration(const std::vector<Generator> &v)
      : multi_filtration_(v.empty() ? _get_default_filtration_value() : v) {};
  /**
   * @brief Constructs filtration value with as many generators than elements in the given vector and moves those
   * elements to initialize the generators.
   * If the vector is empty, then the filtration value is either initialized at infinity if `co` is true or at
   * minus infinity if `co` is false.
   * @pre All generators in the vector have to have the same number of parameters, i.e., size.
   * Furthermore, the generators have to be a minimal generating set.
   *
   * @warning If the set of generators is not minimal or not sorted, the behaviour of most methods is undefined.
   * It is possible to call @ref simplify() after construction if there is a doubt to ensure this property.
   *
   * @param v Vector of generators.
   */
  Multi_critical_filtration(std::vector<Generator> &&v)
      : multi_filtration_(v.empty() ? _get_default_filtration_value() : std::move(v)) {};
  /**
   * @brief Constructs a filtration value with one generator initialzed by the range given by the begin and end
   * iterators.
   *
   * @param it_begin Start of the range.
   * @param it_end End of the range.
   */
  Multi_critical_filtration(typename std::vector<T>::iterator it_begin, typename std::vector<T>::iterator it_end)
      : multi_filtration_(Generators(1, {it_begin, it_end})) {};
  /**
   * @brief Constructs a filtration value with one generator initialzed by the range given by the begin and end
   * const iterators.
   *
   * @param it_begin Start of the range.
   * @param it_end End of the range.
   */
  Multi_critical_filtration(typename std::vector<T>::const_iterator it_begin,
                            typename std::vector<T>::const_iterator it_end)
      : multi_filtration_(Generator(1, {it_begin, it_end})) {};

  // VECTOR-LIKE

  using value_type = T; /**< Entry type. */

  /**
   * @brief Standard operator[].
   */
  Generator &operator[](std::size_t i) { return multi_filtration_[i]; }
  /**
   * @brief Standard operator[] const.
   */
  const Generator &operator[](std::size_t i) const { return multi_filtration_[i]; }

  /**
   * @brief Returns begin iterator of the generator range.
   *
   * @warning If the generator is modified and the new set of generators is not minimal or not sorted, the behaviour
   * of most methods is undefined. It is possible to call @ref simplify() after construction if there is a doubt to
   * ensure this property.
   */
  iterator begin() { return multi_filtration_.begin(); }

  /**
   * @brief Returns end iterator of the generator range.
   *
   * @warning If the generator is modified and the new set of generators is not minimal or not sorted, the behaviour
   * of most methods is undefined. It is possible to call @ref simplify() after construction if there is a doubt to
   * ensure this property.
   */
  iterator end() { return multi_filtration_.end(); }

  /**
   * @brief Returns begin const iterator of the generator range.
   */
  const_iterator begin() const { return multi_filtration_.begin(); }

  /**
   * @brief Returns end const iterator of the generator range.
   */
  const_iterator end() const { return multi_filtration_.end(); }

  /**
   * @brief Reserves space for the given number of generators in the underlying container.
   *
   * @param n Number of generators.
   */
  void reserve(std::size_t n) { multi_filtration_.reserve(n); }

  // CONVERTERS

  /**
   * @brief Casts the object into the type of a generator.
   * @pre The filtration value is 1-critical. If there are more than one generator, only the first will be preserved
   * and if there is no generator, the method will segfault.
   */
  operator Generator() {
    GUDHI_CHECK(num_generators() == 1, "Casting a " + std::to_string(num_generators()) +
                                           "-critical filtration value into an 1-critical filtration value.");
    return multi_filtration_[0];
  }

  // like numpy
  /**
   * @brief Returns a copy with entries casted into the type given as template parameter.
   *
   * @tparam U New type for the entries.
   * @return Copy with new entry type.
   */
  template <typename U>
  Multi_critical_filtration<U> as_type() const {
    std::vector<One_critical_filtration<U>> out(num_generators());
    for (std::size_t i = 0u; i < num_generators(); ++i) {
      out[i] = multi_filtration_[i].template as_type<U>();
    }
    return Multi_critical_filtration<U>(std::move(out));
  }

  // ACCESS

  /**
   * @brief Returns a reference to the underlying container storing the generators.
   *
   * @warning If a generator is modified and the new set of generators is not minimal or not sorted, the behaviour
   * of most methods is undefined. It is possible to call @ref simplify() after construction if there is a doubt to
   * ensure this property.
   */
  const Generators &get_underlying_container() const { return multi_filtration_; }

  /**
   * @brief Returns the number of parameters.
   */
  std::size_t num_parameters() const { return multi_filtration_[0].num_parameters(); }

  /**
   * @brief Returns the number of generators.
   */
  std::size_t num_generators() const { return multi_filtration_.size(); }

  /**
   * @brief Returns a filtration value for which @ref is_plus_inf() returns `true`.
   *
   * @return Infinity.
   */
  constexpr static Multi_critical_filtration inf() { return Multi_critical_filtration(Generator::inf()); }

  /**
   * @brief Returns a filtration value for which @ref is_minus_inf() returns `true`.
   *
   * @return Minus infinity.
   */
  constexpr static Multi_critical_filtration minus_inf() { return Multi_critical_filtration(Generator::minus_inf()); }

  /**
   * @brief Returns a filtration value for which @ref is_nan() returns `true`.
   *
   * @return NaN.
   */
  constexpr static Multi_critical_filtration nan() { return Multi_critical_filtration(Generator::nan()); }

  // DESCRIPTORS

  // TODO: Accept {{-inf, -inf, ...},...} / {{inf, inf, ...},...} / {{NaN, NaN, ...},...} as resp. -inf / inf / NaN.

  /**
   * @brief Returns `true` if and only if the filtration value is considered as infinity.
   */
  bool is_plus_inf() const { return multi_filtration_.size() == 1 && multi_filtration_[0].is_plus_inf(); }

  /**
   * @brief Returns `true` if and only if the filtration value is considered as minus infinity.
   */
  bool is_minus_inf() const { return multi_filtration_.size() == 1 && multi_filtration_[0].is_minus_inf(); }

  /**
   * @brief Returns `true` if and only if the filtration value is considered as NaN.
   */
  bool is_nan() const { return multi_filtration_.size() == 1 && multi_filtration_[0].is_nan(); }

  /**
   * @brief Returns `true` if and only if the filtration value is non-empty and is not considered as infinity,
   * minus infinity or NaN.
   */
  bool is_finite() const {
    if (multi_filtration_.size() > 1) return true;
    return multi_filtration_[0].is_finite();
  }

  // COMPARAISON OPERATORS

  // TODO : this costs a lot... optimize / cheat in some way for python ?

  /**
   * @brief Returns `true` if and only if the positive cones generated by @p b are strictly contained in the
   * positive cones generated by @p a.
   * If @p a and @p b are both not infinite or NaN, they have to have the same number of parameters.
   *
   * Note that not all filtration values are comparable. That is, \f$ a > b \f$ and \f$ b > a \f$ returning both false
   * does **not** imply \f$ a == b \f$.
   */
  friend bool operator<(const Multi_critical_filtration &a, const Multi_critical_filtration &b) {
    for (std::size_t i = 0u; i < b.multi_filtration_.size(); ++i) {
      // for each generator in b, verify if it is strictly in the cone of at least one generator of a
      bool isContained = false;
      for (std::size_t j = 0u; j < a.multi_filtration_.size() && !isContained; ++j) {
        // lexicographical order, so if a[j][0] dom b[j][0], than a[j'] can never strictly contain b[i] for all j' > j.
        if (_first_dominates(a.multi_filtration_[j], b.multi_filtration_[i])) return false;
        isContained = _strictly_contains(a.multi_filtration_[j], b.multi_filtration_[i]);
      }
      if (!isContained) return false;
    }
    return true;
  }

  /**
   * @brief Returns `true` if and only if the positive cones generated by @p a are strictly contained in the
   * positive cones generated by @p b.
   * If @p a and @p b are both not infinite or NaN, they have to have the same number of parameters.
   *
   * Note that not all filtration values are comparable. That is, \f$ a > b \f$ and \f$ b > a \f$ returning both false
   * does **not** imply \f$ a == b \f$.
   */
  friend bool operator>(const Multi_critical_filtration &a, const Multi_critical_filtration &b) { return b < a; }

  /**
   * @brief Returns `true` if and only if the positive cones generated by @p b are contained in or are (partially)
   * equal to the positive cones generated by @p a.
   * If @p a and @p b are both not infinite or NaN, they have to have the same number of parameters.
   *
   * Note that not all filtration values are comparable. That is, \f$ a \le b \f$ and \f$ b \le a \f$ can both return
   * `false`.
   */
  friend bool operator<=(const Multi_critical_filtration &a, const Multi_critical_filtration &b) {
    // check if this curves is below other's curve
    //  ie for each guy in this, check if there is a guy in other that dominates him
    for (std::size_t i = 0u; i < b.multi_filtration_.size(); ++i) {
      // for each generator in b, verify if it is in the cone of at least one generator of a
      bool isContained = false;
      for (std::size_t j = 0u; j < a.multi_filtration_.size() && !isContained; ++j) {
        // lexicographical order, so if a[j][0] strictly dom b[j][0], than a[j'] can never contain b[i] for all j' > j.
        if (_first_strictly_dominates(a.multi_filtration_[j], b.multi_filtration_[i])) return false;
        isContained = _contains(a.multi_filtration_[j], b.multi_filtration_[i]);
      }
      if (!isContained) return false;
    }
    return true;
  }

  /**
   * @brief Returns `true` if and only if the positive cones generated by @p a are contained in or are (partially)
   * equal to the positive cones generated by @p b.
   * If @p a and @p b are both not infinite or NaN, they have to have the same number of parameters.
   *
   * Note that not all filtration values are comparable. That is, \f$ a \ge b \f$ and \f$ b \ge a \f$ can both return
   * `false`.
   */
  friend bool operator>=(const Multi_critical_filtration &a, const Multi_critical_filtration &b) { return b <= a; }

  /**
   * @brief Returns `true` if and only if for each \f$ i \f$, \f$ a[i] \f$ is equal to \f$ b[i] \f$.
   */
  friend bool operator==(const Multi_critical_filtration &a, const Multi_critical_filtration &b) {
    // assumes lexicographical order for both
    return a.multi_filtration_ == b.multi_filtration_;
  }

  /**
   * @brief Returns `true` if and only if \f$ a == b \f$ returns `false`.
   */
  friend bool operator!=(const Multi_critical_filtration &a, const Multi_critical_filtration &b) { return !(a == b); }

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
  void set_num_generators(std::size_t n) {
    if (n == 0) return;
    multi_filtration_.resize(n);
  }

  /**
   * @brief Sets all generators to the least common upper bound between the current generator value and the given value.
   *
   * More formally, it pushes the current filtration value to the cone \f$ \{ y \in \mathbb R^n : y \ge x \} \f$
   * originating in \f$ x \f$. The resulting values corresponds to the generators of the intersection of this cone
   * with the union of positive cones generated by the old generators.
   *
   * @param x The target filtration value towards which to push.
   */
  void push_to_least_common_upper_bound(const Generator &x) {
    if (this->is_plus_inf() || this->is_nan() || x.is_nan() || x.is_minus_inf()) return;

    GUDHI_CHECK(x.is_plus_inf() || x.num_parameters() == multi_filtration_[0].num_parameters() || !is_finite(),
                "Pushing to a filtration value with different number of parameters.");

    if (x.is_plus_inf() || this->is_minus_inf()) {
      multi_filtration_ = {x};
      return;
    }
    for (auto &g : *this) {
      g.push_to_least_common_upper_bound(x);
    }

    simplify();
  }

  /**
   * @brief Sets all generators to the greatest common lower bound between the current generator value and the given
   * value.
   *
   * More formally, it pulls the current filtration value to the cone \f$ \{ y \in \mathbb R^n : y \le x \} \f$
   * originating in \f$ x \f$. The resulting values corresponds to the generators of the intersection of this cone
   * with the union of negative cones generated by the old generators.
   *
   * @param x The target filtration value towards which to pull.
   */
  void pull_to_greatest_common_lower_bound(const Generator &x) {
    if (x.is_plus_inf() || this->is_nan() || x.is_nan() || this->is_minus_inf()) return;

    GUDHI_CHECK(x.is_minus_inf() || x.num_parameters() == multi_filtration_[0].num_parameters() || !is_finite(),
                "Pulling to a filtration value with different number of parameters.");

    if (this->is_plus_inf() || x.is_minus_inf()) {
      multi_filtration_ = {x};
      return;
    }
    for (auto &g : *this) {
      g.pull_to_greatest_common_lower_bound(x);
    }

    simplify();
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
  bool add_generator(const Generator &x) {
    GUDHI_CHECK(x.num_parameters() == multi_filtration_[0].num_parameters() || !is_finite() || !x.is_finite(),
                "Cannot add a generator with different number of parameters.");

    std::size_t end = multi_filtration_.size();

    if (_generator_can_be_added(x, 0, end)) {
      multi_filtration_.resize(end);
      multi_filtration_.push_back(x);
      std::sort(multi_filtration_.begin(), multi_filtration_.end(), Is_strictly_smaller_lexicographically());
      return true;
    }

    return false;
  }

  /**
   * @brief Adds the given generator to the filtration value without any verifications or simplifications.
   *
   * @warning If the resulting set of generators is not minimal after modification, some methods will have an
   * undefined behaviour. Be sure to call @ref simplify() before using them.
   *
   * @param x
   */
  void add_guaranteed_generator(const Generator &x) { multi_filtration_.push_back(x); }

  /*
   * Same as `compute_coordinates_in_grid`, but does the operation in-place.
   */

  /**
   * @brief Projects the filtration value into the given grid. If @p coordinate is false, the entries are set to
   * the nearest upper bound value with the same parameter in the grid and the new generators are simplified and
   * ordered. Otherwise, the entries are set to the indices of those nearest upper bound values. In this case,
   * no simplification or sort are done, such that the new coordinates have a one by one correspondence with the
   * positions of the old generators.
   * The grid has to be represented as a vector of ordered ranges of values convertible into `T`. An index
   * \f$ i \f$ of the vector corresponds to the same parameter as the index \f$ i \f$ in a generator.
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
    GUDHI_CHECK(grid.size() >= num_parameters(),
                "The grid should not be smaller than the number of parameters in the filtration value.");

    for (auto &x : multi_filtration_) {
      x.project_onto_grid(grid, coordinate);
    }

    if (!coordinate) simplify();
  }

  /**
   * @brief Removes all empty generators from the filtration value. If @p include_infinities is true, it also
   * removes the generators at infinity or minus infinity.
   * If the set of generators is empty after removals, it is set to minus infinity if `co` is false or to infinity
   * if `co` is true.
   *
   * @param include_infinities If true, removes also infinity values.
   */
  void remove_empty_generators(bool include_infinities = false) {
    multi_filtration_.erase(std::remove_if(multi_filtration_.begin(), multi_filtration_.end(),
                                           [include_infinities](const Generator &a) {
                                             return a.empty() ||
                                                    ((include_infinities) && (a.is_plus_inf() || a.is_minus_inf()));
                                           }),
                            multi_filtration_.end());
    std::sort(multi_filtration_.begin(), multi_filtration_.end(), Is_strictly_smaller_lexicographically());
    if (multi_filtration_.empty()) multi_filtration_.push_back(Generator{_get_default_value()});
  }

  /**
   * @brief Simplifies the current set of generators such that it becomes minimal. Also orders it in increasing
   * lexicographical order. Only necessary if generators were added "by hand" without verification either trough the
   * constructor or with @ref add_guaranteed_generator "", etc.
   */
  void simplify() {
    std::size_t end = 0;

    for (std::size_t curr = 0; curr < multi_filtration_.size(); ++curr) {
      if (_generator_can_be_added(multi_filtration_[curr], 0, end)) {
        std::swap(multi_filtration_[end], multi_filtration_[curr]);
        ++end;
      }
    }

    multi_filtration_.resize(end);
    std::sort(multi_filtration_.begin(), multi_filtration_.end(), Is_strictly_smaller_lexicographically());
  }

  // FONCTIONNALITIES

  /**
   * @brief Returns a generator with the minimal values of all parameters in any generator of the given filtration
   * value. That is, the greatest lower bound of all generators.
   */
  friend Generator factorize_below(const Multi_critical_filtration &f) {
    if (f.num_generators() == 0) return Generator();
    Generator result(f.num_parameters(), Generator::T_inf);
    for (const auto &g : f) {
      if (g.is_nan() || g.is_minus_inf()) return g;
      if (g.is_plus_inf()) continue;
      for (std::size_t i = 0; i < f.num_parameters(); ++i) {
        result[i] = std::min(result[i], g[i]);
      }
    }
    return result;
  }

  /**
   * @brief Returns a generator with the maximal values of all parameters in any generator of the given filtration
   * value. That is, the least upper bound of all generators.
   */
  friend Generator factorize_above(const Multi_critical_filtration &f) {
    if (f.num_generators() == 0) return Generator();
    Generator result(f.num_parameters(), -Generator::T_inf);
    for (auto &g : f) {
      if (g.is_nan() || g.is_plus_inf()) return g;
      if (g.is_minus_inf()) continue;
      for (std::size_t i = 0; i < g.num_parameters(); ++i) {
        result[i] = std::max(result[i], g[i]);
      }
    }
    return result;
  }

  /**
   * @brief Computes the smallest (resp. the greatest if `co` is true) scalar product of the all generators with the
   * given vector.
   *
   * @tparam U Arithmetic type of the result. Default value: `T`.
   * @param f Filtration value.
   * @param x Vector of coefficients.
   * @return Scalar product of @p f with @p x.
   */
  template <typename U = T>
  friend U compute_linear_projection(const Multi_critical_filtration &f, const std::vector<U> &x) {
    if constexpr (co) {
      U projection = std::numeric_limits<U>::lowest();
      for (const auto &y : f) {
        projection = std::max(projection, compute_linear_projection(y, x));
      }
      return projection;
    } else {
      U projection = std::numeric_limits<U>::max();
      for (const auto &y : f) {
        projection = std::min(projection, compute_linear_projection(y, x));
      }
      return projection;
    }
  }

  /**
   * @brief Computes the coordinates in the given grid, corresponding to the nearest upper bounds of the entries
   * in the given filtration value.
   * The grid has to be represented as a vector of vectors of ordered values convertible into `out_type`. An index
   * \f$ i \f$ of the vector corresponds to the same parameter as the index \f$ i \f$ in a generator.
   * The inner vectors correspond to the possible values of the parameters, ordered by increasing value,
   * forming therefore all together a 2D grid.
   *
   * @tparam out_type Signed arithmetic type. Default value: std::int32_t.
   * @tparam U Type which is convertible into `out_type`.
   * @param f Filtration value to project.
   * @param grid Vector of vectors to project into.
   * @return Filtration value \f$ out \f$ whose entry correspond to the indices of the projected values. That is,
   * the projection of \f$ f[i] \f$ is \f$ grid[i][out[i]] \f$ before simplification (if two generators were
   * projected to the same point, the doubles are removed in the output).
   */
  template <typename out_type = std::int32_t, typename U = T>
  friend Multi_critical_filtration<out_type> compute_coordinates_in_grid(Multi_critical_filtration f,
                                                                         const std::vector<std::vector<U>> &grid) {
    // TODO: by replicating the code of the 1-critical "project_onto_grid", this could be done with just one copy
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
   * @return Filtration value \f$ out \f$ whose entry correspond to \f$ out[i] = grid[i][f[i]] \f$ before
   * simplification (the output is simplified).
   */
  template <typename U>
  friend Multi_critical_filtration<U> evaluate_coordinates_in_grid(const Multi_critical_filtration &f,
                                                                   const std::vector<std::vector<U>> &grid) {
    Multi_critical_filtration<U> out;
    out.set_num_generators(f.num_generators());
    for (std::size_t i = 0; i < f.num_generators(); ++i) {
      out[i] = evaluate_coordinates_in_grid(f[i], grid);
    }
    out.simplify();
    return out;
  }

  // UTILITIES

  /**
   * @brief Outstream operator.
   */
  friend std::ostream &operator<<(std::ostream &stream, const Multi_critical_filtration &f) {
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
    stream << "(k=" << f.multi_filtration_.size() << ")[";
    for (const auto &val : f) {
      stream << val << "; ";
    }
    if (f.multi_filtration_.size() > 0) {
      stream << "\b"
             << "\b";
    }
    stream << "]";
    return stream;
  }

 public:
  /**
   * @brief Indicates if the class manages multi-critical filtration values.
   */
  constexpr static const bool is_multi_critical = true;

 private:
  Generators multi_filtration_; /**< Container for generators. */

  struct Is_strictly_smaller_lexicographically {
    //assumes both generators have the same length if not infinite/nan.
    bool operator()(const Generator &g1, const Generator &g2) {
      if (g1.is_nan() || g2.is_nan()) return !g1.is_nan();
      if (g1.is_plus_inf()) return false;
      if (g2.is_plus_inf()) return true;
      if (g1.is_minus_inf()) return false;
      if (g2.is_minus_inf()) return true;

      // g1 and g2 have to finite and of the same size
      for (std::size_t i = 0; i < g1.size(); ++i) {
        if (g1[i] != g2[i]) return g1[i] < g2[i];
      }
      return false;
    }
  };

  constexpr static T _get_default_value() { return co ? Generator::T_inf : -Generator::T_inf; }

  constexpr static Generators _get_default_filtration_value() { return Generators{Generator{_get_default_value()}}; }

  /**
   * @brief Verifies if @p b is strictly contained in the positive cone originating in `a`.
   */
  static bool _strictly_contains(const Generator &a, const Generator &b) {
    if constexpr (co)
      return a > b;
    else {
      return a < b;
    }
  }
  /**
   * @brief Verifies if @p b is contained in the positive cone originating in `a`.
   */
  static bool _contains(const Generator &a, const Generator &b) {
    if constexpr (co)
      return a >= b;
    else {
      return a <= b;
    }
  }
  
  static bool _first_strictly_dominates(const Generator& a, const Generator& b){
    if constexpr (co){
      return !a.empty() && !b.empty() && a[0] < b[0];
    } else {
      return !a.empty() && !b.empty() && a[0] > b[0];
    }
  }

  static bool _first_dominates(const Generator& a, const Generator& b){
    if constexpr (co){
      return !a.empty() && !b.empty() && a[0] <= b[0];
    } else {
      return !a.empty() && !b.empty() && a[0] >= b[0];
    }
  }

  enum class Rel { EQUAL, DOMINATES, IS_DOMINATED, NONE };

  static Rel _get_domination_relation(const Generator &a, const Generator &b) {
    if (a.is_nan() || b.is_nan()) return Rel::NONE;

    GUDHI_CHECK(a.size() == b.size(),
                "Two generators in the same k-critical value have to have the same numbers of parameters.");

    bool equal = true;
    bool allGreater = true;
    bool allSmaller = true;
    for (unsigned int i = 0; i < a.size(); ++i) {
      if (a[i] < b[i]) {
        if (!allSmaller) return Rel::NONE;
        equal = false;
        allGreater = false;
      } else if (a[i] > b[i]) {
        if (!allGreater) return Rel::NONE;
        equal = false;
        allSmaller = false;
      }
    }
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
  // modifies multi_filtration_ only if true is returned.
  bool _generator_can_be_added(const Generator &x, std::size_t curr, std::size_t &end) {
    if (x.empty() || x.is_nan() || (x.is_plus_inf() && end - curr != 0)) return false;

    if (x.is_minus_inf()) {
      if (end - curr == 1 && multi_filtration_[curr].is_minus_inf()) return false;
      // assumes that everything between curr and end is already simplified
      // so, if end - curr != 1, there can be no minus_inf anymore.
      end = curr;
      return true;
    }

    while (curr != end) {
      Rel res = _get_domination_relation(multi_filtration_[curr], x);
      if (res == Rel::IS_DOMINATED || res == Rel::EQUAL) return false;  // x dominates or is equal
      if (res == Rel::DOMINATES) {                                      // x is dominated
        --end;
        std::swap(multi_filtration_[curr], multi_filtration_[end]);
      } else {                                                          // no relation
        ++curr;
      }
    }
    return true;
  }
};

}  // namespace Gudhi::multi_filtration

namespace std {

template <typename T>
class numeric_limits<Gudhi::multi_filtration::Multi_critical_filtration<T>> {
 public:
  static constexpr bool has_infinity = true;

  static constexpr Gudhi::multi_filtration::Multi_critical_filtration<T> infinity() noexcept {
    return Gudhi::multi_filtration::Multi_critical_filtration<T>::inf();
  };

  // non-standard
  static constexpr Gudhi::multi_filtration::Multi_critical_filtration<T> minus_infinity() noexcept {
    return Gudhi::multi_filtration::Multi_critical_filtration<T>::minus_inf();
  };

  static constexpr Gudhi::multi_filtration::Multi_critical_filtration<T> max() noexcept(false) {
    throw std::logic_error(
        "The maximal value cannot be represented with no finite numbers of generators."
        "Use `max(number_of_generators, number_of_parameters)` instead");
  };

  // non-standard, so I don't want to define default values.
  static constexpr Gudhi::multi_filtration::Multi_critical_filtration<T> max(unsigned int g, unsigned int n) noexcept {
    std::vector<typename Gudhi::multi_filtration::Multi_critical_filtration<T>::Generator> v(
        g, std::vector<T>(n, std::numeric_limits<T>::max()));
    return Gudhi::multi_filtration::Multi_critical_filtration<T>(std::move(v));
  };

  static constexpr Gudhi::multi_filtration::Multi_critical_filtration<T> quiet_NaN() noexcept {
    return Gudhi::multi_filtration::Multi_critical_filtration<T>::nan();
  };
};

}  // namespace std

#endif  // MULTI_CRITICAL_FILTRATIONS_H_
