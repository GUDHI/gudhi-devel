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
 * @file Degree_rips_bifiltration.h
 * @author Hannah Schreiber, David Loiseaux
 * @brief Contains the @ref Gudhi::multi_filtration::Degree_rips_bifiltration class.
 */

#ifndef MF_DEGREE_RIPS_BIFILTRATION_H_
#define MF_DEGREE_RIPS_BIFILTRATION_H_

#include <algorithm>    //std::lower_bound
#include <cmath>        //std::isnan, std::min
#include <cstddef>      //std::size_t
#include <cstdint>      //std::int32_t
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
#include <gudhi/Simplex_tree/filtration_value_utils.h>
#include <gudhi/Multi_parameter_filtration.h>
#include <gudhi/Dynamic_multi_parameter_filtration.h>

namespace Gudhi::multi_filtration {

/**
 * @class Degree_rips_bifiltration Degree_rips_bifiltration.h gudhi/Degree_rips_bifiltration.h
 * @ingroup multi_filtration
 *
 * @brief Class encoding the different generators, i.e., apparition times, of a \f$ k \f$-critical
 * \f$\mathbb R^2\f$-filtration value in a degree-Rips filtration. That is, a \f$ k \f$-critical filtration value is
 * always of the form: \f$ [(v_0,0), (v_1,1), (v_2,2), ..., (v_{k-1},k-1)] \f$, where all pairs \f$ (v_i,i) \f$
 * represent a generator with two parameters: the radius and the degree. More precisely: let \f$ d \f$ be the max
 * degree of a vertex in the complex. A vertex will be \f$ k \f$-critical if it has degree \f$ d - k + 1 \f$ and an
 * edge is \f$ k \f$-critical if one of its end vertices is \f$ k \f$-critical and the other one \f$ j \f$-critical,
 * \f$ j \geq k \f$. The first parameter is the more standard radius parameter of a Rips filtration.
 * Note that the set of generators does not have to be minimal (contrary to @ref Multi_parameter_filtration e.g.),
 * neither ordered lexicographically.
 * Implements the concept @ref FiltrationValue of the @ref Gudhi::Simplex_tree and the concept
 * @ref Gudhi::multi_persistence::MultiFiltrationValue.
 *
 * @details Overloads `std::numeric_limits` such that:
 * - `std::numeric_limits<Degree_rips_bifiltration>::has_infinity` returns `true` if and only if `Co` is false,
 * - `std::numeric_limits<Degree_rips_bifiltration>::has_quiet_NaN` returns `true`,
 * - `std::numeric_limits<Degree_rips_bifiltration>::infinity()` returns
 * @ref Degree_rips_bifiltration::inf() "",
 * - `std::numeric_limits<Degree_rips_bifiltration>::minus_infinity()` returns
 * @ref Degree_rips_bifiltration::minus_inf() "",
 * - `std::numeric_limits<Degree_rips_bifiltration>::max(num_param)` throws if `Co` is true and otherwise returns a
 * @ref Degree_rips_bifiltration with 1 generators with first parameter 0 and second parameter
 *`std::numeric_limits<T>::max()`,
 * - `std::numeric_limits<Degree_rips_bifiltration>::quiet_NaN()` returns @ref Degree_rips_bifiltration::nan().
 *
 * @tparam T Arithmetic type of an entry of the second parameter of a filtration value. Has to be **signed** and
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
class Degree_rips_bifiltration
{
 public:
  using Underlying_container = std::vector<T>; /**< Underlying container for values. */

  // CONSTRUCTORS

  /**
   * @brief Default constructor. Builds filtration value with one generator `(val, 0)`.
   * If Co is false, `val` is -inf, if Co is true, `val` is at +inf.
   *
   * @param number_of_parameters Ignored, the number of parameters is always 2. For interface purposes only.
   */
  Degree_rips_bifiltration([[maybe_unused]] int number_of_parameters = 2) : generators_(1, _get_default_value()) {}

  explicit Degree_rips_bifiltration(Gudhi::simplex_tree::empty_filtration_value_t /*e*/) : generators_(0) {}

  /**
   * @brief Builds a filtration value with one generator `(value, 0)`.
   *
   * @param number_of_parameters Ignored, the number of parameters is always 2. For interface purposes only.
   * @param value Initialization value for the second parameter.
   */
  Degree_rips_bifiltration([[maybe_unused]] int number_of_parameters, T value) : generators_(1, value) {}

  /**
   * @brief Builds filtration value with one generator `(val, i)`, where `val` and `i` are the two first elements
   * of the given range. Note that `i` has to be 0.
   *
   * @tparam ValueRange Range of types convertible to `T`. Should have a begin() method.
   * @param range Values of the generator. The range has to have at least two elements.
   */
  template <class ValueRange = std::initializer_list<T>, class = std::enable_if_t<RangeTraits<ValueRange>::has_begin> >
  Degree_rips_bifiltration(const ValueRange &range) : generators_(1, *(range.begin()))
  {
    GUDHI_CHECK(*(range.begin() + 1) == 0, std::invalid_argument("Second value of the range has to be 0"));
  }

  /**
   * @brief Builds filtration value with one generator `(val, i)`, where `val` and `i` are the two first elements
   * of the given range. Note that `i` has to be 0.
   *
   * @tparam Iterator Iterator type that has to satisfy the requirements of standard LegacyInputIterator and
   * dereferenced elements have to be convertible to `T`.
   * @param it_begin Iterator pointing to the start of the range.
   * @param it_end Iterator pointing to the end of the range.
   */
  template <class Iterator, class = std::enable_if_t<!std::is_arithmetic_v<Iterator> > >
  Degree_rips_bifiltration(Iterator it_begin, [[maybe_unused]] Iterator it_end) : generators_(1, *it_begin)
  {
    GUDHI_CHECK(*(it_begin + 1) == 0, std::invalid_argument("Second value of the range has to be 0"));
  }

  /**
   * @brief Builds a filtration value with given values from the given range. The two first elements of the range have
   * to correspond to the first generator, the two next elements to the second generator and so on... So the length of
   * the range has to be a multiple of 2 and the number of generators will be \f$ k = length / 2 \f$. Note that starting
   * from the second element, every second element has to represent the continuous sequence from 0 to \f$ k \f$.
   * The range is represented by two iterators.
   *
   * @tparam Iterator Iterator type that has to satisfy the requirements of standard LegacyForwardIterator and
   * dereferenced elements have to be convertible to `T`.
   * @param it_begin Iterator pointing to the start of the range.
   * @param it_end Iterator pointing to the end of the range.
   * @param number_of_parameters Ignored, the number of parameters is always 2. For interface purposes only.
   */
  template <class Iterator, class = std::enable_if_t<!std::is_arithmetic_v<Iterator> > >
  Degree_rips_bifiltration(Iterator it_begin, Iterator it_end, [[maybe_unused]] int number_of_parameters)
      : generators_()
  {
    size_type num_gen = std::distance(it_begin, it_end) / 2;
    if constexpr (Ensure1Criticality) {
      if (num_gen > 1) throw std::logic_error("Multiparameter filtration value is not 1-critical.");
    }
    generators_.resize(num_gen);
    Iterator it = it_begin;
    for (size_type i = 0; i < num_gen; ++i) {
      generators_[i] = *it;
      ++it;
      GUDHI_CHECK(
          static_cast<size_type>(*it) == i,
          std::invalid_argument(
              "Every second value of the range has to correspond to a contiguous sequence of integers starting at 0."));
      ++it;
    }
  }

  /**
   * @brief Builds filtration value with given values from the given range. The range only represent the values of the
   * first parameter. So `(generators[0], 0)` is the first generator, `(generators[1], 1)` is the second generator and
   * so on... The range is represented by @ref Degree_rips_bifiltration::Underlying_container "" and copied into the
   * underlying container of the class.
   *
   * @param generators Values for the second parameter.
   * @param number_of_parameters Ignored, the number of parameters is always 2. For interface purposes only.
   */
  Degree_rips_bifiltration(const Underlying_container &generators, [[maybe_unused]] int number_of_parameters)
      : generators_(generators)
  {
    if constexpr (Ensure1Criticality) {
      if (generators_.size() > 1) throw std::logic_error("Multiparameter filtration value is not 1-critical.");
    }
  }

  /**
   * @brief Builds filtration value with given values from the given range. The range only represent the values of the
   * first parameter. So `(generators[0], 0)` is the first generator, `(generators[1], 1)` is the second generator and
   * so on... The range is represented by @ref Degree_rips_bifiltration::Underlying_container "" and **moved** into
   * the underlying container of the class.
   *
   * @param generators Values to move.
   * @param number_of_parameters Ignored, the number of parameters is always 2. For interface purposes only.
   */
  Degree_rips_bifiltration(Underlying_container &&generators, [[maybe_unused]] int number_of_parameters)
      : generators_(std::move(generators))
  {
    if constexpr (Ensure1Criticality) {
      if (generators_.size() > 1) throw std::logic_error("Multiparameter filtration value is not 1-critical.");
    }
  }

  // cannot use = default as it triggers "dummy_g_ may be used uninitialized" compiler warning for nothing
  /**
   * @brief Copy constructor.
   */
  Degree_rips_bifiltration(const Degree_rips_bifiltration &other) : generators_(other.generators_) {}

  // cannot use = default as it triggers "dummy_g_ may be used uninitialized" compiler warning for nothing
  /**
   * @brief Move constructor.
   */
  Degree_rips_bifiltration(Degree_rips_bifiltration &&other) noexcept : generators_(std::move(other.generators_)) {}

  /**
   * @brief Copy constructor.
   *
   * @tparam U Type convertible into `T`.
   */
  template <typename U, bool OtherCo, bool OtherEnsure1Criticality>
  Degree_rips_bifiltration(const Degree_rips_bifiltration<U, OtherCo, OtherEnsure1Criticality> &other)
      : generators_(other.begin(), other.end())
  {
    if constexpr (Ensure1Criticality && !OtherEnsure1Criticality) {
      if (generators_.size() > 1) throw std::logic_error("Multiparameter filtration value is not 1-critical.");
    }
  }

  ~Degree_rips_bifiltration() = default;

  // cannot use = default as it triggers "dummy_g_ may be used uninitialized" compiler warning for nothing
  /**
   * @brief Assign operator.
   */
  Degree_rips_bifiltration &operator=(const Degree_rips_bifiltration &other) {
    generators_ = other.generators_;
    return *this;
  }

  // cannot use = default as it triggers "dummy_g_ may be used uninitialized" compiler warning for nothing
  /**
   * @brief Move assign operator.
   */
  Degree_rips_bifiltration &operator=(Degree_rips_bifiltration &&other) noexcept {
    generators_ = std::move(other.generators_);
    return *this;
  }

  /**
   * @brief Assign operator.
   *
   * @tparam U Type convertible into `T`.
   */
  template <typename U, bool OtherCo, bool OtherEnsure1Criticality>
  Degree_rips_bifiltration &operator=(const Degree_rips_bifiltration<U, OtherCo, OtherEnsure1Criticality> &other)
  {
    if constexpr (Ensure1Criticality && !OtherEnsure1Criticality) {
      if (other.num_generators() > 1) throw std::logic_error("Multiparameter filtration value is not 1-critical.");
    }
    generators_ = Underlying_container(other.begin(), other.end());
    return *this;
  }

  /**
   * @brief Swap operator.
   */
  friend void swap(Degree_rips_bifiltration &f1, Degree_rips_bifiltration &f2) noexcept
  {
    f1.generators_.swap(f2.generators_);
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
   *
   * The reference returned when `p` is 1 can be modified but will have no impact on the value of the second parameter
   * represented on the class and is shared by the second parameter of all generators. If there is a need to store this
   * value, a copy should be stored instead.
   */
  reference operator()(size_type g, size_type p)
  {
    GUDHI_CHECK(g < generators_.size() && p < 2, std::out_of_range("Out of bound index."));
    if (p == 1) {
      dummy_g_ = g;
      return dummy_g_;
    }
    return generators_[g];
  }

  /**
   * @brief Returns const reference to value of parameter `p` of generator `g`.
   */
  const_reference operator()(size_type g, size_type p) const
  {
    GUDHI_CHECK(g < generators_.size() && p < 2, std::out_of_range("Out of bound index."));
    if (p == 1) {
      dummy_g_ = g;
      return dummy_g_;
    }
    return generators_[g];
  }

  /**
   * @brief Let \f$ g \f$ be the first value in `indices` and \f$ p \f$ the second value.
   * Returns reference to value of parameter \f$ p \f$ of generator \f$ g \f$.
   *
   * The reference returned when `p` is 1 can be modified but will have no impact on the value of the second parameter
   * represented on the class and is shared by the second parameter of all generators. If there is a need to store this
   * value, a copy should be stored instead.
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
    return this->operator()(g, *(++it));
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
    return this->operator()(g, *(++it));
  }

  /**
   * @brief Returns an iterator pointing the begining of the underlying container. The element `val_i` at index `i`
   * corresponds to the first parameter of the generator `(val_i, i)`.
   */
  iterator begin() noexcept { return generators_.begin(); }

  /**
   * @brief Returns an iterator pointing the begining of the underlying container. The element `val_i` at index `i`
   * corresponds to the first parameter of the generator `(val_i, i)`.
   */
  const_iterator begin() const noexcept { return generators_.begin(); }

  /**
   * @brief Returns an iterator pointing the begining of the underlying container. The element `val_i` at index `i`
   * corresponds to the first parameter of the generator `(val_i, i)`.
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
   * The element `val_i` at index `i` corresponds to the first parameter of the generator `(val_i, i)`.
   */
  reverse_iterator rbegin() noexcept { return generators_.rbegin(); }

  /**
   * @brief Returns a reverse iterator pointing to the first element from the back of the underlying container.
   * The element `val_i` at index `i` corresponds to the first parameter of the generator `(val_i, i)`.
   */
  const_reverse_iterator rbegin() const noexcept { return generators_.rbegin(); }

  /**
   * @brief Returns a reverse iterator pointing to the first element from the back of the underlying container.
   * The element `val_i` at index `i` corresponds to the first parameter of the generator `(val_i, i)`.
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
   * @brief Returns the size of the underlying container. Corresponds exactly to @ref num_generators(), but enables
   * to use the class as a classic range with a `begin`, `end` and `size` method.
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
      generators_.reserve(number_of_generators);
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
  Degree_rips_bifiltration<U, OCo, OEns> as_type() const
  {
    std::vector<U> out(generators_.begin(), generators_.end());
    return Degree_rips_bifiltration<U, OCo, OEns>(std::move(out), num_parameters());
  }

  /**
   * @brief Converts the filtration value to @ref Multi_parameter_filtration without any set simplification.
   * @warning The filtration value is converted one to one and is not simplified to a minimal set of generators or
   * ordered by lexicographical order, that undefines the behaviour of some methods of the class.
   * Use @ref as_type(const Degree_rips_bifiltration&) instead if needed.
   */
  Multi_parameter_filtration<T, Co, Ensure1Criticality> convert_to_non_simplified_multi_parameter_filtration() const
  {
    std::vector<T> out(generators_.size() * 2);
    size_type i = 0;
    for (size_type g = 0; g < generators_.size(); ++g) {
      out[i] = generators_[g];
      out[i + 1] = g;
      i += 2;
    }
    return Multi_parameter_filtration<T, Co, Ensure1Criticality>(std::move(out), 2);
  }

  /**
   * @brief Converts the filtration value to @ref Dynamic_multi_parameter_filtration without any set simplification.
   * @warning The filtration value is converted one to one and is not simplified to a minimal set of generators or
   * ordered by lexicographical order, that undefines the behaviour of some methods of the class.
   * Use @ref as_type(const Degree_rips_bifiltration&) instead if needed.
   */
  Dynamic_multi_parameter_filtration<T, Co, Ensure1Criticality>
  convert_to_non_simplified_dynamic_multi_parameter_filtration() const
  {
    std::vector<Multi_parameter_generator<T> > out;
    out.reserve(generators_.size());
    for (size_type g = 0; g < generators_.size(); ++g) {
      std::vector<T> v = {generators_[g], static_cast<T>(g)};
      out.emplace_back(std::move(v));
    }
    return Dynamic_multi_parameter_filtration<T, Co, Ensure1Criticality>(std::move(out), 2);
  }

  // ACCESS

  /**
   * @brief Returns the number of parameters in the filtration value.
   */
  static constexpr size_type num_parameters() { return 2; }

  /**
   * @brief Returns the number of generators in the filtration value, i.e. the criticality of the element.
   */
  size_type num_generators() const { return generators_.size(); }

  /**
   * @brief Returns the total number of values in the filtration value, that is,
   * @ref num_parameters() * @ref num_generators().
   */
  size_type num_entries() const { return generators_.size() * 2; }

  /**
   * @brief Returns a filtration value for which @ref is_plus_inf() returns `true`. Throws if `Co` is true.
   */
  static Degree_rips_bifiltration inf(int number_of_parameters = 2)
  {
    if constexpr (Co) {
      throw std::logic_error("No biggest value possible for Co-filtrations yet.");
    } else {
      return Degree_rips_bifiltration(number_of_parameters, T_inf);
    }
  }

  /**
   * @brief Returns a filtration value for which @ref is_minus_inf() returns `true`.
   */
  static Degree_rips_bifiltration minus_inf(int number_of_parameters = 2)
  {
    return Degree_rips_bifiltration(number_of_parameters, T_m_inf);
  }

  /**
   * @brief Returns a filtration value for which @ref is_nan() returns `true`.
   */
  static constexpr Degree_rips_bifiltration nan([[maybe_unused]] int number_of_parameters = 2)
  {
    return Degree_rips_bifiltration(Gudhi::simplex_tree::empty_filtration_value_t());
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
  constexpr bool is_plus_inf() const
  {
    if constexpr (Co) {
      return false;
    } else {
      if (generators_.empty()) return false;
      for (const T &v : generators_) {
        if (v != T_inf) return false;
      }
      return true;
    }
  }

  /**
   * @brief Returns `true` if and only if the filtration value is considered as minus infinity.
   */
  constexpr bool is_minus_inf() const
  {
    if constexpr (Co) {
      return generators_.size() == 1 && generators_[0] == T_m_inf;
    } else {
      return !generators_.empty() && generators_[0] == T_m_inf;
    }
  }

  /**
   * @brief Returns `true` if and only if the filtration value is considered as NaN.
   */
  constexpr bool is_nan() const { return generators_.empty(); }

  /**
   * @brief Returns `true` if and only if the filtration value is non-empty and is not considered as plus infinity,
   * minus infinity or NaN.
   */
  bool is_finite() const
  {
    if constexpr (Co) {
      return !generators_.empty() && (generators_.size() != 1 || generators_[0] != T_m_inf);
    } else {
      if (generators_.empty() || generators_[0] == T_m_inf) return false;
      for (const T &v : generators_) {
        if (v != T_inf) return true;
      }
      return false;
    }
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
  friend bool is_strict_less_than_lexicographically(const Degree_rips_bifiltration &a,
                                                    const Degree_rips_bifiltration &b)
  {
    if (&a == &b) return false;
    if (a.is_nan()) return false;
    if (b.is_nan()) return true;

    if constexpr (inverse) {
      if (a.num_generators() != b.num_generators()) {
        if (a.num_generators() == 0) return true;
        if (b.num_generators() == 0) return false;
        if (a.generators_[0] < b.generators_[0]) return true;
        if (b.generators_[0] < a.generators_[0]) return false;
        return a.num_generators() < b.num_generators();
      }
    }

    for (std::size_t i = 0U; i < std::min(a.num_generators(), b.num_generators()); ++i) {
      if constexpr (inverse) i = std::min(a.num_generators(), b.num_generators()) - 1 - i;
      if (_is_nan(a.generators_[i]) && !_is_nan(b.generators_[i])) return false;
      if (_is_nan(b.generators_[i])) return true;
      if (a.generators_[i] < b.generators_[i]) return true;
      if (b.generators_[i] < a.generators_[i]) return false;
      if constexpr (inverse) i = std::min(a.num_generators(), b.num_generators()) - 1 - i;
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
  friend bool is_less_or_equal_than_lexicographically(const Degree_rips_bifiltration &a,
                                                      const Degree_rips_bifiltration &b)
  {
    if (&a == &b) return true;
    if (b.is_nan()) return true;
    if (a.is_nan()) return false;

    if constexpr (inverse) {
      if (a.num_generators() != b.num_generators()) {
        if (a.num_generators() == 0) return true;
        if (b.num_generators() == 0) return false;
        if (a.generators_[0] < b.generators_[0]) return true;
        if (b.generators_[0] < a.generators_[0]) return false;
        return a.num_generators() < b.num_generators();
      }
    }

    for (std::size_t i = 0U; i < std::min(a.num_generators(), b.num_generators()); ++i) {
      if constexpr (inverse) i = std::min(a.num_generators(), b.num_generators()) - 1 - i;
      if (_is_nan(a.generators_[i]) && !_is_nan(b.generators_[i])) return false;
      if (_is_nan(b.generators_[i])) return true;
      if (a.generators_[i] < b.generators_[i]) return true;
      if (b.generators_[i] < a.generators_[i]) return false;
      if constexpr (inverse) i = std::min(a.num_generators(), b.num_generators()) - 1 - i;
    }
    return a.num_generators() <= b.num_generators();
  }

  /**
   * @brief Returns `true` if and only if the cones generated by @p b are strictly contained in the
   * cones generated by @p a (recall that the cones are positive if `Co` is false and negative if `Co` is true).
   *
   * Note that not all filtration values are comparable. That is, \f$ a < b \f$ and \f$ b < a \f$ returning both false
   * does **not** imply \f$ a == b \f$. If a total order is needed, use @ref is_strict_less_than_lexicographically
   * instead.
   */
  friend bool operator<(const Degree_rips_bifiltration &a, const Degree_rips_bifiltration &b)
  {
    if (&a == &b) return false;
    if (a.generators_.size() == 0 || b.generators_.size() == 0) return false;
    return _compare_strict(0, a.generators_, b.generators_, a.generators_[0]);
  }

  /**
   * @brief Returns `true` if and only if the cones generated by @p a are strictly contained in the
   * cones generated by @p b (recall that the cones are positive if `Co` is false and negative if `Co` is true).
   *
   * Note that not all filtration values are comparable. That is, \f$ a \le b \f$ and \f$ b \le a \f$ can both return
   * `false`. If a total order is needed, use @ref is_less_or_equal_than_lexicographically instead.
   */
  friend bool operator<=(const Degree_rips_bifiltration &a, const Degree_rips_bifiltration &b)
  {
    if (a.generators_.size() == 0 || b.generators_.size() == 0) return false;
    if (&a == &b) return true;
    return _compare(0, a.generators_, b.generators_, a.generators_[0]);
  }

  /**
   * @brief Returns `true` if and only if the cones generated by @p b are contained in or are (partially)
   * equal to the cones generated by @p a (recall that the cones are positive if `Co` is false and negative if `Co` is
   * true).
   *
   * Note that not all filtration values are comparable. That is, \f$ a > b \f$ and \f$ b > a \f$ returning both false
   * does **not** imply \f$ a == b \f$. If a total order is needed, use @ref is_strict_less_than_lexicographically
   * instead.
   */
  friend bool operator>(const Degree_rips_bifiltration &a, const Degree_rips_bifiltration &b) { return b < a; }

  /**
   * @brief Returns `true` if and only if the cones generated by @p a are contained in or are (partially)
   * equal to the cones generated by @p b (recall that the cones are positive if `Co` is false and negative if `Co` is
   * true).
   *
   * Note that not all filtration values are comparable. That is, \f$ a \ge b \f$ and \f$ b \ge a \f$ can both return
   * `false`. If a total order is needed, use @ref is_less_or_equal_than_lexicographically instead.
   */
  friend bool operator>=(const Degree_rips_bifiltration &a, const Degree_rips_bifiltration &b) { return b <= a; }

  /**
   * @brief Returns `true` if and only if for each \f$ i,j \f$, \f$ a(i,j) \f$ is equal to \f$ b(i,j) \f$.
   *
   * @warning The method considers different two filtration values with different generators, but as the set of
   * generators is rarely minimal, it is still possible that those two values are equivalent.
   */
  friend bool operator==(const Degree_rips_bifiltration &a, const Degree_rips_bifiltration &b)
  {
    if (a.is_nan() || b.is_nan()) return false;
    if (&a == &b) return true;
    return a.generators_ == b.generators_;
  }

  /**
   * @brief Returns `true` if and only if \f$ a == b \f$ returns `false`.
   */
  friend bool operator!=(const Degree_rips_bifiltration &a, const Degree_rips_bifiltration &b) { return !(a == b); }

  // ARITHMETIC OPERATORS

  // opposite
  /**
   * @brief Returns a filtration value such that an entry at index \f$ i,0 \f$ is equal to \f$ -f(i,0) \f$.
   *
   * Used conventions:
   * - \f$ -NaN = NaN \f$.
   *
   * @param f Value to opposite.
   * @return The opposite of @p f.
   */
  friend Degree_rips_bifiltration operator-(const Degree_rips_bifiltration &f)
  {
    using F = Degree_rips_bifiltration;

    Underlying_container result(f.generators_);
    std::for_each(result.begin(), result.end(), [](T &v) {
      if (v == F::T_inf)
        v = F::T_m_inf;
      else if (v == F::T_m_inf)
        v = F::T_inf;
      else
        v = -v;
    });
    return Degree_rips_bifiltration(std::move(result), Degree_rips_bifiltration::num_parameters());
  }

  // subtraction
  /**
   * @brief Returns a filtration value such that an entry at index \f$ (g,0) \f$, with
   * \f$ 0 \leq g \leq num_generators \f$ is equal to \f$ f(g,0) - r(0) \f$
   * if \f$ 0 < length_r \f$ and to \f$ f(g,0) \f$ otherwise.
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
  friend Degree_rips_bifiltration operator-(Degree_rips_bifiltration f, const ValueRange &r)
  {
    f -= r;
    return f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ (g,0) \f$, with
   * \f$ 0 \leq g \leq num_generators \f$ is equal to \f$ r(0) - f(g,0) \f$
   * if \f$ 0 < length_r \f$ and to -\f$ f(g,0) \f$ otherwise.
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
                                     !std::is_same_v<ValueRange, Degree_rips_bifiltration> > >
  friend Degree_rips_bifiltration operator-(const ValueRange &r, Degree_rips_bifiltration f)
  {
    if (r.begin() == r.end()) return -f;
    return *(r.begin()) - f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ (g,0) \f$ is equal to \f$ f(g,0) - val \f$.
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
  friend Degree_rips_bifiltration operator-(Degree_rips_bifiltration f, const T &val)
  {
    f -= val;
    return f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ (g,0) \f$ is equal to \f$ val - f(g,0) \f$.
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
  friend Degree_rips_bifiltration operator-(const T &val, Degree_rips_bifiltration f)
  {
    f._apply_operation(val, [](T &valF, const T &valR) {
      valF = -valF;
      _add(valF, valR);
    });
    return f;
  }

  /**
   * @brief Modifies the first parameter such that an entry at index \f$ (g,0) \f$, with
   * \f$ 0 \leq g \leq num_generators \f$ is equal to \f$ f(g,0) - r(0) \f$
   * if \f$ 0 < length_r \f$ and to \f$ f(g,0) \f$ otherwise.
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
  friend Degree_rips_bifiltration &operator-=(Degree_rips_bifiltration &f, const ValueRange &r)
  {
    if (r.begin() == r.end()) return f;
    f -= *(r.begin());
    return f;
  }

  /**
   * @brief Modifies the first parameter such that an entry at index \f$ (g,0) \f$ is equal to \f$ f(g,0) - val \f$.
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
  friend Degree_rips_bifiltration &operator-=(Degree_rips_bifiltration &f, const T &val)
  {
    f._apply_operation(val, [](T &valF, const T &valR) { _subtract(valF, valR); });
    return f;
  }

  // addition
  /**
   * @brief Returns a filtration value such that an entry at index \f$ (g,0) \f$, with
   * \f$ 0 \leq g \leq num_generators \f$ is equal to \f$ f(g,0) + r(0) \f$
   * if \f$ 0 < length_r \f$ and to \f$ f(g,0) \f$ otherwise.
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
  friend Degree_rips_bifiltration operator+(Degree_rips_bifiltration f, const ValueRange &r)
  {
    f += r;
    return f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ (g,0) \f$, with
   * \f$ 0 \leq g \leq num_generators \f$ is equal to \f$ r(0) + f(g,0) \f$
   * if \f$ 0 < length_r \f$ and to \f$ f(g,0) \f$ otherwise.
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
                                     !std::is_same_v<ValueRange, Degree_rips_bifiltration> > >
  friend Degree_rips_bifiltration operator+(const ValueRange &r, Degree_rips_bifiltration f)
  {
    f += r;
    return f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ (g,0) \f$ is equal to \f$ f(g,0) + val \f$.
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
  friend Degree_rips_bifiltration operator+(Degree_rips_bifiltration f, const T &val)
  {
    f += val;
    return f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ (g,0) \f$ is equal to \f$ val + f(g,0) \f$.
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
  friend Degree_rips_bifiltration operator+(const T &val, Degree_rips_bifiltration f)
  {
    f += val;
    return f;
  }

  /**
   * @brief Modifies the first parameter such that an entry at index \f$ (g,0) \f$, with
   * \f$ 0 \leq g \leq num_generators \f$ is equal to \f$ f(g,0) + r(0) \f$
   * if \f$ 0 < length_r \f$ and to \f$ f(g,0) \f$ otherwise.
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
  friend Degree_rips_bifiltration &operator+=(Degree_rips_bifiltration &f, const ValueRange &r)
  {
    if (r.begin() == r.end()) return f;
    f += *(r.begin());
    return f;
  }

  /**
   * @brief Modifies the first parameter such that an entry at index \f$ (g,0) \f$ is equal to \f$ f(g,0) + val \f$.
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
  friend Degree_rips_bifiltration &operator+=(Degree_rips_bifiltration &f, const T &val)
  {
    f._apply_operation(val, [](T &valF, const T &valR) { _add(valF, valR); });
    return f;
  }

  // multiplication
  /**
   * @brief Returns a filtration value such that an entry at index \f$ (g,0) \f$, with
   * \f$ 0 \leq g \leq num_generators \f$ is equal to \f$ f(g,0) * r(0) \f$
   * if \f$ 0 < length_r \f$ and to \f$ f(g,0) \f$ otherwise.
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
  friend Degree_rips_bifiltration operator*(Degree_rips_bifiltration f, const ValueRange &r)
  {
    f *= r;
    return f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ (g,0) \f$, with
   * \f$ 0 \leq g \leq num_generators \f$ is equal to \f$ r(0) * f(g,0) \f$
   * if \f$ 0 < length_r \f$ and to \f$ f(g,0) \f$ otherwise.
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
                                     !std::is_same_v<ValueRange, Degree_rips_bifiltration> > >
  friend Degree_rips_bifiltration operator*(const ValueRange &r, Degree_rips_bifiltration f)
  {
    f *= r;
    return f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ (g,0) \f$ is equal to \f$ f(g,0) * val \f$.
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
  friend Degree_rips_bifiltration operator*(Degree_rips_bifiltration f, const T &val)
  {
    f *= val;
    return f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ (g,0) \f$ is equal to \f$ val * f(g,0) \f$.
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
  friend Degree_rips_bifiltration operator*(const T &val, Degree_rips_bifiltration f)
  {
    f *= val;
    return f;
  }

  /**
   * @brief Modifies the first parameter such that an entry at index \f$ (g,0) \f$, with
   * \f$ 0 \leq g \leq num_generators \f$ is equal to \f$ f(g,0) + r(0) \f$
   * if \f$ 0 < length_r \f$ and to \f$ f(g,0) \f$ otherwise.
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
  friend Degree_rips_bifiltration &operator*=(Degree_rips_bifiltration &f, const ValueRange &r)
  {
    if (r.begin() == r.end()) return f;
    f *= *(r.begin());
    return f;
  }

  /**
   * @brief Modifies the first parameter such that an entry at index \f$ (g,0) \f$ is equal to \f$ f(g,0) * val \f$.
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
  friend Degree_rips_bifiltration &operator*=(Degree_rips_bifiltration &f, const T &val)
  {
    f._apply_operation(val, [](T &valF, const T &valR) { _multiply(valF, valR); });
    return f;
  }

  // division
  /**
   * @brief Returns a filtration value such that an entry at index \f$ (g,0) \f$, with
   * \f$ 0 \leq g \leq num_generators \f$ is equal to \f$ f(g,0) / r(0) \f$
   * if \f$ 0 < length_r \f$ and to \f$ f(g,0) \f$ otherwise.
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
  friend Degree_rips_bifiltration operator/(Degree_rips_bifiltration f, const ValueRange &r)
  {
    f /= r;
    return f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ (g,0) \f$, with
   * \f$ 0 \leq g \leq num_generators \f$ is equal to \f$ r(0) / f(g,0) \f$
   * if \f$ 0 < length_r \f$ and to \f$ 1 / f(g,0) \f$ otherwise.
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
                                     !std::is_same_v<ValueRange, Degree_rips_bifiltration> > >
  friend Degree_rips_bifiltration operator/(const ValueRange &r, Degree_rips_bifiltration f)
  {
    if (r.begin() == r.end()) return 1 / f;
    return *(r.begin()) / f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ (g,0) \f$ is equal to \f$ f(g,0) / val \f$.
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
  friend Degree_rips_bifiltration operator/(Degree_rips_bifiltration f, const T &val)
  {
    f /= val;
    return f;
  }

  /**
   * @brief Returns a filtration value such that an entry at index \f$ (g,0) \f$ is equal to \f$ val / f(g,0) \f$.
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
  friend Degree_rips_bifiltration operator/(const T &val, Degree_rips_bifiltration f)
  {
    f._apply_operation(val, [](T &valF, const T &valR) {
      T tmp = valF;
      valF = valR;
      _divide(valF, tmp);
    });
    return f;
  }

  /**
   * @brief Modifies the first parameter such that an entry at index \f$ (g,0) \f$, with
   * \f$ 0 \leq g \leq num_generators \f$ is equal to \f$ f(g,0) / r(0) \f$
   * if \f$ 0 < length_r \f$ and to \f$ f(g,0) \f$ otherwise.
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
  friend Degree_rips_bifiltration &operator/=(Degree_rips_bifiltration &f, const ValueRange &r)
  {
    if (r.begin() == r.end()) return f;
    f /= *(r.begin());
    return f;
  }

  /**
   * @brief Modifies the first parameter such that an entry at index \f$ (g,0) \f$ is equal to \f$ f(g,0) / val \f$.
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
  friend Degree_rips_bifiltration &operator/=(Degree_rips_bifiltration &f, const T &val)
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
   * @param g New number of generators.
   */
  void set_num_generators(size_type g)
  {
    static_assert(!Ensure1Criticality, "Number of generators cannot be set for a 1-critical only filtration value.");

    if (g == 0) return;
    generators_.resize(g, _get_default_value());
  }

  /**
   * @brief Adds the given generator to the filtration value.
   *
   * It is possible that the generator is ignored if the first parameter is overshadowed by an already existing
   * generator with same second parameter. This would mean that adding the given generator will not span more
   * "lifetime" and therefore there is no need to store it.
   *
   * Let \f$ max_idx \$f be the highest second parameter stored so far. If the given second parameter \f$ i \$f to add
   * is strictly higher than \f$ max_idx + 1 \$f, all possible values between \f$ max_idx \$f and \f$ i \$f will also
   * be added and the corresponding first parameters will be initialized with -inf if `Co` is false and with +inf
   * if `Co` is true.
   *
   * @tparam GeneratorRange Range of elements convertible to `T`. Must have a begin(), end() method and the iterator
   * type should satisfy the requirements of the standard `LegacyForwardIterator`.
   * @param x New generator to add. Has to have the 2 parameters.
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
   * @brief Adds the given generator to the filtration value.
   *
   * It is possible that the generator is ignored if the first parameter is overshadowed by an already existing
   * generator with same second parameter. This would mean that adding the given generator will not span more
   * "lifetime" and therefore there is no need to store it.
   *
   * Let \f$ max_idx \$f be the highest second parameter stored so far. If the given second parameter \f$ i \$f to add
   * is strictly higher than \f$ max_idx + 1 \$f, all possible values between \f$ max_idx \$f and \f$ i \$f will also
   * be added and the corresponding first parameters will be initialized with inf if `Co` is false and with -inf
   * if `Co` is true.
   *
   * @tparam Iterator Iterator class satisfying the requirements of the standard `LegacyForwardIterator`.
   * The dereferenced type has to be convertible to `T`.
   * @param genStart Iterator pointing to the begining of the range of two elements.
   * @param genEnd Iterator pointing to the end of the range.
   * @return true If and only if the generator is actually added to the set of generators.
   * @return false Otherwise.
   */
  template <class Iterator>
  bool add_generator(Iterator genStart, [[maybe_unused]] Iterator genEnd)
  {
    GUDHI_CHECK(std::distance(genStart, genEnd) == 2,
                std::invalid_argument("Wrong range size. Should correspond to the number of parameters."));

    const T val = *genStart;
    ++genStart;
    const size_type index = *genStart;

    GUDHI_CHECK(index >= 0, std::invalid_argument("Second parameter has to be a positive index."));

    if (_is_nan(val)) return false;

    if (index < generators_.size()) {
      if (_dominates(val, generators_[index])) return false;
      generators_[index] = val;
      return true;
    }

    if constexpr (Ensure1Criticality) {
      if (index != 0) throw std::logic_error("Cannot add additional generator to a 1-critical only filtration value.");
    }

    generators_.resize(index + 1, _get_default_minus_value());
    generators_[index] = val;
    return true;
  }

  /**
   * @brief Same as @ref add_generator "".
   */
  template <class GeneratorRange = std::initializer_list<T>,
            class = std::enable_if_t<RangeTraits<GeneratorRange>::has_begin> >
  void add_guaranteed_generator(const GeneratorRange &x)
  {
    add_generator(x.begin(), x.end());
  }

  /**
   * @brief Does nothing, only for interface purposes.
   */
  void simplify() {}

  /**
   * @brief Does nothing, only for interface purposes.
   */
  void remove_empty_generators([[maybe_unused]] bool include_infinities = false) {}

  /**
   * @brief Let \f$ g_c \f$ be the positive cone generated by the generator \f$ g \f$ of the filtration value. And let
   * \f$ x_c \f$ be the positive cone generated by the given point \f$ x \f$. For each generator \f$ g \f$ of the
   * filtration value such that \f$ g_c \f$ strictly contains \f$ x_c \f$, pushes it to the border of \f$ x_c \f$.
   *
   * Say \f$ x = (val, idx) \f$. The second parameter \f$ idx \f$ is assumed to be positive.
   * Note also that, if `Co` is true, the filtration value cannot be pushed to \f$ (*,inf) \f$, and that,
   * if `Ensure1Criticality` is true, \f$ idx \f$ has to be \f$ 0 \f$.
   *
   * @tparam GeneratorRange Either a range of into `T` convertible elements with a begin(), end() and size() method,
   * or @ref Degree_rips_bifiltration<U,...> with `U` convertible into `T`.
   * @param x Range towards to push. Has to have 2 elements (all other elements above that count are ignored). If the
   * range is a @ref Degree_rips_bifiltration, the generator at second parameter 0 is chosen.
   * @param exclude_infinite_values If true, values at infinity or minus infinity are not affected.
   * @return true If the filtration value was actually modified.
   * @return false Otherwise.
   */
  template <class GeneratorRange = std::initializer_list<value_type>,
            class = std::enable_if_t<RangeTraits<GeneratorRange>::has_begin> >
  bool push_to_least_common_upper_bound(const GeneratorRange &x, bool exclude_infinite_values = false)
  {
    if (is_nan()) return false;

    size_type start = 0;
    T startVal = 0;
    T xVal;

    if constexpr (RangeTraits<GeneratorRange>::is_multi_filtration) {
      if (x.size() == 0) return false;  // nan
      xVal = x(0,0);
    } else {
      GUDHI_CHECK(x.size() == 2,
                  std::invalid_argument("Wrong range size. Should correspond to the number of parameters."));

      auto it = x.begin();
      xVal = *it;
      ++it;
      startVal = *it;
      if (startVal < 0) throw std::invalid_argument("Second parameter has to be a positive integer.");
    }

    if (startVal == T_inf) {
      if constexpr (Co) {
        throw std::invalid_argument("Filtration values with negative cones cannot be pushed to (*,inf).");
      }
      if (is_plus_inf()) return false;
      // not exactly true, but the closest feasible representation
      *this = inf();
      return true;
    }
    start = startVal;
    if ((start == 0 && xVal == T_m_inf) || _is_nan(xVal)) return false;

    if constexpr (Ensure1Criticality) {
      // TODO: we could at some point allow [(inf,0),...,(inf,n-1),(v,n)] as 1-critical, but it would make everything
      // more complicated and this class does not seem to be the right option to handle those cases anyway...
      if (start != 0)
        throw std::invalid_argument(
            "Pushing to an index other than 0 or inf is not permitted with Ensure1Criticality at true.");
    }

    bool modified = false;
    T newVal = _get_default_minus_value();

    for (size_type i = 0; i < std::min(start, generators_.size()); ++i){
      auto& v = generators_[i];
      if (!_is_nan(v) && (!exclude_infinite_values || (v != T_inf && v != T_m_inf))){
        T threshold = std::max(v, xVal);
        if (_dominates(newVal, threshold)) newVal = threshold;
        if (v != _get_default_minus_value()) {
          modified = true;
          v = _get_default_minus_value();
        }
      }
    }

    if (start < generators_.size()){
      auto& vs = generators_[start];
      T threshold = std::max(vs, xVal);
      if (_dominates(newVal, threshold)) newVal = threshold;
      if (newVal != vs) {
        modified = true;
        vs = newVal;
      }

      for (size_type i = start + 1; i < generators_.size(); ++i) {
        auto& v = generators_[i];
        if (!_is_nan(v) && v < xVal && (!exclude_infinite_values || (v != T_inf && v != T_m_inf))) {
          modified = true;
          v = xVal;
        }
      }
    } else {
      modified = true;
      generators_.resize(start + 1, _get_default_minus_value());
      generators_[start] = newVal;
    }

    if (is_plus_inf()) *this = inf();

    return modified;
  }

  /**
   * @brief Let \f$ g_c \f$ be the negative cone generated by the generator \f$ g \f$ of the filtration value. And let
   * \f$ x_c \f$ be the negative cone generated by the given point \f$ x \f$. For each generator \f$ g \f$ of the
   * filtration value such that \f$ g_c \f$ strictly contains \f$ x_c \f$, pulls it to the border of \f$ x_c \f$.
   *
   * Say \f$ x = (val, idx) \f$. The second parameter \f$ idx \f$ is assumed to be positive.
   *
   * @tparam GeneratorRange Either a range of into `T` convertible elements with a begin(), end() and size() method,
   * or @ref Degree_rips_bifiltration<U,...> with `U` convertible into `T`.
   * @param x Range towards to push. Has to have 2 elements (all other elements above that count are ignored). If the
   * range is a @ref Degree_rips_bifiltration, the generator at second parameter 0 is chosen.
   * @param exclude_infinite_values If true, values at infinity or minus infinity are not affected.
   * @return true If the filtration value was actually modified.
   * @return false Otherwise.
   */
  template <class GeneratorRange = std::initializer_list<value_type>,
            class = std::enable_if_t<RangeTraits<GeneratorRange>::has_begin> >
  bool pull_to_greatest_common_lower_bound(const GeneratorRange &x, bool exclude_infinite_values = false)
  {
    if (is_nan()) return false;

    size_type end;
    T endVal = 0;
    T xVal;

    if constexpr (RangeTraits<GeneratorRange>::is_multi_filtration) {
      if (x.size() == 0) return false;  // nan
      xVal = x(0,0);
    } else {
      GUDHI_CHECK(x.size() == 2,
                  std::invalid_argument("Wrong range size. Should correspond to the number of parameters."));

      auto it = x.begin();
      xVal = *it;
      ++it;
      endVal = *it;
    }

    if (endVal < 0) throw std::invalid_argument("Second parameter has to be a positive integer.");
    if (endVal == T_inf) end = generators_.size();
    else end = endVal;

    if ((end >= generators_.size() - 1 && xVal == T_inf) || _is_nan(xVal)) return false;

    bool modified = false;
    T newVal = _get_default_minus_value();

    for (size_type i = 0; i < std::min(end, generators_.size()); ++i) {
      auto &v = generators_[i];
      if (!_is_nan(v) && v > xVal && (!exclude_infinite_values || (v != T_inf && v != T_m_inf))) {
        modified = true;
        v = xVal;
      }
    }

    if (end < generators_.size()) {
      for (size_type i = end; i < generators_.size(); ++i) {
        auto v = generators_[i];
        if (!_is_nan(v) && (!exclude_infinite_values || (v != T_inf && v != T_m_inf))){
          T threshold = std::min(v, xVal);
          if (_dominates(newVal, threshold)) newVal = threshold;
        }
      }

      modified |= (end + 1) != generators_.size();
      generators_.resize(end + 1);

      modified |= generators_[end] != newVal;
      generators_[end] = newVal;
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
   * therefore all together a 2D grid. The range of the second parameter has to start at 0 and continue continuously.
   *
   * @tparam OneDimArray A range of values convertible into `T` ordered by increasing value. Has to implement
   * a begin, end and operator[] method.
   * @param grid Vector of @p OneDimArray with size at least 2.
   * @param coordinate If true, the values are set to the coordinates of the projection in the grid. If false,
   * the values are set to the values at the coordinates of the projection.
   */
  template <typename OneDimArray>
  void project_onto_grid(const std::vector<OneDimArray> &grid, bool coordinate = true)
  {
    GUDHI_CHECK(
        grid.size() >= 2,
        std::invalid_argument("The grid should not be smaller than the number of parameters in the filtration value."));

    GUDHI_CHECK_code(const OneDimArray &indices = grid[1]);
    const OneDimArray &values = grid[0];
    for (size_type g = 0; g < num_generators(); ++g) {
      GUDHI_CHECK_code(GUDHI_CHECK(static_cast<size_type>(indices[g]) == g, std::invalid_argument("Unvalid grid.")));

      auto v = static_cast<typename OneDimArray::value_type>(generators_[g]);
      auto d = std::distance(values.begin(), std::lower_bound(values.begin(), values.end(), v));
      if (d != 0 && std::abs(v - values[d]) > std::abs(v - values[d - 1])) {
        --d;
      }
      generators_[g] = coordinate ? static_cast<T>(d) : static_cast<T>(values[d]);
    }
  }

  // FONCTIONNALITIES

  /**
   * @brief Returns the filtration value that is the greatest lower bound of all generators.
   */
  friend Degree_rips_bifiltration factorize_below(const Degree_rips_bifiltration &f)
  {
    if (f.num_generators() <= 1) return f;

    bool nan = true;
    size_type idx = f.num_generators();
    T val = T_inf;

    if constexpr (Co) {
      // -inf are "non existing" generators if not all of them are at -inf
      bool inf = true;
      for (size_type g = 0; g < f.num_generators(); ++g) {
        T v = f.generators_[g];
        if (!_is_nan(v) && v != T_m_inf) {
          nan = false;
          inf = false;
          val = v < val ? v : val;
          idx = g < idx ? g : idx;
        } else if (_is_nan(v)) {
          inf = false;
        } else {
          nan = false;
        }
      }
      if (inf) return Degree_rips_bifiltration::minus_inf();
    } else {
      idx = 0;
      for (const T &v : f.generators_) {
        if (!_is_nan(v)) {
          nan = false;
          val = v < val ? v : val;
        }
      }
    }

    if (nan) return Degree_rips_bifiltration::nan();

    Underlying_container result(idx + 1, _get_default_minus_value());
    result[idx] = val;
    return Degree_rips_bifiltration(std::move(result), 2);
  }

  /**
   * @brief Returns the filtration value that is the least upper bound of all generators.
   */
  friend Degree_rips_bifiltration factorize_above(const Degree_rips_bifiltration &f)
  {
    if (f.num_generators() <= 1) return f;

    bool nan = true;
    size_type idx = 0;
    T val = T_m_inf;

    if constexpr (Co) {
      idx = f.num_generators() - 1;
      for (const T &v : f.generators_) {
        if (!_is_nan(v)) {
          nan = false;
          val = v > val ? v : val;
        }
      }
    } else {
      // inf are "non existing" generators if not all of them are at inf
      bool inf = true;
      for (size_type g = 0; g < f.num_generators(); ++g) {
        T v = f.generators_[g];
        if (!_is_nan(v) && v != T_inf) {
          nan = false;
          inf = false;
          val = v > val ? v : val;
          idx = g > idx ? g : idx;
        } else if (_is_nan(v)) {
          inf = false;
        } else {
          nan = false;
        }
      }
      if (inf) return Degree_rips_bifiltration::inf();
    }

    if (nan) return Degree_rips_bifiltration::nan();

    Underlying_container result(idx + 1, _get_default_minus_value());
    result[idx] = val;
    return Degree_rips_bifiltration(std::move(result), 2);
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
  friend U compute_linear_projection(const Degree_rips_bifiltration &f, const std::vector<U> &x)
  {
    auto project_generator = [&](size_type g) -> U {
      U projection = 0;
      std::size_t size = std::min(x.size(), Degree_rips_bifiltration::num_parameters());
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
  friend U compute_euclidean_distance_to(const Degree_rips_bifiltration &f, const Degree_rips_bifiltration &other)
  {
    // TODO: verify if this really makes a differences in the 1-critical case, otherwise just keep the general case
    if constexpr (Ensure1Criticality) {
      return _compute_frobenius_norm(Degree_rips_bifiltration::num_parameters(),
                                     [&](size_type p) -> T { return f(0, p) - other(0, p); });
    } else {
      U res = std::numeric_limits<U>::max();
      for (size_type g1 = 0; g1 < f.num_generators(); ++g1) {
        for (size_type g2 = 0; g2 < other.num_generators(); ++g2) {
          // Euclidean distance as a Frobenius norm with matrix 1 x p and values 'f(g1, p) - other(g2, p)'
          // Order in the min important to spread possible NaNs
          res = std::min(_compute_frobenius_norm(Degree_rips_bifiltration::num_parameters(),
                                                 [&](size_type p) -> T { return f(g1, p) - other(g2, p); }),
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
  friend U compute_norm(const Degree_rips_bifiltration &f)
  {
    // Frobenius norm with matrix g x p based on Euclidean norm

    if (f.num_generators() == 1) return f.generators_[0];

    U out = 0;
    for (size_type g = 0; g < f.num_generators(); ++g) {
      out += g * g;
      out += f.generators_[g] * f.generators_[g];
    }

    if constexpr (std::is_integral_v<U>) {
      // to avoid Windows issue that don't know how to cast integers for cmath methods
      return std::sqrt(static_cast<double>(out));
    } else {
      return std::sqrt(out);
    }
  }

  /**
   * @brief Computes the coordinates in the given grid, corresponding to the nearest upper bounds of the entries
   * in the given filtration value.
   * The grid has to be represented as a vector of vectors of ordered values convertible into `OutValue`. An index
   * \f$ i \f$ of the vector corresponds to the same parameter as the index \f$ i \f$ in a generator of the filtration
   * value. The ranges correspond to the possible values of the parameters, ordered by increasing value, forming
   * therefore all together a 2D grid. The range of the second parameter has to start at 0 and continue continuously.
   *
   * @tparam OutValue Signed arithmetic type. Default value: std::int32_t.
   * @tparam U Type which is convertible into `OutValue`.
   * @param f Filtration value to project.
   * @param grid Vector of vectors to project into.
   * @return Filtration value \f$ out \f$ whose entry correspond to the indices of the projected values. That is,
   * the projection of \f$ f(g,p) \f$ is \f$ grid[p][out(g,p)] \f$.
   */
  template <typename OutValue = std::int32_t, typename U = T>
  friend Degree_rips_bifiltration<OutValue, Co, Ensure1Criticality> compute_coordinates_in_grid(
      Degree_rips_bifiltration f,
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
  friend Degree_rips_bifiltration<U, Co, Ensure1Criticality> evaluate_coordinates_in_grid(
      const Degree_rips_bifiltration &f,
      const std::vector<std::vector<U> > &grid)
  {
    GUDHI_CHECK(grid.size() >= f.num_parameters(),
                std::invalid_argument(
                    "The size of the grid should correspond to the number of parameters in the filtration value."));

    U grid_inf = Degree_rips_bifiltration<U, Co, Ensure1Criticality>::T_inf;
    std::vector<U> outVec(f.num_generators());

    GUDHI_CHECK_code(const std::vector<U> &indices = grid[1]);
    const std::vector<U> &values = grid[0];
    for (size_type g = 0; g < f.num_generators(); ++g) {
      GUDHI_CHECK_code(GUDHI_CHECK(static_cast<size_type>(indices[g]) == g, std::invalid_argument("Unvalid grid.")));

      const T &c = f.generators_[g];
      outVec[g] = (c == T_inf ? grid_inf : values[c]);
    }

    return Degree_rips_bifiltration<U, Co, Ensure1Criticality>(std::move(outVec),
                                                               Degree_rips_bifiltration::num_parameters());
  }

  // UTILITIES

  /**
   * @brief Outstream operator.
   */
  friend std::ostream &operator<<(std::ostream &stream, const Degree_rips_bifiltration &f)
  {
    const size_type num_gen = f.num_generators();
    const size_type num_param = Degree_rips_bifiltration::num_parameters();

    stream << "( k = " << num_gen << " ) [ ";
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
  friend std::istream &operator>>(std::istream &stream, Degree_rips_bifiltration &f)
  {
    size_type num_gen;
    char delimiter;
    stream >> delimiter;  // (
    stream >> delimiter;  // k
    stream >> delimiter;  // =
    if (delimiter != '=') throw std::invalid_argument("Invalid incoming stream format for Multi_parameter_generator.");
    stream >> num_gen;
    if (!stream.good()) throw std::invalid_argument("Invalid incoming stream format for Multi_parameter_generator.");
    f.generators_.resize(num_gen);
    stream >> delimiter;  // )
    stream >> delimiter;  // [
    if (delimiter != '[') throw std::invalid_argument("Invalid incoming stream format for Multi_parameter_generator.");
    if (num_gen == 0) return stream;
    T val;
    for (size_type i = 0; i < num_gen; ++i) {
      stream >> delimiter;  // [
      val = _get_value<T>(stream);
      f.generators_[i] = val;
      stream >> delimiter;  // ,
      stream >> val;
      if (val != static_cast<T>(i))
        throw std::invalid_argument("Invalid incoming stream format for Multi_parameter_generator.");
      if (!stream.good()) throw std::invalid_argument("Invalid incoming stream format for Multi_parameter_generator.");
      stream >> delimiter;  // ]
      stream >> delimiter;  // ; or last ]
    }
    if (delimiter != ']') throw std::invalid_argument("Invalid incoming stream format for Multi_parameter_generator.");

    return stream;
  }

  /**
   * @brief Returns true if and only if the given filtration value is at plus infinity.
   */
  friend bool is_positive_infinity(const Degree_rips_bifiltration &f)
  {
    return f.is_plus_inf();
  }

  /**
   * @brief Adds the generators of the second argument to the first argument.
   *
   * @param f1 Filtration value to modify.
   * @param f2 Filtration value to merge with the first one. Should have the same number of parameters than the other.
   * @return true If the first argument was actually modified.
   * @return false Otherwise.
   */
  friend bool unify_lifetimes(Degree_rips_bifiltration &f1, const Degree_rips_bifiltration &f2)
  {
    bool modified = false;
    for (size_type g = 0; g < f2.num_generators(); ++g) {
      modified |= f1.add_generator({f2.generators_[g], static_cast<T>(g)});
    }
    return modified;
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
  friend bool intersect_lifetimes(Degree_rips_bifiltration &f1, const Degree_rips_bifiltration &f2)
  {
    if (f2.num_generators() == 0 || f1.num_generators() == 0) {
      Degree_rips_bifiltration res(2, _get_default_minus_value());
      swap(f1, res);
      return f1 != res;
    }

    bool modified = false;
    T threshold1 = f1.generators_[0];
    T threshold2 = f2.generators_[0];
    for (size_type g = 0; g < std::max(f1.num_generators(), f2.num_generators()); ++g) {
      if (g < f1.num_generators())
        threshold1 = _strictly_dominates(threshold1, f1.generators_[g]) ? f1.generators_[g] : threshold1;
      else {
        f1.generators_.push_back(0);
        modified = true;
      }
      if (g < f2.num_generators())
        threshold2 = _strictly_dominates(threshold2, f2.generators_[g]) ? f2.generators_[g] : threshold2;
      if (_strictly_dominates(threshold2, threshold1)) {
        if (f1.generators_[g] != threshold2) modified = true;
        f1.generators_[g] = threshold2;
      } else {
        if (f1.generators_[g] != threshold1) modified = true;
        f1.generators_[g] = threshold1;
      }
    }

    return modified;
  }

  /**
   * @brief Serialize given value into the buffer at given pointer.
   *
   * @param value Value to serialize.
   * @param start Pointer to the start of the space in the buffer where to store the serialization.
   * @return End position of the serialization in the buffer.
   */
  friend char *serialize_value_to_char_buffer(const Degree_rips_bifiltration &value, char *start)
  {
    const size_type length = value.generators_.size();
    const std::size_t arg_size = sizeof(T) * length;
    const std::size_t type_size = sizeof(size_type);
    memcpy(start, &length, type_size);
    memcpy(start + type_size, value.generators_.data(), arg_size);
    return start + arg_size + type_size;
  }

  /**
   * @brief Deserialize the value from a buffer at given pointer and stores it in given value.
   *
   * @param value Value to fill with the deserialized filtration value.
   * @param start Pointer to the start of the space in the buffer where the serialization is stored.
   * @return End position of the serialization in the buffer.
   */
  friend const char *deserialize_value_from_char_buffer(Degree_rips_bifiltration &value, const char *start)
  {
    const std::size_t type_size = sizeof(size_type);
    size_type length;
    memcpy(&length, start, type_size);
    std::size_t arg_size = sizeof(T) * length;
    value.generators_.resize(length);
    memcpy(value.generators_.data(), start + type_size, arg_size);
    return start + arg_size + type_size;
  }

  /**
   * @brief Returns the serialization size of the given filtration value.
   */
  friend std::size_t get_serialization_size_of(const Degree_rips_bifiltration &value)
  {
    return sizeof(size_type) + (sizeof(T) * value.generators_.size());
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
  mutable T dummy_g_;               /**< Dummy value for the first parameter, such that a reference can be returned. */

  /**
   * @brief Default value of an element in the filtration value.
   */
  constexpr static T _get_default_value() { return Co ? T_inf : T_m_inf; }
  constexpr static T _get_default_minus_value() { return Co ? T_m_inf : T_inf; }

  // <
  static bool _compare_strict(size_type i, const Underlying_container &a, const Underlying_container &b, T threshold)
  {
    if (i >= b.size()) return true;
    if (_is_nan(b[i])) return false;

    if (i >= a.size() && _strictly_dominates(threshold, b[i])) return false;
    if (i < a.size()) {
      if (_is_nan(a[i])) return false;
      threshold = _strictly_dominates(threshold, a[i]) ? a[i] : threshold;
      if (a[i] == threshold && b[i] == threshold) return false;
    }
    if (_strictly_dominates(threshold, b[i])) return false;

    return _compare_strict(i + 1, a, b, threshold);
  }

  // <=
  static bool _compare(size_type i, const Underlying_container &a, const Underlying_container &b, T threshold)
  {
    if (i >= b.size()) return true;
    if (_is_nan(b[i])) return false;

    if (i >= a.size() && _strictly_dominates(threshold, b[i])) return false;
    if (i < a.size()) {
      if (_is_nan(a[i])) return false;
      threshold = _strictly_dominates(threshold, a[i]) ? a[i] : threshold;
    }
    if (_strictly_dominates(threshold, b[i])) return false;

    return _compare(i + 1, a, b, threshold);
  }

  static bool _strictly_dominates(T a, T b)
  {
    if constexpr (Co) {
      return a < b;
    } else {
      return a > b;
    }
  }

  static bool _dominates(T a, T b)
  {
    if constexpr (Co) {
      return a <= b;
    } else {
      return a >= b;
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

  template <class F, typename U = T>
  static U _compute_frobenius_norm(size_type number_of_elements, F &&norm)
  {
    if (number_of_elements == 1) return norm(0);

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
class numeric_limits<Gudhi::multi_filtration::Degree_rips_bifiltration<T, Co, Ensure1Criticality> >
{
 public:
  using Filtration_value = Gudhi::multi_filtration::Degree_rips_bifiltration<T, Co, Ensure1Criticality>;

  static constexpr bool has_infinity = !Co;
  static constexpr bool has_quiet_NaN = true;

  static constexpr Filtration_value infinity(std::size_t p = 2) noexcept(!Co) { return Filtration_value::inf(p); };

  // non-standard
  static constexpr Filtration_value minus_infinity(std::size_t p = 2) noexcept
  {
    return Filtration_value::minus_inf(p);
  };

  static constexpr Filtration_value max(std::size_t p = 2) noexcept(!Co)
  {
    if constexpr (Co) {
      throw std::logic_error("No biggest value possible for Co-filtrations yet.");
    } else {
      return Filtration_value(p, std::numeric_limits<T>::max());
    }
  };

  static constexpr Filtration_value lowest(std::size_t p = 2) noexcept { return Filtration_value::minus_inf(p); };

  static constexpr Filtration_value quiet_NaN(std::size_t p = 2) noexcept { return Filtration_value::nan(p); };
};

}  // namespace std

#endif  // MF_DEGREE_RIPS_BIFILTRATION_H_
