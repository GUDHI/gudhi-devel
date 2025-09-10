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
 * @file Dynamic_multi_parameter_filtration.h
 * @author Hannah Schreiber, David Loiseaux
 * @brief Contains the @ref Gudhi::multi_filtration::Dynamic_multi_parameter_filtration class.
 */

#ifndef MF_DYNAMIC_MULTI_PARAMETER_FILTRATION_H_
#define MF_DYNAMIC_MULTI_PARAMETER_FILTRATION_H_

#include <algorithm>    //std::lower_bound
#include <cmath>        //std::isnan, std::min
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

#include <gudhi/Debug_utils.h>
#include <gudhi/Multi_filtration/Multi_parameter_generator.h>
#include <gudhi/Multi_filtration/multi_filtration_utils.h>

namespace Gudhi::multi_filtration {

/**
 * @class Dynamic_multi_parameter_filtration Dynamic_multi_parameter_filtration.h
 * gudhi/Dynamic_multi_parameter_filtration.h
 * @ingroup multi_filtration
 *
 * @brief Class encoding the different generators, i.e., apparition times, of a \f$ k \f$-critical
 * \f$\mathbb R^n\f$-filtration value. E.g., the filtration value of a simplex, or, of the algebraic generator of a
 * module presentation. Different from @ref Multi_parameter_filtration, the underlying container is a vector of vectors
 * and therefore less memory efficient, but much more flexible when modifying the filtration value. So, this class is
 * preferable if a lot of generators need to be added on the fly or removed. But when the filtration value is more or
 * less fixed, e.g. for 1-critical filtrations, we recommend @ref Multi_parameter_filtration instead. Implements
 * the concept @ref FiltrationValue of the @ref Gudhi::Simplex_tree and the concept
 * @ref Gudhi::multi_persistence::MultiFiltrationValue.
 *
 * @details Overloads `std::numeric_limits` such that:
 * - `std::numeric_limits<Dynamic_multi_parameter_filtration>::has_infinity` returns `true`,
 * - `std::numeric_limits<Dynamic_multi_parameter_filtration>::has_quiet_NaN` returns `true`,
 * - `std::numeric_limits<Dynamic_multi_parameter_filtration>::infinity(int)` returns
 * @ref Dynamic_multi_parameter_filtration::inf(int) "",
 * - `std::numeric_limits<Dynamic_multi_parameter_filtration>::minus_infinity(int)` returns
 * @ref Dynamic_multi_parameter_filtration::minus_inf(int) "",
 * - `std::numeric_limits<Dynamic_multi_parameter_filtration>::max(int num_param)` returns a @ref
 * Dynamic_multi_parameter_filtration with one generator of `num_param` parameters evaluated at value
 * `std::numeric_limits<T>::max()`,
 * - `std::numeric_limits<Dynamic_multi_parameter_filtration>::quiet_NaN(int)` returns
 * @ref Dynamic_multi_parameter_filtration::nan(int).
 *
 * Multi-critical filtrations are filtrations such that the lifetime of each object is a union of positive cones in
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
class Dynamic_multi_parameter_filtration
{
 public:
  using Generator = Multi_parameter_generator<T>; /**< Generator type. */
  using Underlying_container = std::vector<Generator>; /**< Underlying container for values. */

  // CONSTRUCTORS

  /**
   * @brief Default constructor. Builds filtration value with one generator and given number of parameters.
   * If Co is false, initializes at -inf, if Co is true, at +inf.
   *
   * @param number_of_parameters If negative, takes the default value instead. Default value: 2.
   */
  Dynamic_multi_parameter_filtration(int number_of_parameters = 2)
      : number_of_parameters_(number_of_parameters < 0 ? 2 : number_of_parameters),
        generators_(1, Generator(1, Co ? T_inf : T_m_inf))
  {}

  /**
   * @brief Builds filtration value with one generator and given number of parameters.
   * All values are initialized at the given value.
   *
   * @warning The generator `{-inf, -inf, ...}`/`{inf, inf, ...}` with \f$ number_of_parameters > 1 \f$ entries is
   * valid but will not benefit from possible optimizations. If those values are not planed to be replaced, it is
   * recommended to use the static methods @ref minus_inf() or @ref inf(), or set `number_of_parameters` to 1, instead.
   *
   * @param number_of_parameters If negative, is set to 2 instead.
   * @param value Initialization value for every value in the generator.
   */
  Dynamic_multi_parameter_filtration(int number_of_parameters, T value)
      : number_of_parameters_(number_of_parameters < 0 ? 2 : number_of_parameters),
        generators_(1, Generator(number_of_parameters_, value))
  {}

  /**
   * @brief Builds filtration value with one generator that is initialized with the given range. The number of
   * parameters are therefore deduced from the length of the range.
   *
   * @tparam ValueRange Range of types convertible to `T`. Should have a begin() and end() method.
   * @param range Values of the generator.
   */
  template <class ValueRange = std::initializer_list<T>, class = std::enable_if_t<RangeTraits<ValueRange>::has_begin> >
  Dynamic_multi_parameter_filtration(const ValueRange &range) : generators_(1, Generator(range.begin(), range.end()))
  {
    number_of_parameters_ = generators_[0].size();
  }

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
  Dynamic_multi_parameter_filtration(Iterator it_begin, Iterator it_end) : generators_(1, Generator(it_begin, it_end))
  {
    number_of_parameters_ = generators_[0].size();
  }

  /**
   * @brief Builds filtration value with given number of parameters and values from the given range. Lets \f$ p \f$
   * be the number of parameters. The \f$ p \f$ first elements of the range have to correspond to the first generator,
   * the \f$ p \f$ next elements to the second generator and so on... So the length of the range has to be a multiple
   * of \f$ p \f$ and the number of generators will be \f$ length / p \f$. The range is represented by two iterators.
   *
   * @tparam Iterator Iterator type that has to satisfy the requirements of standard LegacyForwardIterator and
   * dereferenced elements have to be convertible to `T`.
   * @param it_begin Iterator pointing to the start of the range.
   * @param it_end Iterator pointing to the end of the range.
   * @param number_of_parameters Negative values are associated to 0.
   */
  template <class Iterator, class = std::enable_if_t<!std::is_arithmetic_v<Iterator> > >
  Dynamic_multi_parameter_filtration(Iterator it_begin, Iterator it_end, int number_of_parameters)
      : number_of_parameters_(number_of_parameters < 0 ? 0 : number_of_parameters), generators_()
  {
    // Will discard any value at the end which does not fit into a complete generator.
    const size_type num_gen = std::distance(it_begin, it_end) / number_of_parameters;

    if constexpr (Ensure1Criticality) {
      if (num_gen != 1) throw std::logic_error("Multiparameter filtration value is not 1-critical.");
    }

    Iterator it = it_begin;
    for (size_type i = 0; i < num_gen; ++i) {
      Iterator endIt = it;
      std::advance(endIt, number_of_parameters);
      generators_.emplace_back(it, endIt);
      it = endIt;
    }
  }

  /**
   * @brief Builds filtration value with given number of parameters and values copied from the given
   * @ref Dynamic_multi_parameter_filtration::Underlying_container container.
   *
   * @param generators Values.
   * @param number_of_parameters Negative values are associated to 0.
   */
  Dynamic_multi_parameter_filtration(const Underlying_container &generators, int number_of_parameters)
      : number_of_parameters_(number_of_parameters < 0 ? 0 : number_of_parameters), generators_(generators)
  {
    GUDHI_CHECK(number_of_parameters > 0 || generators_.empty(),
                std::invalid_argument("Number of parameters cannot be 0 if the container is not empty."));

    if constexpr (Ensure1Criticality) {
      if (generators_.size() != 1) throw std::logic_error("Multiparameter filtration value is not 1-critical.");
    }
  }

  /**
   * @brief Builds filtration value with given number of parameters and values moved from the given
   * @ref Dynamic_multi_parameter_filtration::Underlying_container container.
   *
   * @param generators Values to move.
   * @param number_of_parameters Negative values are associated to 0.
   */
  Dynamic_multi_parameter_filtration(Underlying_container &&generators, int number_of_parameters)
      : number_of_parameters_(number_of_parameters < 0 ? 0 : number_of_parameters), generators_(std::move(generators))
  {
    GUDHI_CHECK(number_of_parameters > 0 || generators_.empty(),
                std::invalid_argument("Number of parameters cannot be 0 if the container is not empty."));

    if constexpr (Ensure1Criticality) {
      if (generators_.size() != 1) throw std::logic_error("Multiparameter filtration value is not 1-critical.");
    }
  }

  /**
   * @brief Copy constructor.
   */
  Dynamic_multi_parameter_filtration(const Dynamic_multi_parameter_filtration &other) = default;

  /**
   * @brief Copy constructor.
   *
   * @tparam U Type convertible into `T`.
   */
  template <typename U, bool OtherCo, bool OtherEnsure1Criticality>
  Dynamic_multi_parameter_filtration(
      const Dynamic_multi_parameter_filtration<U, OtherCo, OtherEnsure1Criticality> &other)
      : number_of_parameters_(other.num_parameters()), generators_(other.begin(), other.end())
  {
    if constexpr (Ensure1Criticality && !OtherEnsure1Criticality) {
      if (generators_.size() != 1) throw std::logic_error("Multiparameter filtration value is not 1-critical.");
    }
  }

  /**
   * @brief Move constructor.
   */
  Dynamic_multi_parameter_filtration(Dynamic_multi_parameter_filtration &&other) noexcept = default;

  ~Dynamic_multi_parameter_filtration() = default;

  /**
   * @brief Assign operator.
   */
  Dynamic_multi_parameter_filtration &operator=(const Dynamic_multi_parameter_filtration &other) = default;

  /**
   * @brief Move assign operator.
   */
  Dynamic_multi_parameter_filtration &operator=(Dynamic_multi_parameter_filtration &&other) noexcept = default;

  /**
   * @brief Assign operator.
   *
   * @tparam U Type convertible into `T`.
   */
  template <typename U, bool OtherCo, bool OtherEnsure1Criticality>
  Dynamic_multi_parameter_filtration &operator=(
      const Dynamic_multi_parameter_filtration<U, OtherCo, OtherEnsure1Criticality> &other)
  {
    if constexpr (Ensure1Criticality && !OtherEnsure1Criticality) {
      if (other.num_generators() != 1) throw std::logic_error("Multiparameter filtration value is not 1-critical.");
    }
    generators_ = Underlying_container(other.begin(), other.end());
    number_of_parameters_ = other.num_parameters();
    return *this;
  }

  /**
   * @brief Swap operator.
   */
  friend void swap(Dynamic_multi_parameter_filtration &f1, Dynamic_multi_parameter_filtration &f2) noexcept
  {
    f1.generators_.swap(f2.generators_);
    std::swap(f1.number_of_parameters_, f2.number_of_parameters_);
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
   * If the value is at +/- inf or NaN and it needs to be modified, @ref force_generator_size_to_number_of_parameters
   * needs potentially to be called first such that this methods returns the right reference.
   */
  reference operator()(size_type g, size_type p)
  {
    GUDHI_CHECK(g < generators_.size() && p < number_of_parameters_, std::out_of_range("Out of bound index."));
    if (generators_[g].size() < number_of_parameters_) return generators_[g][0];
    return generators_[g][p];
  }

  /**
   * @brief Returns const reference to value of parameter `p` of generator `g`.
   */
  const_reference operator()(size_type g, size_type p) const
  {
    GUDHI_CHECK(g < generators_.size() && p < number_of_parameters_, std::out_of_range("Out of bound index."));
    if (generators_[g].size() < number_of_parameters_) return generators_[g][0];
    return generators_[g][p];
  }

  /**
   * @brief Returns const reference to the requested generator.
   *
   * @param g Index of the generator.
   */
  const Generator &operator[](size_type g) const
  {
    GUDHI_CHECK(g < generators_.size(), std::out_of_range("Out of bound index."));
    return generators_[g];
  }

  /**
   * @brief Returns reference to the requested generator.
   *
   * @param g Index of the generator.
   */
  Generator &operator[](size_type g)
  {
    GUDHI_CHECK(g < generators_.size(), std::out_of_range("Out of bound index."));
    return generators_[g];
  }

  /**
   * @brief Let \f$ g \f$ be the first value in `indices` and \f$ p \f$ the second value.
   * Returns reference to value of parameter \f$ p \f$ of generator \f$ g \f$.
   * If the value is at +/- inf or NaN and it needs to be modified, @ref force_generator_size_to_number_of_parameters
   * needs potentially to be called first such that this methods returns the right reference.
   *
   * @tparam IndexRange Range with a begin() and size() method.
   * @param indices Range with at least two elements. The first element should correspond to the generator number and
   * the second element to the parameter number.
   */
  template <class IndexRange = std::initializer_list<size_type>,
            class = std::enable_if_t<RangeTraits<IndexRange>::has_begin> >
  reference operator[](const IndexRange &indices)
  {
    GUDHI_CHECK(indices.size() >= 2,
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
    GUDHI_CHECK(indices.size() >= 2,
                std::invalid_argument(
                    "Exactly 2 indices allowed only: first the generator number, second the parameter number."));
    auto it = indices.begin();
    size_type g = *it;
    return this->operator()(g, *(++it));
  }

  /**
   * @brief Returns an iterator pointing the begining of the underlying container. Each element stored is a whole
   * generator with a size() and a operator[].
   *
   * @warning If a generator is modified and the new set of generators is not minimal or not sorted, the behaviour
   * of most methods is undefined. It is possible to call @ref simplify() after construction if there is a doubt to
   * ensure this property.
   */
  iterator begin() noexcept { return generators_.begin(); }

  /**
   * @brief Returns an iterator pointing the begining of the underlying container. Each element stored is a whole
   * generator with a size() and a operator[].
   */
  const_iterator begin() const noexcept { return generators_.begin(); }

  /**
   * @brief Returns an iterator pointing the begining of the underlying container. Each element stored is a whole
   * generator with a size() and a operator[].
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
   * Each element stored is a whole generator with a size() and a operator[].
   *
   * @warning If a generator is modified and the new set of generators is not minimal or not sorted, the behaviour
   * of most methods is undefined. It is possible to call @ref simplify() after construction if there is a doubt to
   * ensure this property.
   */
  reverse_iterator rbegin() noexcept { return generators_.rbegin(); }

  /**
   * @brief Returns a reverse iterator pointing to the first element from the back of the underlying container.
   * Each element stored is a whole generator with a size() and a operator[].
   */
  const_reverse_iterator rbegin() const noexcept { return generators_.rbegin(); }

  /**
   * @brief Returns a reverse iterator pointing to the first element from the back of the underlying container.
   * Each element stored is a whole generator with a size() and a operator[].
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
   * @brief Returns the size of the underlying container. Corresponds exactly to @ref num_generators(), but enables to
   * use the class as a classic range with a `begin`, `end` and `size` method.
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
  Dynamic_multi_parameter_filtration<U, OCo, OEns> as_type() const
  {
    std::vector<Multi_parameter_generator<U> > out(num_generators());
    for (size_type g = 0; g < num_generators(); ++g) {
      out[g] = generators_[g].template as_type<U>();
    }
    return Dynamic_multi_parameter_filtration<U, OCo, OEns>(std::move(out), num_parameters());
  }

  // ACCESS

  /**
   * @brief Returns the number of parameters in the filtration value.
   */
  size_type num_parameters() const { return number_of_parameters_; }

  /**
   * @brief Returns the number of generators in the filtration value, i.e. the criticality of the element.
   */
  size_type num_generators() const
  {
    if constexpr (Ensure1Criticality) {
      return 1;  // for possible optimizations? If there is none, we can just keep the other version
    } else {
      return generators_.size();
    }
  }

  /**
   * @brief Returns the total number of values in the filtration value, that is,
   * @ref num_parameters() * @ref num_generators().
   */
  size_type num_entries() const { return generators_.size() * number_of_parameters_; }

  /**
   * @brief Returns a filtration value with given number of parameters for which @ref is_plus_inf() returns `true`
   * or an empty filtration value if `number_of_parameters` is 0.
   */
  static Dynamic_multi_parameter_filtration inf(int number_of_parameters)
  {
    if (number_of_parameters == 0) return Dynamic_multi_parameter_filtration();
    Underlying_container out(1, Generator::inf());
    return Dynamic_multi_parameter_filtration(std::move(out), number_of_parameters);
  }

  /**
   * @brief Returns a filtration value with given number of parameters for which @ref is_minus_inf() returns `true`
   * or an empty filtration value if `number_of_parameters` is 0.
   */
  static Dynamic_multi_parameter_filtration minus_inf(int number_of_parameters)
  {
    if (number_of_parameters == 0) return Dynamic_multi_parameter_filtration();
    Underlying_container out(1, Generator::minus_inf());
    return Dynamic_multi_parameter_filtration(std::move(out), number_of_parameters);
  }

  /**
   * @brief Returns a filtration value with given number of parameters for which @ref is_nan() returns `true`
   * or an empty filtration value if `number_of_parameters` is 0.
   */
  static Dynamic_multi_parameter_filtration nan(int number_of_parameters)
  {
    if (number_of_parameters == 0) return Dynamic_multi_parameter_filtration();
    Underlying_container out(1, Generator::nan());
    return Dynamic_multi_parameter_filtration(std::move(out), number_of_parameters);
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
    for (const Generator &g : generators_) {
      if (!g.is_plus_inf()) return false;
    }
    return !generators_.empty();
  }

  /**
   * @brief Returns `true` if and only if the filtration value is considered as minus infinity.
   */
  [[nodiscard]] bool is_minus_inf() const
  {
    for (const Generator &g : generators_) {
      if (!g.is_minus_inf()) return false;
    }
    return !generators_.empty();
  }

  /**
   * @brief Returns `true` if and only if the filtration value is considered as NaN.
   */
  [[nodiscard]] bool is_nan() const
  {
    if (generators_.empty()) return false;
    for (const Generator &g : generators_) {
      if (!g.is_nan()) return false;
    }
    return true;
  }

  /**
   * @brief Returns `true` if and only if the filtration value is non-empty and is not considered as plus infinity,
   * minus infinity or NaN.
   */
  [[nodiscard]] bool is_finite() const
  {
    bool isInf = true, isMinusInf = true, isNan = true;
    for (const Generator &g : generators_) {
      if (!g.is_plus_inf()) isInf = false;
      if (!g.is_minus_inf()) isMinusInf = false;
      if (!g.is_nan()) isNan = false;
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
  friend bool is_strict_less_than_lexicographically(const Dynamic_multi_parameter_filtration &a,
                                                    const Dynamic_multi_parameter_filtration &b)
  {
    if (&a == &b) return false;

    GUDHI_CHECK(a.num_parameters() == b.num_parameters(),
                std::invalid_argument("Only filtration values with same number of parameters can be compared."));

    // line order matters
    if (a.is_nan()) return false;
    if (b.is_nan()) return true;
    if (a.is_plus_inf() || b.is_minus_inf()) return false;
    if (a.is_minus_inf() || b.is_plus_inf()) return true;

    // TODO: verify if this really makes a differences in the 1-critical case, otherwise just keep the general case
    if constexpr (Ensure1Criticality) {
      for (std::size_t p = 0U; p < a.num_parameters(); ++p) {
        if constexpr (inverse) p = a.num_parameters() - 1 - p;
        if (_is_nan(a.generators_[0][p]) && !_is_nan(b.generators_[0][p])) return false;
        if (_is_nan(b.generators_[0][p])) return true;
        if (a.generators_[0][p] < b.generators_[0][p]) return true;
        if (b.generators_[0][p] < a.generators_[0][p]) return false;
        if constexpr (inverse) p = a.num_parameters() - 1 - p;
      }
      return false;
    } else {
      for (std::size_t g = 0U; g < std::min(a.num_generators(), b.num_generators()); ++g) {
        std::size_t gA = g;
        std::size_t gB = g;
        if constexpr (inverse) {
          gA = a.num_generators() - 1 - g;
          gB = b.num_generators() - 1 - g;
        }
        for (std::size_t p = 0U; p < a.num_parameters(); ++p) {
          if constexpr (inverse) p = a.num_parameters() - 1 - p;
          if (_is_nan(a.generators_[gA][p]) && !_is_nan(b.generators_[gB][p])) return false;
          if (_is_nan(b.generators_[gB][p])) return true;
          if (a.generators_[gA][p] < b.generators_[gB][p]) return true;
          if (b.generators_[gB][p] < a.generators_[gA][p]) return false;
          if constexpr (inverse) p = a.num_parameters() - 1 - p;
        }
      }
      return a.num_generators() < b.num_generators();
    }
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
  friend bool is_less_or_equal_than_lexicographically(const Dynamic_multi_parameter_filtration &a,
                                                      const Dynamic_multi_parameter_filtration &b)
  {
    if (&a == &b) return true;

    GUDHI_CHECK(a.num_parameters() == b.num_parameters(),
                std::invalid_argument("Only filtration values with same number of parameters can be compared."));

    // line order matters
    if (b.is_nan()) return true;
    if (a.is_nan()) return false;
    if (a.is_minus_inf() || b.is_plus_inf()) return true;
    if (a.is_plus_inf() || b.is_minus_inf()) return false;

    // TODO: verify if this really makes a differences in the 1-critical case, otherwise just keep the general case
    if constexpr (Ensure1Criticality) {
      for (std::size_t p = 0U; p < a.num_parameters(); ++p) {
        if constexpr (inverse) p = a.num_parameters() - 1 - p;
        if (_is_nan(a.generators_[0][p]) && !_is_nan(b.generators_[0][p])) return false;
        if (_is_nan(b.generators_[0][p])) return true;
        if (a.generators_[0][p] < b.generators_[0][p]) return true;
        if (b.generators_[0][p] < a.generators_[0][p]) return false;
        if constexpr (inverse) p = a.num_parameters() - 1 - p;
      }
      return true;
    } else {
      for (std::size_t g = 0U; g < std::min(a.num_generators(), b.num_generators()); ++g) {
        std::size_t gA = g;
        std::size_t gB = g;
        if constexpr (inverse) {
          gA = a.num_generators() - 1 - g;
          gB = b.num_generators() - 1 - g;
        }
        for (std::size_t p = 0U; p < a.num_parameters(); ++p) {
          if constexpr (inverse) p = a.num_parameters() - 1 - p;
          if (_is_nan(a.generators_[gA][p]) && !_is_nan(b.generators_[gB][p])) return false;
          if (_is_nan(b.generators_[gB][p])) return true;
          if (a.generators_[gA][p] < b.generators_[gB][p]) return true;
          if (b.generators_[gB][p] < a.generators_[gA][p]) return false;
          if constexpr (inverse) p = a.num_parameters() - 1 - p;
        }
      }
      return a.num_generators() <= b.num_generators();
    }
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
  friend bool operator<(const Dynamic_multi_parameter_filtration &a, const Dynamic_multi_parameter_filtration &b)
  {
    if (&a == &b) return false;

    GUDHI_CHECK(a.num_parameters() == b.num_parameters(),
                std::invalid_argument("Only filtration values with same number of parameters can be compared."));

    if (a.num_generators() == 0 || b.num_generators() == 0) return false;

    // TODO: verify if this really makes a differences in the 1-critical case, otherwise just keep the general case
    if constexpr (Ensure1Criticality) {
      if (_first_dominates(a.generators_[0], b.generators_[0])) return false;
      return _strictly_contains(a.generators_[0], b.generators_[0]);
    } else {
      for (std::size_t i = 0U; i < b.num_generators(); ++i) {
        // for each generator in b, verify if it is strictly in the cone of at least one generator of a
        bool isContained = false;
        for (std::size_t j = 0U; j < a.num_generators() && !isContained; ++j) {
          // lexicographical order, so if a[j][0] dom b[j][0], than a[j'] can never strictly contain b[i] for all
          // j' > j.
          if (_first_dominates(a.generators_[j], b.generators_[i])) return false;
          isContained = _strictly_contains(a.generators_[j], b.generators_[i]);
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
  friend bool operator<=(const Dynamic_multi_parameter_filtration &a, const Dynamic_multi_parameter_filtration &b)
  {
    GUDHI_CHECK(a.num_parameters() == b.num_parameters(),
                std::invalid_argument("Only filtration values with same number of parameters can be compared."));

    if (a.num_generators() == 0 || b.num_generators() == 0) return false;
    if (a.is_nan() || b.is_nan()) return false;
    if (&a == &b) return true;

    // TODO: verify if this really makes a differences in the 1-critical case, otherwise just keep the general case
    if constexpr (Ensure1Criticality) {
      if (_first_strictly_dominates(a.generators_[0], b.generators_[0])) return false;
      return _contains(a.generators_[0], b.generators_[0]);
    } else {
      // check if this curves is below other's curve
      //  ie for each guy in this, check if there is a guy in other that dominates him
      for (std::size_t i = 0U; i < b.num_generators(); ++i) {
        // for each generator in b, verify if it is in the cone of at least one generator of a
        bool isContained = false;
        for (std::size_t j = 0U; j < a.num_generators() && !isContained; ++j) {
          // lexicographical order, so if a[j][0] strictly dom b[j][0], than a[j'] can never contain b[i] for all
          // j' > j.
          if (_first_strictly_dominates(a.generators_[j], b.generators_[i])) return false;
          isContained = _contains(a.generators_[j], b.generators_[i]);
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
  friend bool operator>(const Dynamic_multi_parameter_filtration &a, const Dynamic_multi_parameter_filtration &b)
  {
    return b < a;
  }

  /**
   * @brief Returns `true` if and only if the cones generated by @p a are contained in or are (partially)
   * equal to the cones generated by @p b (recall that the cones are positive if `Co` is false and negative if `Co` is
   * true).
   * Both @p a and @p b  have to have the same number of parameters.
   *
   * Note that not all filtration values are comparable. That is, \f$ a \ge b \f$ and \f$ b \ge a \f$ can both return
   * `false`. If a total order is needed, use @ref is_less_or_equal_than_lexicographically instead.
   */
  friend bool operator>=(const Dynamic_multi_parameter_filtration &a, const Dynamic_multi_parameter_filtration &b)
  {
    return b <= a;
  }

  /**
   * @brief Returns `true` if and only if for each \f$ i,j \f$, \f$ a(i,j) \f$ is equal to \f$ b(i,j) \f$.
   */
  friend bool operator==(const Dynamic_multi_parameter_filtration &a, const Dynamic_multi_parameter_filtration &b)
  {
    if (a.is_nan() || b.is_nan()) return false;
    if (&a == &b) return true;
    if (a.number_of_parameters_ != b.number_of_parameters_) return false;
    // assumes lexicographical order for both
    return a.generators_ == b.generators_;
  }

  /**
   * @brief Returns `true` if and only if \f$ a == b \f$ returns `false`.
   */
  friend bool operator!=(const Dynamic_multi_parameter_filtration &a, const Dynamic_multi_parameter_filtration &b)
  {
    return !(a == b);
  }

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
  friend Dynamic_multi_parameter_filtration operator-(const Dynamic_multi_parameter_filtration &f)
  {
    Underlying_container result(f.generators_);
    std::for_each(result.begin(), result.end(), [](Generator &v) { v = -v; });
    return Dynamic_multi_parameter_filtration(std::move(result), f.num_parameters());
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
   * @warning The operator accepts @ref Dynamic_multi_parameter_filtration with the same or different template
   * parameters as `ValueRange`. But if the number of generators is higher than 1, only the first generator will be
   * used for the operation.
   *
   * @tparam ValueRange Either a range of into `T` convertible elements with a begin() and end() method,
   * or @ref Dynamic_multi_parameter_filtration<U,...> with `U` convertible into `T`.
   * @param f First element of the subtraction.
   * @param r Second element of the subtraction.
   */
  template <class ValueRange, class = std::enable_if_t<RangeTraits<ValueRange>::has_begin> >
  friend Dynamic_multi_parameter_filtration operator-(Dynamic_multi_parameter_filtration f, const ValueRange &r)
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
   * @warning The operator accepts @ref Dynamic_multi_parameter_filtration with the same or different template
   * parameters as `ValueRange`. But if the number of generators is higher than 1, only the first generator will be
   * used for the operation.
   *
   * @tparam ValueRange Either a range of into `T` convertible elements with a begin() and end() method,
   * or @ref Dynamic_multi_parameter_filtration<U,...> with `U` convertible into `T`.
   * @param r First element of the subtraction.
   * @param f Second element of the subtraction.
   */
  template <class ValueRange,
            class = std::enable_if_t<RangeTraits<ValueRange>::has_begin &&
                                     !std::is_same_v<ValueRange, Dynamic_multi_parameter_filtration> > >
  friend Dynamic_multi_parameter_filtration operator-(const ValueRange &r, Dynamic_multi_parameter_filtration f)
  {
    for (Generator &g : f.generators_) {
      g = r - g;
    }

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
  friend Dynamic_multi_parameter_filtration operator-(Dynamic_multi_parameter_filtration f, const T &val)
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
  friend Dynamic_multi_parameter_filtration operator-(const T &val, Dynamic_multi_parameter_filtration f)
  {
    for (Generator &g : f.generators_) {
      g = val - g;
    }
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
   * @warning The operator accepts @ref Dynamic_multi_parameter_filtration with the same or different template
   * parameters as `ValueRange`. But if the number of generators is higher than 1, only the first generator will be
   * used for the operation.
   *
   * @tparam ValueRange Either a range of into `T` convertible elements with a begin() and end() method,
   * or @ref Dynamic_multi_parameter_filtration<U,...> with `U` convertible into `T`.
   * @param f First element of the subtraction.
   * @param r Second element of the subtraction.
   */
  template <class ValueRange, class = std::enable_if_t<RangeTraits<ValueRange>::has_begin> >
  friend Dynamic_multi_parameter_filtration &operator-=(Dynamic_multi_parameter_filtration &f, const ValueRange &r)
  {
    for (Generator &g : f.generators_) {
      g -= r;
    }
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
  friend Dynamic_multi_parameter_filtration &operator-=(Dynamic_multi_parameter_filtration &f, const T &val)
  {
    for (Generator &g : f.generators_) {
      g -= val;
    }
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
   * @warning The operator accepts @ref Dynamic_multi_parameter_filtration with the same or different template
   * parameters as `ValueRange`. But if the number of generators is higher than 1, only the first generator will be
   * used for the operation.
   *
   * @tparam ValueRange Either a range of into `T` convertible elements with a begin() and end() method,
   * or @ref Dynamic_multi_parameter_filtration<U,...> with `U` convertible into `T`.
   * @param f First element of the addition.
   * @param r Second element of the addition.
   */
  template <class ValueRange, class = std::enable_if_t<RangeTraits<ValueRange>::has_begin> >
  friend Dynamic_multi_parameter_filtration operator+(Dynamic_multi_parameter_filtration f, const ValueRange &r)
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
   * @warning The operator accepts @ref Dynamic_multi_parameter_filtration with the same or different template
   * parameters as `ValueRange`. But if the number of generators is higher than 1, only the first generator will be
   * used for the operation.
   *
   * @tparam ValueRange Either a range of into `T` convertible elements with a begin() and end() method,
   * or @ref Dynamic_multi_parameter_filtration<U,...> with `U` convertible into `T`.
   * @param r First element of the addition.
   * @param f Second element of the addition.
   */
  template <class ValueRange,
            class = std::enable_if_t<RangeTraits<ValueRange>::has_begin &&
                                     !std::is_same_v<ValueRange, Dynamic_multi_parameter_filtration> > >
  friend Dynamic_multi_parameter_filtration operator+(const ValueRange &r, Dynamic_multi_parameter_filtration f)
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
  friend Dynamic_multi_parameter_filtration operator+(Dynamic_multi_parameter_filtration f, const T &val)
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
  friend Dynamic_multi_parameter_filtration operator+(const T &val, Dynamic_multi_parameter_filtration f)
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
   * @warning The operator accepts @ref Dynamic_multi_parameter_filtration with the same or different template
   * parameters as `ValueRange`. But if the number of generators is higher than 1, only the first generator will be
   * used for the operation.
   *
   * @tparam ValueRange Either a range of into `T` convertible elements with a begin() and end() method,
   * or @ref Dynamic_multi_parameter_filtration<U,...> with `U` convertible into `T`.
   * @param f First element of the addition.
   * @param r Second element of the addition.
   */
  template <class ValueRange, class = std::enable_if_t<RangeTraits<ValueRange>::has_begin> >
  friend Dynamic_multi_parameter_filtration &operator+=(Dynamic_multi_parameter_filtration &f, const ValueRange &r)
  {
    for (Generator &g : f.generators_) {
      g += r;
    }
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
  friend Dynamic_multi_parameter_filtration &operator+=(Dynamic_multi_parameter_filtration &f, const T &val)
  {
    for (Generator &g : f.generators_) {
      g += val;
    }
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
   * @warning The operator accepts @ref Dynamic_multi_parameter_filtration with the same or different template
   * parameters as `ValueRange`. But if the number of generators is higher than 1, only the first generator will be
   * used for the operation.
   *
   * @tparam ValueRange Either a range of into `T` convertible elements with a begin() and end() method,
   * or @ref Dynamic_multi_parameter_filtration<U,...> with `U` convertible into `T`.
   * @param f First element of the multiplication.
   * @param r Second element of the multiplication.
   */
  template <class ValueRange, class = std::enable_if_t<RangeTraits<ValueRange>::has_begin> >
  friend Dynamic_multi_parameter_filtration operator*(Dynamic_multi_parameter_filtration f, const ValueRange &r)
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
   * @warning The operator accepts @ref Dynamic_multi_parameter_filtration with the same or different template
   * parameters as `ValueRange`. But if the number of generators is higher than 1, only the first generator will be
   * used for the operation.
   *
   * @tparam ValueRange Either a range of into `T` convertible elements with a begin() and end() method,
   * or @ref Dynamic_multi_parameter_filtration<U,...> with `U` convertible into `T`.
   * @param r First element of the multiplication.
   * @param f Second element of the multiplication.
   */
  template <class ValueRange,
            class = std::enable_if_t<RangeTraits<ValueRange>::has_begin &&
                                     !std::is_same_v<ValueRange, Dynamic_multi_parameter_filtration> > >
  friend Dynamic_multi_parameter_filtration operator*(const ValueRange &r, Dynamic_multi_parameter_filtration f)
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
  friend Dynamic_multi_parameter_filtration operator*(Dynamic_multi_parameter_filtration f, const T &val)
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
  friend Dynamic_multi_parameter_filtration operator*(const T &val, Dynamic_multi_parameter_filtration f)
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
   * @warning The operator accepts @ref Dynamic_multi_parameter_filtration with the same or different template
   * parameters as `ValueRange`. But if the number of generators is higher than 1, only the first generator will be
   * used for the operation.
   *
   * @tparam ValueRange Either a range of into `T` convertible elements with a begin() and end() method,
   * or @ref Dynamic_multi_parameter_filtration<U,...> with `U` convertible into `T`.
   * @param f First element of the multiplication.
   * @param r Second element of the multiplication.
   */
  template <class ValueRange, class = std::enable_if_t<RangeTraits<ValueRange>::has_begin> >
  friend Dynamic_multi_parameter_filtration &operator*=(Dynamic_multi_parameter_filtration &f, const ValueRange &r)
  {
    for (Generator &g : f.generators_) {
      g *= r;
    }
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
  friend Dynamic_multi_parameter_filtration &operator*=(Dynamic_multi_parameter_filtration &f, const T &val)
  {
    for (Generator &g : f.generators_) {
      g *= val;
    }
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
   * @warning The operator accepts @ref Dynamic_multi_parameter_filtration with the same or different template
   * parameters as `ValueRange`. But if the number of generators is higher than 1, only the first generator will be
   * used for the operation.
   *
   * @tparam ValueRange Either a range of into `T` convertible elements with a begin() and end() method,
   * or @ref Dynamic_multi_parameter_filtration<U,...> with `U` convertible into `T`.
   * @param f First element of the division.
   * @param r Second element of the division.
   */
  template <class ValueRange, class = std::enable_if_t<RangeTraits<ValueRange>::has_begin> >
  friend Dynamic_multi_parameter_filtration operator/(Dynamic_multi_parameter_filtration f, const ValueRange &r)
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
   * @warning The operator accepts @ref Dynamic_multi_parameter_filtration with the same or different template
   * parameters as `ValueRange`. But if the number of generators is higher than 1, only the first generator will be
   * used for the operation.
   *
   * @tparam ValueRange Either a range of into `T` convertible elements with a begin() and end() method,
   * or @ref Dynamic_multi_parameter_filtration<U,...> with `U` convertible into `T`.
   * @param r First element of the division.
   * @param f Second element of the division.
   */
  template <class ValueRange,
            class = std::enable_if_t<RangeTraits<ValueRange>::has_begin &&
                                     !std::is_same_v<ValueRange, Dynamic_multi_parameter_filtration> > >
  friend Dynamic_multi_parameter_filtration operator/(const ValueRange &r, Dynamic_multi_parameter_filtration f)
  {
    for (Generator &g : f.generators_) {
      g = r / g;
    }
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
  friend Dynamic_multi_parameter_filtration operator/(Dynamic_multi_parameter_filtration f, const T &val)
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
  friend Dynamic_multi_parameter_filtration operator/(const T &val, Dynamic_multi_parameter_filtration f)
  {
    for (Generator &g : f.generators_) {
      g = val / g;
    }
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
   * @warning The operator accepts @ref Dynamic_multi_parameter_filtration with the same or different template
   * parameters as `ValueRange`. But if the number of generators is higher than 1, only the first generator will be
   * used for the operation.
   *
   * @tparam ValueRange Either a range of into `T` convertible elements with a begin() and end() method,
   * or @ref Dynamic_multi_parameter_filtration<U,...> with `U` convertible into `T`.
   * @param f First element of the division.
   * @param r Second element of the division.
   */
  template <class ValueRange, class = std::enable_if_t<RangeTraits<ValueRange>::has_begin> >
  friend Dynamic_multi_parameter_filtration &operator/=(Dynamic_multi_parameter_filtration &f, const ValueRange &r)
  {
    for (Generator &g : f.generators_) {
      g /= r;
    }
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
  friend Dynamic_multi_parameter_filtration &operator/=(Dynamic_multi_parameter_filtration &f, const T &val)
  {
    for (Generator &g : f.generators_) {
      g /= val;
    }
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

    generators_.resize(g);
  }

  /**
   * @brief If a generator is at +/- infinity or NaN, the underlying container can potentially (surely if just
   * constructed with default values) only contain one element, representing the value, instead of number of parameters
   * elements. This method forces the underlying container of the given generator to contain explicitly
   * an elements for each parameter, if it was not already the case.
   */
  void force_generator_size_to_number_of_parameters(size_type g)
  {
    if (g > num_generators()) return;
    generators_[g].force_size_to_number_of_parameters(number_of_parameters_);
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
    if constexpr (std::is_same_v<GeneratorRange, Generator>) {
      if (x.is_nan()) return false;
      bool xIsPlusInf = x.is_plus_inf();
      bool xIsMinusInf = x.is_minus_inf();
      if (!Co && xIsPlusInf) return false;
      if (Co && xIsMinusInf) return false;
      if (xIsPlusInf) {
        if (is_plus_inf()) return false;
        *this = inf(number_of_parameters_);
        return true;
      }
      if (xIsMinusInf) {
        if (is_minus_inf()) return false;
        *this = minus_inf(number_of_parameters_);
        return true;
      }
    }
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

    size_type end = num_generators();

    if (_generator_can_be_added(genStart, 0, end)) {
      generators_.resize(end);
      generators_.emplace_back(genStart, genEnd);
      if constexpr (Ensure1Criticality) {
        if (generators_.size() != 1)
          throw std::logic_error("Multiparameter filtration value is not 1-critical anymore.");
      }
      std::sort(generators_.begin(), generators_.end(), Is_strictly_smaller_lexicographically());
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

    if constexpr (std::is_same_v<GeneratorRange, Generator>) {
      generators_.push_back(x);
    } else {
      GUDHI_CHECK(x.size() == num_parameters(),
                  std::invalid_argument("Wrong range size. Should correspond to the number of parameters."));
      generators_.emplace_back(x.begin(), x.end());
    }
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
      size_type end = 0;

      for (std::size_t curr = 0; curr < generators_.size(); ++curr) {
        if (!generators_[curr].is_finite()) {
          if constexpr (Co) {
            if (generators_[curr].is_plus_inf()) {
              *this = inf(number_of_parameters_);
              return;
            }
          } else {
            if (generators_[curr].is_minus_inf()) {
              *this = minus_inf(number_of_parameters_);
              return;
            }
          }
          if (end == 0) ++end;  // if first element is +/-inf or nan, it should be kept at first
        } else {
          if (_generator_can_be_added(generators_[curr].begin(), 0, end)) {
            swap(generators_[end], generators_[curr]);
            ++end;
          }
        }
      }

      generators_.resize(end);
      std::sort(generators_.begin(), generators_.end(), Is_strictly_smaller_lexicographically());
    }
  }

  /**
   * @brief Removes all empty generators from the filtration value. If @p include_infinities is true, it also
   * removes the generators at infinity or minus infinity (or with NaN value).
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
      const Generator &g = generators_[0];
      if (g.size() == 0 || (include_infinities && !g.is_finite())) generators_.clear();
    } else {
      generators_.erase(std::remove_if(generators_.begin(),
                                       generators_.end(),
                                       [include_infinities](const Generator &a) {
                                         return a.size() == 0 || (include_infinities && !a.is_finite());
                                       }),
                        generators_.end());
      std::sort(generators_.begin(), generators_.end(), Is_strictly_smaller_lexicographically());
    }

    if (generators_.empty()) {
      generators_ = {Generator(1, Co ? T_inf : T_m_inf)};
    }
  }

  /**
   * @brief Sets each generator of the filtration value to the least common upper bound between it and the given value.
   *
   * More formally, it pushes the current generator to the cone \f$ \{ y \in \mathbb R^n : y \ge x \} \f$
   * originating in \f$ x \f$. The resulting value corresponds to the intersection of both
   * cones: \f$ \mathrm{this} = \min \{ y \in \mathbb R^n : y \ge this \} \cap \{ y \in \mathbb R^n : y \ge x \} \f$.
   *
   * @warning The operator accepts @ref Dynamic_multi_parameter_filtration with the same or different template
   * parameters as `GeneratorRange`. But if the number of generators is higher than 1, only the first generator will be
   * used for the operation.
   *
   * @tparam GeneratorRange Either a range of into `T` convertible elements with a begin(), end() and size() method,
   * or @ref Dynamic_multi_parameter_filtration<U,...> with `U` convertible into `T`.
   * @param x Range towards to push. Has to have as many elements than @ref num_parameters().
   * @param exclude_infinite_values If true, values at infinity or minus infinity are not affected.
   * @return true If the filtration value was actually modified.
   * @return false Otherwise.
   */
  template <class GeneratorRange = std::initializer_list<value_type>,
            class = std::enable_if_t<RangeTraits<GeneratorRange>::has_begin> >
  bool push_to_least_common_upper_bound(const GeneratorRange &x, bool exclude_infinite_values = false)
  {
    bool modified = false;

    for (Generator &g : generators_) {
      modified |= g.push_to_least_common_upper_bound(x, exclude_infinite_values);
    }

    if (modified && num_generators() > 1) simplify();

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
   * @warning The operator accepts @ref Dynamic_multi_parameter_filtration with the same or different template
   * parameters as `GeneratorRange`. But if the number of generators is higher than 1, only the first generator will be
   * used for the operation.
   *
   * @tparam GeneratorRange Either a range of into `T` convertible elements with a begin(), end() and size() method,
   * or @ref Dynamic_multi_parameter_filtration<U,...> with `U` convertible into `T`.
   * @param x Range towards to pull. Has to have as many elements than @ref num_parameters().
   * @param exclude_infinite_values If true, values at infinity or minus infinity are not affected.
   * @return true If the filtration value was actually modified.
   * @return false Otherwise.
   */
  template <class GeneratorRange = std::initializer_list<value_type>,
            class = std::enable_if_t<RangeTraits<GeneratorRange>::has_begin> >
  bool pull_to_greatest_common_lower_bound(const GeneratorRange &x, bool exclude_infinite_values = false)
  {
    bool modified = false;

    for (Generator &g : generators_) {
      modified |= g.pull_to_greatest_common_lower_bound(x, exclude_infinite_values);
    }

    if (modified && num_generators() > 1) simplify();

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

    for (Generator &g : generators_) {
      g.project_onto_grid(grid, coordinate);
    }

    if (!coordinate && num_generators() > 1) simplify();
  }

  // FONCTIONNALITIES

  /**
   * @brief Returns a generator with the minimal values of all parameters in any generator of the given filtration
   * value. That is, the greatest lower bound of all generators.
   */
  friend Dynamic_multi_parameter_filtration factorize_below(const Dynamic_multi_parameter_filtration &f)
  {
    if (f.num_generators() <= 1) return f;

    bool isMinusInf = true;
    bool nan = true;
    Underlying_container result(1, Generator(f.num_parameters(), T_inf));
    for (size_type p = 0; p < f.num_parameters(); ++p) {
      for (size_type g = 0; g < f.num_generators(); ++g) {
        T val = f(g, p);
        if (!_is_nan(val)) {
          nan = false;
          result[0][p] = val < result[0][p] ? val : result[0][p];
        }
      }
      if (nan)
        result[0][p] = std::numeric_limits<T>::quiet_NaN();
      else
        nan = true;
      if (result[0][p] != T_m_inf) isMinusInf = false;
    }

    if (isMinusInf) result = {Generator::minus_inf()};

    return Dynamic_multi_parameter_filtration(std::move(result), f.num_parameters());
  }

  /**
   * @brief Returns a generator with the maximal values of all parameters in any generator of the given filtration
   * value. That is, the least upper bound of all generators.
   */
  friend Dynamic_multi_parameter_filtration factorize_above(const Dynamic_multi_parameter_filtration &f)
  {
    if (f.num_generators() <= 1) return f;

    bool isPlusInf = true;
    bool nan = true;
    Underlying_container result(1, Generator(f.num_parameters(), T_m_inf));
    for (size_type p = 0; p < f.num_parameters(); ++p) {
      for (size_type g = 0; g < f.num_generators(); ++g) {
        T val = f(g, p);
        if (!_is_nan(val)) {
          nan = false;
          result[0][p] = val > result[0][p] ? val : result[0][p];
        }
      }
      if (nan)
        result[0][p] = std::numeric_limits<T>::quiet_NaN();
      else
        nan = true;
      if (result[0][p] != T_inf) isPlusInf = false;
    }

    if (isPlusInf) result = {Generator::inf()};

    return Dynamic_multi_parameter_filtration(std::move(result), f.num_parameters());
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
  friend U compute_linear_projection(const Dynamic_multi_parameter_filtration &f, const std::vector<U> &x)
  {
    if (f.num_generators() == 1) return compute_linear_projection(f.generators_[0], x);

    if constexpr (Co) {
      U projection = std::numeric_limits<U>::lowest();
      for (const Generator &g : f.generators_) {
        // Order in the max important to spread possible NaNs
        projection = std::max(compute_linear_projection(g, x), projection);
      }
      return projection;
    } else {
      U projection = std::numeric_limits<U>::max();
      for (const Generator &g : f.generators_) {
        // Order in the min important to spread possible NaNs
        projection = std::min(compute_linear_projection(g, x), projection);
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
  friend U compute_euclidean_distance_to(const Dynamic_multi_parameter_filtration &f,
                                         const Dynamic_multi_parameter_filtration &other)
  {
    GUDHI_CHECK(f.num_parameters() == other.num_parameters(),
                std::invalid_argument("We cannot compute the distance between two points of different dimensions."));

    // TODO: verify if this really makes a differences in the 1-critical case, otherwise just keep the general case
    if constexpr (Ensure1Criticality) {
      return compute_euclidean_distance_to(f.generators_[0], other.generators_[0]);
    } else {
      U res = std::numeric_limits<U>::max();
      for (const Generator &g1 : f.generators_) {
        for (const Generator &g2 : other.generators_) {
          // Order in the min important to spread possible NaNs
          res = std::min(compute_euclidean_distance_to(g1, g2), res);
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
  friend U compute_norm(const Dynamic_multi_parameter_filtration &f)
  {
    // Frobenius norm with matrix g x p based on Euclidean norm
    U out = 0;
    for (const Generator &g : f.generators_) {
      out += compute_squares(g);
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
  friend Dynamic_multi_parameter_filtration<OutValue, Co, Ensure1Criticality> compute_coordinates_in_grid(
      Dynamic_multi_parameter_filtration f,
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
  friend Dynamic_multi_parameter_filtration<U, Co, Ensure1Criticality> evaluate_coordinates_in_grid(
      const Dynamic_multi_parameter_filtration &f,
      const std::vector<std::vector<U> > &grid)
  {
    GUDHI_CHECK(grid.size() >= f.num_parameters(),
                std::invalid_argument(
                    "The size of the grid should correspond to the number of parameters in the filtration value."));

    std::vector<Multi_parameter_generator<U> > outVec(f.num_generators());

    size_type i = 0;
    for (const Generator &g : f.generators_) {
      outVec[i] = evaluate_coordinates_in_grid(g, grid);
      ++i;
    }

    Dynamic_multi_parameter_filtration<U, Co, Ensure1Criticality> out(std::move(outVec), f.num_parameters());
    if constexpr (!Ensure1Criticality)
      if (out.num_generators() > 1) out.simplify();
    return out;
  }

  // UTILITIES

  /**
   * @brief Outstream operator.
   */
  friend std::ostream &operator<<(std::ostream &stream, const Dynamic_multi_parameter_filtration &f)
  {
    const size_type num_gen = f.num_generators();
    const size_type num_param = f.num_parameters();

    stream << "( k = " << num_gen << " ) ( p = " << num_param << " ) [ ";
    for (size_type g = 0; g < num_gen; ++g) {
      stream << f.generators_[g];
      if (g < num_gen - 1) stream << "; ";
    }
    stream << " ]";

    return stream;
  }

  /**
   * @brief Instream operator.
   */
  friend std::istream &operator>>(std::istream &stream, Dynamic_multi_parameter_filtration &f)
  {
    size_type num_gen;
    size_type num_param;
    char delimiter;
    stream >> delimiter;  // (
    stream >> delimiter;  // k
    stream >> delimiter;  // =
    stream >> num_gen;
    if (!stream.good())
      throw std::invalid_argument("Invalid incoming stream format for Dynamic_multi_parameter_filtration (num_gen).");
    f.generators_.resize(num_gen);
    stream >> delimiter;  // )
    stream >> delimiter;  // (
    stream >> delimiter;  // p
    stream >> delimiter;  // =
    stream >> num_param;
    if (!stream.good())
      throw std::invalid_argument("Invalid incoming stream format for Dynamic_multi_parameter_filtration (num_param).");
    f.number_of_parameters_ = num_param;
    stream >> delimiter;  // )
    stream >> delimiter;  // [
    if (delimiter != '[')
      throw std::invalid_argument("Invalid incoming stream format for Dynamic_multi_parameter_filtration ([).");
    if (num_gen == 0) return stream;
    for (size_type i = 0; i < num_gen; ++i) {
      stream >> f.generators_[i];
      stream >> delimiter;  // ; or last ]
    }
    if (delimiter != ']')
      throw std::invalid_argument("Invalid incoming stream format for Dynamic_multi_parameter_filtration (]).");

    return stream;
  }

  /**
   * @brief Returns true if and only if the given filtration value is at plus infinity.
   */
  friend bool is_positive_infinity(const Dynamic_multi_parameter_filtration &f)
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
  friend bool unify_lifetimes(Dynamic_multi_parameter_filtration &f1, const Dynamic_multi_parameter_filtration &f2)
  {
    GUDHI_CHECK(f1.num_parameters() == f2.num_parameters() || !f1.is_finite() || !f2.is_finite(),
                "Cannot unify two filtration values with different number of parameters.");

    // TODO: verify if this really makes a differences in the 1-critical case, otherwise just keep the general case
    // if general case is kept: add (num_gen == 1) test to throw if unification is not 1-critical anymore.
    if constexpr (Ensure1Criticality) {
      // WARNING: costly check
      GUDHI_CHECK(f1 <= f2 || f2 <= f1,
                  "When 1-critical only, two non-comparable filtration values cannot be unified.");

      if constexpr (Co) {
        return f1.push_to_least_common_upper_bound(f2);
      } else {
        return f1.pull_to_greatest_common_lower_bound(f2);
      }
    } else {
      bool modified = false;
      for (const Generator &g : f2.generators_) {
        modified |= f1.add_generator(g);
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
  friend bool intersect_lifetimes(Dynamic_multi_parameter_filtration &f1, const Dynamic_multi_parameter_filtration &f2)
  {
    if (f1.is_nan() || f2.is_nan()) return false;

    if constexpr (Co) {
      if (f1.is_plus_inf()) {
        if (f2.is_plus_inf()) return false;
        f1 = f2;
        return true;
      }
      if (f1.is_minus_inf()) {
        return false;
      }
    } else {
      if (f1.is_minus_inf()) {
        if (f2.is_minus_inf()) return false;
        f1 = f2;
        return true;
      }
      if (f1.is_plus_inf()) {
        return false;
      }
    }

    GUDHI_CHECK(f1.num_parameters() == f2.num_parameters(),
                "Cannot intersect two filtration values with different number of parameters.");

    if constexpr (Ensure1Criticality) {
      if constexpr (Co) {
        return f1.pull_to_greatest_common_lower_bound(f2);
      } else {
        return f1.push_to_least_common_upper_bound(f2);
      }
    } else {
      const size_type num_param = f1.num_parameters();
      Dynamic_multi_parameter_filtration res = Co ? minus_inf(num_param) : inf(num_param);
      std::vector<T> newGen(num_param);
      // TODO: see if the order can be used to avoid g1 * g2 add_generator and
      // perhaps even to replace add_generator by add_guaranteed_generator
      for (size_type g1 = 0; g1 < f1.num_generators(); ++g1) {
        for (size_type g2 = 0; g2 < f2.num_generators(); ++g2) {
          GUDHI_CHECK(f1.generators_[g1].size() == num_param, "Filtration value f1 is not minimal.");
          GUDHI_CHECK(f2.generators_[g2].size() == num_param, "Filtration value f2 is not minimal.");
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
  friend char *serialize_value_to_char_buffer(const Dynamic_multi_parameter_filtration &value, char *start)
  {
    const std::size_t nberOfGenerators = value.generators_.size();
    const std::size_t type_size = sizeof(std::size_t);
    memcpy(start, &value.number_of_parameters_, type_size);
    memcpy(start + type_size, &nberOfGenerators, type_size);
    char *curr = start + type_size + type_size;
    for (const Generator &g : value) {
      curr = serialize_value_to_char_buffer(g, curr);
    }
    return curr;
  }

  /**
   * @brief Deserialize the value from a buffer at given pointer and stores it in given value.
   *
   * @param value Value to fill with the deserialized filtration value.
   * @param start Pointer to the start of the space in the buffer where the serialization is stored.
   * @return End position of the serialization in the buffer.
   */
  friend const char *deserialize_value_from_char_buffer(Dynamic_multi_parameter_filtration &value, const char *start)
  {
    const std::size_t type_size = sizeof(std::size_t);
    std::size_t nberOfGenerators;
    memcpy(&value.number_of_parameters_, start, type_size);
    memcpy(&nberOfGenerators, start + type_size, type_size);
    value.generators_.resize(nberOfGenerators);
    const char *curr = start + type_size + type_size;
    for (Generator &g : value) {
      curr = deserialize_value_from_char_buffer(g, curr);
    }
    return curr;
  }

  /**
   * @brief Returns the serialization size of the given filtration value.
   */
  friend std::size_t get_serialization_size_of(const Dynamic_multi_parameter_filtration &value)
  {
    std::size_t genSizes = sizeof(std::size_t) * 2;
    for (const Generator &g : value) {
      genSizes += get_serialization_size_of(g);
    }
    return genSizes;
  }

  /**
   * @brief Plus infinity value of an entry of the filtration value.
   */
  constexpr static const T T_inf = Generator::T_inf;

  /**
   * @brief Minus infinity value of an entry of the filtration value.
   */
  constexpr static const T T_m_inf = Generator::T_m_inf;

 private:
  size_type number_of_parameters_;  /**< Number of parameters. */
  Underlying_container generators_; /**< Container of the filtration value elements. */

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
  static bool _strictly_contains(const Generator &a, const Generator &b)
  {
    if constexpr (Co)
      return a > b;
    else {
      return a < b;
    }
  }

  /**
   * @brief Verifies if @p b is contained in the positive cone originating in `a`.
   */
  static bool _contains(const Generator &a, const Generator &b)
  {
    if constexpr (Co)
      return a >= b;
    else {
      return a <= b;
    }
  }

  /**
   * @brief Verifies if the first element of @p b strictly dominates the first element of `a`.
   */
  static bool _first_strictly_dominates(const Generator &a, const Generator &b)
  {
    if constexpr (Co) {
      return a.size() != 0 && b.size() != 0 && a[0] < b[0];
    } else {
      return a.size() != 0 && b.size() != 0 && a[0] > b[0];
    }
  }

  /**
   * @brief Verifies if the first element of @p b dominates the first element of `a`.
   */
  static bool _first_dominates(const Generator &a, const Generator &b)
  {
    if constexpr (Co) {
      return a.size() != 0 && b.size() != 0 && a[0] <= b[0];
    } else {
      return a.size() != 0 && b.size() != 0 && a[0] >= b[0];
    }
  }

  enum class Rel : std::uint8_t { EQUAL, DOMINATES, IS_DOMINATED, NONE };

  template <class Iterator>
  static Rel _get_domination_relation(const Generator &a, Iterator itB)
  {
    if (a.is_nan()) return Rel::NONE;

    bool equal = true;
    bool allGreater = true;
    bool allSmaller = true;
    bool allNaNA = true;
    bool allNaNB = true;
    for (unsigned int i = 0; i < a.size(); ++i) {
      if (a[i] < *itB) {
        if (!allSmaller) return Rel::NONE;
        equal = false;
        allGreater = false;
      } else if (a[i] > *itB) {
        if (!allGreater) return Rel::NONE;
        equal = false;
        allSmaller = false;
      }
      if (!_is_nan(a[i])) allNaNA = false;
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
   * by `generators_[curr]`. If x is dominated by or is equal to `generators_[curr]`, it cannot be added. If it
   * dominates `generators_[curr]`, it has to replace `generators_[curr]`. If there is no relation between both,
   * `generators_[curr]` has no influence on the addition of x.
   *
   * Assumes between 'curr' and 'end' everything is simplified:
   * no nan values and if there is an inf/-inf, then 'end - curr == 1'.
   */
  template <class Iterator>
  bool _generator_can_be_added(Iterator x, size_type curr, size_type &end)
  {
    // assumes that everything between curr and end is simplified
    // so, only generators_[curr] can be at inf or -inf.
    if constexpr (Co) {
      if (generators_[curr].is_plus_inf()) {
        return false;
      }
      if (generators_[curr].is_minus_inf()) {
        end = curr;
        return true;
      }
    } else {
      if (generators_[curr].is_minus_inf()) {
        return false;
      }
      if (generators_[curr].is_plus_inf()) {
        end = curr;
        return true;
      }
    }

    while (curr != end) {
      Rel res = _get_domination_relation(generators_[curr], x);
      if (res == Rel::IS_DOMINATED || res == Rel::EQUAL) return false;  // x dominates or is equal
      if (res == Rel::DOMINATES) {                                      // x is dominated
        --end;
        swap(generators_[curr], generators_[end]);
      } else {  // no relation
        ++curr;
      }
    }
    return true;
  }

  struct Is_strictly_smaller_lexicographically {
    // assumes both generators have the same length if not infinite/nan.
    bool operator()(const Generator &g1, const Generator &g2)
    {
      // orders such that -inf < 'finite values'  < inf < NaN.

      if (g1.is_nan() || g2.is_nan()) return !g1.is_nan();
      if (g1.is_plus_inf()) return false;
      if (g2.is_plus_inf()) return true;
      if (g2.is_minus_inf()) return false;
      if (g1.is_minus_inf()) return true;

      // g1 and g2 have to be finite and of the same size
      for (std::size_t i = 0; i < g1.size(); ++i) {
        if (g1[i] != g2[i]) return g1[i] < g2[i];
      }
      return false;
    }
  };
};

}  // namespace Gudhi::multi_filtration

namespace std {

template <typename T, bool Co, bool Ensure1Criticality>
class numeric_limits<Gudhi::multi_filtration::Dynamic_multi_parameter_filtration<T, Co, Ensure1Criticality> >
{
 public:
  using Filtration_value = Gudhi::multi_filtration::Dynamic_multi_parameter_filtration<T, Co, Ensure1Criticality>;

  static constexpr bool has_infinity = true;
  static constexpr bool has_quiet_NaN = true;

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

  static constexpr Filtration_value quiet_NaN(std::size_t p = 1) noexcept { return Filtration_value::nan(p); };
};

}  // namespace std

#endif  // MF_DYNAMIC_MULTI_PARAMETER_FILTRATION_H_
