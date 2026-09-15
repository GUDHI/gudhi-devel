/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2026 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file StoragePolicy.h
 * @author Hannah Schreiber
 * @brief Contains the concept for the multiparameter filtration value class.
 */

/// Gudhi namespace.
namespace Gudhi {
/// Multi-filtration namespace.
namespace multi_filtration {

/**
 * @ingroup multi_filtration
 *
 * @brief Concept for the first template parameter of @ref Multi_parameter_filtration_value.
 * The internal order of generators and parameters has to be fixed, whatever order it is,
 * i.e. if two instances of @ref StoragePolicy are equal, the generators and parameters are ordered the same for both.
 */
class StoragePolicy {
 public:
  /**
   * @brief Type of an element of the filtration value at a parameter (in a generator).
   */
  using value_type = unspecified;
  /**
   * @brief Type of the underlying element container.
   */
  using Underlying_container = unspecified;
  /**
   * @brief Type of the container size.
   */
  using size_type = unspecified;
  /**
   * @brief Reference type for an element in the filtration value.
   */
  using reference = unspecified;
  /**
   * @brief Const reference type for an element in the filtration value.
   * Does not really have to be a r-value (e.g., can just be the same than `value_type`, that is a copy).
   */
  using const_reference = unspecified;
  /**
   * @brief Iterator type for a generator of the filtration value. Has to be at least LegacyForwardIterator.
   */
  using iterator = unspecified;
  /**
   * @brief Const iterator type for a generator of the filtration value. Has to be at least LegacyForwardIterator.
   */
  using const_iterator = unspecified;
  /**
   * @brief Returns the type of @ref StoragePolicy if @ref value_type is equal to the template parameter `U`.
   */
  template <typename U>
  using As_type = unspecified;

  /**
   * @brief Value considered +infinity for an element in the filtration value.
   * @note `-T_inf` does **not** have to be equal to `T_m_inf`.
   *
   * @tparam U Type of the element. Default value: @ref value_type.
   */
  template <typename U = value_type>
  constexpr static const U T_inf;

  /**
   * @brief Value considered -infinity for an element in the filtration value.
   * @note `-T_m_inf` does **not** have to be equal to `T_inf`.
   *
   * @tparam U Type of the element. Default value: @ref value_type.
   */
  template <typename U = value_type>
  constexpr static const U T_m_inf;

  /**
   * @brief True if and only if exactly one of the parameter values is just implicitly "stored" and is (therefore)
   * fixed/not updatable. See for example @ref Gudhi::multi_filtration::Degree_bifiltration.
   * For the comparaison operators `>`, `>=`, `<=` and `<` of @ref Multi_parameter_filtration_value to properly work
   * if needed, additional constrains are: the implicit values are integers (even if @ref value_type is not), the
   * generators are internally always ordered by increasing implicit values and no two distinct generators can have
   * the same value for this parameter.
   */
  constexpr static const bool has_an_implicit_axis;
  /**
   * @brief True if and only if @ref sort sorts lexicographically.
   */
  constexpr static const bool has_lexicographical_storage;
  /**
   * @brief True if and only if the storage strategy allows to simplify the set of generators such that it becomes
   * minimal.
   */
  constexpr static const bool has_minimal_set_representation;
  /**
   * @brief True if and only if @ref inf does not throw with the given template parameter.
   */
  template <bool Co>
  constexpr static const bool has_infinity;
  /**
   * @brief True if and only if @ref minus_inf does not throw with the given template parameter.
   */
  template <bool Co>
  constexpr static const bool has_minus_infinity;
  /**
   * @brief True if and only if @ref nan does not throw.
   */
  constexpr static const bool has_quiet_NaN;

  /**
   * @brief Constructs an empty container.
   */
  StoragePolicy();

  /**
   * @brief Constructs a single generator with given number of parameters and all values at given value.
   */
  StoragePolicy(size_type numberOfParameters, T value);

  /**
   * @brief Constructs a single generator from the given range.
   *
   * @tparam Iterator LegacyForwardIterator.
   * @param itBegin Begin iterator of the range.
   * @param itEnd End iterator of the range.
   */
  template <class Iterator>
  StoragePolicy(Iterator itBegin, Iterator itEnd);

  /**
   * @brief Builds filtration value with given number of parameters and values from the given range. Lets \f$ p \f$
   * be the number of parameters. The \f$ p \f$ first elements of the range have to correspond to the first generator,
   * the \f$ p \f$ next elements to the second generator and so on... So the length of the range has to be a multiple
   * of \f$ p \f$ and the number of generators will be \f$ length / p \f$. The range is represented by two iterators.
   *
   * @tparam Iterator LegacyForwardIterator.
   * @param itBegin Begin iterator of the range.
   * @param itEnd End iterator of the range.
   * @param numberOfParameters Number of parameters.
   */
  template <class Iterator, class = std::enable_if_t<!std::is_arithmetic_v<Iterator> > >
  StoragePolicy(Iterator itBegin, Iterator itEnd, size_type numberOfParameters);

  /**
   * @brief Builds filtration value with given number of parameters and values copied from the given
   * @ref Underlying_container container.
   */
  StoragePolicy(const Underlying_container &generators, size_type numberOfParameters);

  /**
   * @brief Builds filtration value with given number of parameters and values moved from the given
   * @ref Underlying_container container.
   */
  StoragePolicy(Underlying_container &&generators, size_type numberOfParameters);

  /**
   * @brief Copy constructor.
   */
  StoragePolicy(const StoragePolicy &other);

  /**
   * @brief Move constructor.
   */
  StoragePolicy(StoragePolicy &&other) noexcept;

  /**
   * @brief Assign operator.
   */
  StoragePolicy &operator=(const StoragePolicy &other);

  /**
   * @brief Move assign operator.
   */
  StoragePolicy &operator=(StoragePolicy &&other) noexcept;

  /**
   * @brief Swap operator.
   */
  friend void swap(StoragePolicy &f1, StoragePolicy &f2) noexcept;

  /**
   * @brief Returns a copy of this with given number of parameters and generators. If the given number of parameters
   * \f$ p \f$ is smaller than the current one, only the \f$ p \f$ first parameters are copied and if \f$ p \f$ is
   * greater than the current one, then additional parameters are added at the end using the given default value.
   * Same for the given number of generators.
   */
  template <typename U = value_type>
  As_type<U> copy(size_type numberOfParameters, size_type numberOfGenerators, U defaultValue) const;

  /**
   * @brief Only necessary if @ref has_an_implicit_axis is true. Returns the parameter index which is implicit.
   */
  static constexpr size_type implicit_axis() noexcept;

  /**
   * @brief Only necessary if @ref has_an_implicit_axis is true. Returns the generator index for which the implicit
   * parameter has the given implicit value. If the generator is not stored or the given value is invalid, has to
   * return @ref null_value with the output type as template parameter. If the given value is +/-infinity (and not
   * invalid), returns @ref T_inf, resp. @ref T_m_inf with the output type as template parameter.
   *
   * @param val Implicit value of the fixed parameter at the desired generator.
   */
  auto get_generator_of_implicit_value(value_type val) const;

  /**
   * @brief Only necessary if @ref has_an_implicit_axis is true. Returns the failing value of
   * @ref get_generator_of_implicit_value.
   *
   * @tparam U
   * @return constexpr U
   */
  template <typename U>
  static constexpr U null_value() noexcept;

  /**
   * @brief Only necessary if @ref has_an_implicit_axis is true. Returns true if and only if
   * `(*this)(g, implicit_axis()) == other(g, implicit_axis())` for every `g`.
   * This is meant to be a quick test.
   */
  bool has_same_implicit_values(const StoragePolicy& other) const;

  /**
   * @brief Returns reference to underlying container.
   */
  Underlying_container &get_underlying_container();

  /**
   * @brief Returns const reference to underlying container.
   */
  const Underlying_container &get_underlying_container() const;

  /**
   * @brief Returns a reference to element of parameter `p` of generator `g`.
   * If @ref has_an_implicit_axis is true and `p` is equal to @ref implicit_axis, the returned reference
   * can be pointing to a shared object or similar. The only thing to guarantee in that case is that the value
   * directly after the call is right. Changing the value of that object by the user does not have to have any
   * impact on the values of the stored generators.
   */
  reference operator()(size_type g, size_type p);

  /**
   * @brief Returns the value of parameter `p` of generator `g`.
   */
  const_reference operator()(size_type g, size_type p) const;

  /**
   * @brief Returns begin iterator to generator `g`.
   */
  iterator begin(size_type g);

  /**
   * @brief Returns end iterator to generator `g`.
   */
  iterator end(size_type g);

  /**
   * @brief Returns begin const iterator to generator `g`.
   */
  const_iterator begin(size_type g);

  /**
   * @brief Returns end const iterator to generator `g`.
   */
  const_iterator end(size_type g) const;

  /**
   * @brief Reserves space for the given number of generators in the underlying container.
   */
  void reserve(size_type number_of_generators);

  /**
   * @brief Returns the number of parameters in the filtration value.
   */
  size_type num_parameters() const;

  /**
   * @brief Returns the number of generators in the filtration value.
   */
  size_type num_generators() const;

  /**
   * @brief Returns the total number of values in the filtration value.
   */
  size_type num_entries() const;

  /**
   * @brief Returns a filtration value with given number of parameters for which @ref is_plus_inf() returns `true`
   * or an empty filtration value if `numberOfParameters` is 0.
   * Throws if @ref has_infinity is false.
   *
   * @tparam Co If `true`, reverses the poset order, i.e., the order \f$ \le \f$  in \f$ \mathbb R^n \f$ becomes
   * \f$ \ge \f$. That is, the positive cones representing a lifetime become all negative instead. Can be ignored
   * if there is no differences between the "two infinities".
   */
  template <bool Co>
  static StoragePolicy inf(size_type numberOfParameters);

  /**
   * @brief Returns a filtration value with given number of parameters for which @ref is_minus_inf() returns `true`
   * or an empty filtration value if `numberOfParameters` is 0.
   * Throws if @ref has_minus_infinity is false.
   *
   * @tparam Co If `true`, reverses the poset order, i.e., the order \f$ \le \f$  in \f$ \mathbb R^n \f$ becomes
   * \f$ \ge \f$. That is, the positive cones representing a lifetime become all negative instead. Can be ignored
   * if there is no differences between the "two -infinities".
   */
  template <bool Co>
  static StoragePolicy minus_inf(size_type numberOfParameters);

  /**
   * @brief Returns a filtration value with given number of parameters for which @ref is_nan() returns `true`
   * or an empty filtration value if `numberOfParameters` is 0 (for which @ref is_nan() sometimes can also evaluate
   * to true).
   * Throws if @ref has_quiet_NaN is false.
   */
  static StoragePolicy nan(size_type numberOfParameters);

  /**
   * @brief Returns `true` if and only if the filtration value is considered as +infinity with the given template
   * parameter.
   */
  template <bool Co>
  [[nodiscard]] bool is_plus_inf() const;

  /**
   * @brief Returns `true` if and only if the filtration value is considered as -infinity with the given template
   * parameter.
   */
  template <bool Co>
  [[nodiscard]] bool is_minus_inf() const;

  /**
   * @brief Returns `true` if and only if the filtration value is considered as NaN.
   */
  [[nodiscard]] bool is_nan() const;

  /**
   * @brief Returns `true` if and only if the filtration value is not considered as +/-infinity (with given
   * template parameter) nor NaN.
   */
  template <bool Co>
  [[nodiscard]] bool is_finite() const;

  /**
   * @brief Sets the number of generators. If there were less generators before, new generators with given value are
   * constructed. If there were more generators before, the exceed of generators is destroyed (any generator with index
   * higher or equal to @p g to be more precise). If `g` is zero, empties completely the container.
   *
   * @param g New number of generators.
   * @param defaultValue Value for each parameter for newly added generators.
   */
  void set_num_generators(size_type g, value_type defaultValue);

  /**
   * @brief Only necessary if @ref has_minimal_set_representation is **false**. Emplaces the given generator in the
   * stored set either at the end of the set or at the right position (w.r.t. the internal order defined by
   * @ref sort "" or by a forced internal order) if this can be done trivially. In the latter case, no insertion
   * can happen, if there is a generator overshadowing it (but it does not have to be tested, i.e. if the insertion
   * happens, it does not mean that there is no generator overshadowing it).
   *
   * @tparam Co If `true`, reverses the poset order, i.e., the order \f$ \le \f$  in \f$ \mathbb R^n \f$ becomes
   * \f$ \ge \f$. That is, the positive cones representing a lifetime become all negative instead.
   * @tparam Iterator LegacyForwardIterator.
   * @param startIt Begin iterator of the new generator.
   * @param endIt End iterator of the new generator.
   * @param nullValue Default value for eventual additional generators which existence is forced by the insertion.
   * @return true If the insertion actually happened.
   * @return false Otherwise.
   */
  template <bool Co, class Iterator>
  bool emplace(Iterator startIt, Iterator endIt, value_type nullValue);

  /**
   * @brief Only necessary if @ref has_minimal_set_representation is **true**. Emplaces the given generator in the
   * stored set at the end without any verification.
   *
   * @tparam Iterator LegacyForwardIterator.
   */
  template <class Iterator>
  void emplace_back(Iterator startIt, Iterator endIt);

  /**
   * @brief Only necessary if @ref has_minimal_set_representation is **true**.
   * Swaps position of the two given generators.
   */
  void swap_generators(size_type g1, size_type g2);

  /**
   * @brief Sorts the set of generators if the order of the generators is not always fixed (in which case it does
   * nothing). The order has to be total, i.e., for a same set of generators, the new order has to be always the same.
   * If @ref has_lexicographical_storage is true, the sort has also to be lexicographical. The "words" considered for
   * the lexicographical order are all the generators concatenated together in order of generator index and then in
   * order of parameter index.
   */
  void sort();

  /**
   * @brief Serialize given value into the buffer at given pointer.
   *
   * @param value Value to serialize.
   * @param start Pointer to the start of the space in the buffer where to store the serialization.
   * @return End position of the serialization in the buffer.
   */
  friend char *serialize_value_to_char_buffer(const StoragePolicy &value, char *start);

  /**
   * @brief Deserialize the value from a buffer at given pointer and stores it in given value.
   *
   * @param value Value to fill with the deserialized filtration value.
   * @param start Pointer to the start of the space in the buffer where the serialization is stored.
   * @return End position of the serialization in the buffer.
   */
  friend const char *deserialize_value_from_char_buffer(StoragePolicy &value, const char *start);

  /**
   * @brief Returns the serialization size of the given filtration value.
   */
  friend std::size_t get_serialization_size_of(const StoragePolicy &value);
};

/**
 * @brief Only necessary if @ref Multi_parameter_filtration_value::as_type has to be used. Returns a copy of the given
 * @ref StoragePolicy but in the format of another class (or same with different template parameters) respecting the
 * concept. Has to exist in the `Gudhi::multi_filtration` namespace.
 * The possibilities of `OutStoragePolicy` only has to cover the users need of
 * @ref Multi_parameter_filtration_value::as_type.
 *
 * @tparam OutStoragePolicy Class respecting the @ref StoragePolicy concept.
 * @tparam Co The second template parameter of the desired @ref Multi_parameter_filtration_value returned by
 * @ref Multi_parameter_filtration_value::as_type. Default value: false.
 */
template <class OutStoragePolicy, bool Co = false>
inline OutStoragePolicy as_type(const StoragePolicy &f);

}  // namespace multi_filtration
}  // namespace Gudhi
