/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - 2024/08 Hannah Schreiber: Update of the concept after several additions to the Simplex tree.
 *      - YYYY/MM Author: Description of the modification
 */

/** \brief Value type for a filtration function on a cell complex.
 *
 * Needs to implement `std::numeric_limits<FiltrationValue>::has_infinity`,
 * `std::numeric_limits<FiltrationValue>::infinity()` and `std::numeric_limits<FiltrationValue>::max()`.
 * But when `std::numeric_limits<FiltrationValue>::has_infinity` returns `true`,
 * `std::numeric_limits<FiltrationValue>::max()` can simply throw when called, as well as,
 * `std::numeric_limits<FiltrationValue>::infinity()` if `std::numeric_limits<FiltrationValue>::has_infinity`
 * returns `false`.
 *
 * A <EM>filtration</EM> of a cell complex (see FilteredComplex) is
 * a function \f$f:\mathbf{K} \rightarrow \mathbb{R}\f$ satisfying \f$f(\tau)\leq
 * f(\sigma)\f$ whenever \f$\tau \subseteq \sigma\f$. Ordering the simplices
 * by increasing filtration values (breaking ties so as a simplex appears after
 * its subsimplices of same filtration value) provides an indexing scheme
 * (see IndexingTag).
 */
struct FiltrationValue {
  /**
   * @brief Has to construct the default value of FiltrationValue.
   */
  FiltrationValue();
  /**
   * @brief Only necessary if `FiltrationValue(0)` does not correspond to an "empty" filtration value, or,
   * if there is an ambiguity at compile time between several constructors if
   * @ref Gudhi::simplex_tree::empty_filtration_value_t is passed.
   * Should construct a representative of what is an "empty" object. E.g., 0 for a numerical value or {} for a vector.
   * Two empty objects should evaluate to true if compared for equality.
   */
  explicit FiltrationValue(Gudhi::simplex_tree::empty_filtration_value_t e);

  // only for default ordering of filtration_vect_ in initialize_filtration and for prune_above_filtration
  /**
   * @brief Strictly smaller operator. If the filtration values are totally ordered, should be a StrictWeakOrdering.
   */
  friend bool operator<(const FiltrationValue& f1, const FiltrationValue& f2);
  /**
   * @brief Equality operator
   */
  friend bool operator==(const FiltrationValue& f1, const FiltrationValue& f2);

  /**
   * @brief Given two filtration values at which a simplex exists, computes the minimal union of births generating
   * a lifetime including those two values. The result is stored in the first parameter.
   * The overload for arithmetic types like `double` or `int` is already implemented as the minimum of the
   * two given values and can also be used for non native arithmetic types like `CGAL::Gmpq` as long as it has an
   * `operator<`. The overload is available with @ref Gudhi::unify_lifetimes "".
   *
   * For a k-critical filtration, FiltrationValue should be able to store an union of values (corresponding to the
   * different births of a same simplex) and this method adds the values of @p f2 in @p f1 and removes the values
   * from @p f1 which are comparable and greater than other values.
   * In the special case of 1-critical filtration, as the union should not contain more than one birth element,
   * this method is expected to throw if the two given elements in the filtration values are not comparable.
   * If they are comparable, the union is simply the minimum of both.
   *
   * @return True if and only if the values in @p f1 were actually modified.
   */
  friend bool unify_lifetimes(FiltrationValue& f1, const FiltrationValue& f2);

  /**
   * @brief Given two filtration values, stores in the first value the lowest common upper bound of the two values.
   * The overload for arithmetic types like `double` or `int` is already implemented as the maximum of the two
   * given values and can also be used for non native arithmetic types like `CGAL::Gmpq` as long as it has an
   * `operator<`. The overload is available with @ref Gudhi::intersect_lifetimes "".
   *
   * @return True if and only if the values in @p f1 were actually modified.
   */
  friend bool intersect_lifetimes(FiltrationValue& f1, const FiltrationValue& f2);

  /**
   * @private
   * @brief Only necessary when serializing the simplex tree. Serialize the given value and insert it at start position.
   * Overloads for native arithmetic types or other simple types are already implemented.
   * 
   * @param value The value to serialize.
   * @param start Start position where the value is serialized.
   * @return The new position in the array of char for the next serialization.
   */
  friend char* serialize_value_to_char_buffer(const FiltrationValue& value, char* start);

  /**
   * @private
   * @brief Only necessary when deserializing the simplex tree. Deserializes the object at the given start position in
   * an array of char, which was previously serialized by @ref serialize_value_to_char_buffer "".
   * Then, sets the value with it.
   * Overloads for native arithmetic types or other simple types are already implemented.
   * 
   * @param value The value where to deserialize based on its type.
   * @param start Start position where the value is serialized.
   * @return The new position in the array of char for the next deserialization.
   */
  friend const char* deserialize_value_from_char_buffer(FiltrationValue& value, const char* start);

  /**
   * @private
   * @brief Only necessary when serializing the simplex tree. Returns the serialization size of the given object.
   * Overloads for native arithmetic types or other simple types are already implemented.
   */
  friend std::size_t get_serialization_size_of(const FiltrationValue& value);
};
