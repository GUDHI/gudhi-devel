/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2024 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file FieldOperators.h
 * @author Hannah Schreiber
 * @brief Contains the concept for the matrix field operators.
 */

namespace Gudhi {
namespace persistence_matrix {

/** 
 * @ingroup persistence_matrix
 *
 * @brief Concept of the field operator classes needed for the class @ref Matrix.
 *
 * Implementations of this concept are @ref Gudhi::persistence_fields::Zp_field_operators,
 * @ref Gudhi::persistence_fields::Z2_field_operators,
 * @ref Gudhi::persistence_fields::Multi_field_operators and
 * @ref Gudhi::persistence_fields::Multi_field_operators_with_small_characteristics.
 */
class FieldOperators 
{
 public:
  using element_type = unspecified;         /**< Type for the elements in the field. */
  using characteristic_type = unspecified;  /**< Type for the field characteristic. */

  /**
   * @brief Default constructor. If a non-zero characteristic is given, initializes the field with it.
   * The characteristic can later be changed again or initialized with @ref set_characteristic.
   *
   * @param characteristic Prime number corresponding to the desired characteristic of the field.
   */
  FieldOperators(characteristic_type characteristic = 0);

  /**
   * @brief Sets the characteristic of the field. Can eventually be omitted if the characteristic of the class
   * is fixed.
   * 
   * @param characteristic Prime number corresponding to the desired characteristic of the field.
   */
  void set_characteristic(characteristic_type characteristic);
  /**
   * @brief Returns the current characteristic.
   * 
   * @return The value of the current characteristic.
   */
  characteristic_type get_characteristic() const;

  /**
   * @brief Returns the value of an integer in the field.
   * That is the positive value of the integer modulo the current characteristic.
   * 
   * @tparam Integer_type A native integer type: int, unsigned int, long int, bool, etc.
   * @param e Integer to return the value from.
   * @return @p e modulo the current characteristic, such that the result is positive.
   */
  template <typename Integer_type>
  element_type get_value(Integer_type e) const;

  /**
   * @brief Returns the sum of two elements in the field.
   * 
   * @param e1 First element.
   * @param e2 Second element.
   * @return `(e1 + e2) % characteristic`, such that the result is positive.
   */
  element_type add(element_type e1, element_type e2) const;

  // /**
  //  * @brief Returns the substraction in the field of the first element by the second element.
  //  * 
  //  * @param e1 First element.
  //  * @param e2 Second element.
  //  * @return `(e1 - e2) % characteristic`, such that the result is positive.
  //  */
  // element_type substract(element_type e1, element_type e2) const;

  /**
   * @brief Returns the multiplication of two elements in the field.
   * 
   * @param e1 First element.
   * @param e2 Second element.
   * @return `(e1 * e2) % characteristic`, such that the result is positive.
   */
  element_type multiply(element_type e1, element_type e2) const;

  /**
   * @brief Multiplies the first element with the second one and adds the third one. Returns the result in the field.
   * 
   * @param e First element.
   * @param m Second element.
   * @param a Third element.
   * @return `(e * m + a) % characteristic`, such that the result is positive.
   */
  element_type multiply_and_add(element_type e, element_type m, element_type a) const;

  // /**
  //  * @brief Adds the first element to the second one and multiplies the third one with it.
  //  * Returns the result in the field.
  //  * 
  //  * @param e First element.
  //  * @param a Second element.
  //  * @param m Third element.
  //  * @return `((e + a) * m) % characteristic`, such that the result is positive.
  //  */
  // element_type add_and_multiply(element_type e, element_type a, element_type m) const;

  // /**
  //  * @brief Returns true if the two given elements are equal in the field, false otherwise.
  //  * 
  //  * @param e1 First element to compare.
  //  * @param e2 Second element to compare.
  //  * @return true If `e1 % characteristic == e2 % characteristic`.
  //  * @return false Otherwise.
  //  */
  // bool are_equal(element_type e1, element_type e2) const;

  /**
   * @brief Returns the inverse of the given element in the field.
   * 
   * @param e Element to get the inverse from.
   * @return Inverse in the current field of `e % characteristic`.
   */
  element_type get_inverse(element_type e) const;
  // /**
  //  * @brief In the case the field is a multi-field, returns the inverse of the given element in the fields
  //  * corresponding to the given sub-product of the product of all characteristics in the multi-field.
  //  * See @cite boissonnat:hal-00922572 for more details.
  //  * If the field is a usual field, simply returns the inverse in the field.
  //  * 
  //  * @param e Element to get the inverse from.
  //  * @param productOfCharacteristics Product of the different characteristics to take into account in the
  //  * multi-field.
  //  * @return If a multi-field: pair of the inverse of @p e and the characteristic the inverse is coming from.
  //  * If a normal field: pair of the inverse of @p e and @p productOfCharacteristics.
  //  */
  // std::pair<element_type, characteristic_type> get_partial_inverse(
  //   element_type e, characteristic_type productOfCharacteristics) const;

  /**
   * @brief Returns the additive identity of the field.
   * 
   * @return The additive identity of the field.
   */
  static constexpr element_type get_additive_identity();
  /**
   * @brief Returns the multiplicative identity of the field.
   * 
   * @return The multiplicative identity of the field.
   */
  static constexpr element_type get_multiplicative_identity();
  // /**
  //  * @brief For multi-fields, returns the partial multiplicative identity of the field from the given product.
  //  * See @cite boissonnat:hal-00922572 for more details.
  //  * Otherwise, simply returns the multiplicative identity of the field.
  //  *  
  //  * @param productOfCharacteristics Product of the different characteristics to take into account in the multi-field.
  //  * @return The partial multiplicative identity of the field
  //  */
  // static constexpr element_type get_partial_multiplicative_identity(
  //     [[maybe_unused]] characteristic_type productOfCharacteristics);

  /**
   * @brief Assign operator.
   */
  FieldOperators& operator=(FieldOperators other);
  /**
   * @brief Swap operator.
   */
  friend void swap(FieldOperators& f1, FieldOperators& f2);
};

}  // namespace persistence_matrix
}  // namespace Gudhi
