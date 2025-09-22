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
  using Element = unspecified;         /**< Type for the elements in the field. */
  using Characteristic = unspecified;  /**< Type for the field characteristic. */

  /**
   * @brief Indicates the value of the characteristic when not initialized.
   */
  inline static const Characteristic nullCharacteristic;

  /**
   * @brief Default constructor. If a non-zero characteristic is given, initializes the field with it.
   * The characteristic can later be changed again or initialized with @ref set_characteristic.
   *
   * @param characteristic Prime number corresponding to the desired characteristic of the field.
   */
  FieldOperators(Characteristic characteristic = 0);

  /**
   * @brief Sets the characteristic of the field. Can eventually be omitted if the characteristic of the class
   * is fixed.
   * 
   * @param characteristic Prime number corresponding to the desired characteristic of the field.
   */
  void set_characteristic(const Characteristic& characteristic);
  /**
   * @brief Returns the current characteristic.
   * 
   * @return The value of the current characteristic.
   */
  const Characteristic& get_characteristic() const;

  /**
   * @brief Returns the value of an integer in the field.
   * That is the positive value of the integer modulo the current characteristic.
   * 
   * @tparam Integer_type A native integer type: int, unsigned int, long int, bool, etc.
   * @param e Integer to return the value from.
   * @return @p e modulo the current characteristic, such that the result is positive.
   */
  template <typename Integer_type>
  Element get_value(Integer_type e) const;

  // void get_value(Element& e) const;

  /**
   * @brief Stores in the first element the sum of two given elements in the field, that is
   * `(e1 + e2) % characteristic`, such that the result is positive.
   * 
   * @param e1 First element.
   * @param e2 Second element.
   */
  void add_inplace(Element& e1, const Element& e2) const;

  /**
   * @brief Stores in the first element the multiplication of two given elements in the field,
   * that is `(e1 * e2) % characteristic`, such that the result is positive.
   * 
   * @param e1 First element.
   * @param e2 Second element.
   */
  void multiply_inplace(Element& e1, const Element& e2) const;

  /**
   * @brief Multiplies the first element with the second one and adds the third one, that is
   * `(e * m + a) % characteristic`, such that the result is positive. Stores the result in the first element.
   * 
   * @param e First element.
   * @param m Second element.
   * @param a Third element.
   */
  void multiply_and_add_inplace_front(Element& e, const Element& m, const Element& a) const;
  /**
   * @brief Multiplies the first element with the second one and adds the third one, that is
   * `(e * m + a) % characteristic`, such that the result is positive. Stores the result in the third element.
   * 
   * @param e First element.
   * @param m Second element.
   * @param a Third element.
   */
  void multiply_and_add_inplace_back(const Element& e, const Element& m, Element& a) const;

  /**
   * @brief Returns the inverse of the given element in the field.
   * 
   * @param e Element to get the inverse from.
   * @return Inverse in the current field of `e % characteristic`.
   */
  Element get_inverse(const Element& e) const;

  /**
   * @brief Returns the additive identity of the field.
   * 
   * @return The additive identity of the field.
   */
  static const Element& get_additive_identity();
  /**
   * @brief Returns the multiplicative identity of the field.
   * 
   * @return The multiplicative identity of the field.
   */
  static const Element& get_multiplicative_identity();

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
