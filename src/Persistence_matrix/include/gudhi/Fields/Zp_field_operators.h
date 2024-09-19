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
 * @file Zp_field_operators.h
 * @author Hannah Schreiber
 * @brief Contains the @ref Gudhi::persistence_fields::Zp_field_operators class.
 */

#ifndef MATRIX_FIELD_ZP_OPERATOR_H_
#define MATRIX_FIELD_ZP_OPERATOR_H_

#include <stdexcept>
#include <utility>
#include <vector>
#include <limits.h>

namespace Gudhi {
/// Field namespace
namespace persistence_fields {

/**
 * @class Zp_field_operators Zp_field_operators.h gudhi/Fields/Zp_field_operators.h
 * @ingroup persistence_fields
 *
 * @brief Class defining operators for the \f$ \mathbb{F}_p \f$ field for any prime number \f$ p \f$.
 *
 * @tparam Unsigned_integer_type A native unsigned integer type: unsigned int, long unsigned int, etc.
 * Will be used as the field element type.
 */
template <typename Unsigned_integer_type = unsigned int,
          class = std::enable_if_t<std::is_unsigned_v<Unsigned_integer_type> > >
class Zp_field_operators
{
 public:
  using Element = Unsigned_integer_type; /**< Type for the elements in the field. */
  using Characteristic = Element;   /**< Type for the field characteristic. */
  template <class T>
  using isSignedInteger = std::enable_if_t<std::is_signed_v<T> >;

  /**
   * @brief Default constructor. If a non-zero characteristic is given, initializes the field with it.
   * The characteristic can later be changed again or initialized with @ref set_characteristic.
   *
   * @param characteristic Prime number corresponding to the desired characteristic of the field.
   */
  Zp_field_operators(Characteristic characteristic = 0) : characteristic_(0) {
    if (characteristic != 0) set_characteristic(characteristic);
  }
  /**
   * @brief Copy constructor.
   * 
   * @param toCopy Operators to copy.
   */
  Zp_field_operators(const Zp_field_operators& toCopy)
      : characteristic_(toCopy.characteristic_), inverse_(toCopy.inverse_) {}
  /**
   * @brief Move constructor.
   * 
   * @param toMove Operators to move.
   */
  Zp_field_operators(Zp_field_operators&& toMove) noexcept
      : characteristic_(std::exchange(toMove.characteristic_, 0)), inverse_(std::move(toMove.inverse_)) {}

  /**
   * @brief Sets the characteristic of the field.
   * 
   * @param characteristic Prime number corresponding to the desired characteristic of the field.
   */
  void set_characteristic(Characteristic characteristic) {
    if (characteristic <= 1)
      throw std::invalid_argument("Characteristic must be strictly positive and a prime number.");

    inverse_.resize(characteristic);
    inverse_[0] = 0;
    for (unsigned int i = 1; i < characteristic; ++i) {
      unsigned int inv = 1;
      unsigned int mult = inv * i;
      while ((mult % characteristic) != 1) {
        ++inv;
        if (mult == characteristic) throw std::invalid_argument("Characteristic must be a prime number.");
        mult = inv * i;
      }
      inverse_[i] = inv;
    }

    characteristic_ = characteristic;
  }
  /**
   * @brief Returns the current characteristic.
   * 
   * @return The value of the current characteristic.
   */
  const Characteristic& get_characteristic() const { return characteristic_; }

  /**
   * @brief Returns the value of an integer in the field.
   * That is the positive value of the integer modulo the current characteristic.
   * 
   * @param e Unsigned integer to return the value from.
   * @return @p e modulo the current characteristic, such that the result is positive.
   */
  Element get_value(Element e) const { return e < characteristic_ ? e : e % characteristic_; }
  /**
   * @brief Returns the value of an integer in the field.
   * That is the positive value of the integer modulo the current characteristic.
   * 
   * @tparam Signed_integer_type A native signed integer type: int, long int, etc.
   * @param e Integer to return the value from.
   * @return @p e modulo the current characteristic, such that the result is positive.
   */
  template <typename Signed_integer_type, class = isSignedInteger<Signed_integer_type> >
  Element get_value(Signed_integer_type e) const {
    if (e < -static_cast<Signed_integer_type>(characteristic_)) e = e % characteristic_;
    if (e < 0) return e += characteristic_;
    return e < static_cast<Signed_integer_type>(characteristic_) ? e : e % characteristic_;
  }

  /**
   * @brief Returns the sum of two elements in the field.
   * 
   * @param e1 First element.
   * @param e2 Second element.
   * @return `(e1 + e2) % characteristic`, such that the result is positive.
   */
  Element add(Element e1, Element e2) const {
    return _add(get_value(e1), get_value(e2), characteristic_);
  }

  /**
   * @brief Stores in the first element the sum of two given elements in the field, that is
   * `(e1 + e2) % characteristic`, such that the result is positive.
   * 
   * @param e1 First element.
   * @param e2 Second element.
   */
  void add_inplace(Element& e1, Element e2) const {
    e1 = _add(get_value(e1), get_value(e2), characteristic_);
  }

  /**
   * @brief Returns the subtraction in the field of the first element by the second element.
   * 
   * @param e1 First element.
   * @param e2 Second element.
   * @return `(e1 - e2) % characteristic`, such that the result is positive.
   */
  Element subtract(Element e1, Element e2) const {
    return _subtract(get_value(e1), get_value(e2), characteristic_);
  }

  /**
   * @brief Stores in the first element the subtraction in the field of the first element by the second element,
   * that is `(e1 - e2) % 2`, such that the result is positive.
   * 
   * @param e1 First element.
   * @param e2 Second element.
   */
  void subtract_inplace_front(Element& e1, Element e2) const {
    e1 = _subtract(get_value(e1), get_value(e2), characteristic_);
  }
  /**
   * @brief Stores in the second element the subtraction in the field of the first element by the second element,
   * that is `(e1 - e2) % 2`, such that the result is positive.
   * 
   * @param e1 First element.
   * @param e2 Second element.
   */
  void subtract_inplace_back(Element e1, Element& e2) const {
    e2 = _subtract(get_value(e1), get_value(e2), characteristic_);
  }

  /**
   * @brief Returns the multiplication of two elements in the field.
   * 
   * @param e1 First element.
   * @param e2 Second element.
   * @return `(e1 * e2) % characteristic`, such that the result is positive.
   */
  Element multiply(Element e1, Element e2) const {
    return _multiply(get_value(e1), get_value(e2), characteristic_);
  }

  /**
   * @brief Stores in the first element the multiplication of two given elements in the field,
   * that is `(e1 * e2) % characteristic`, such that the result is positive.
   * 
   * @param e1 First element.
   * @param e2 Second element.
   */
  void multiply_inplace(Element& e1, Element e2) const {
    e1 = _multiply(get_value(e1), get_value(e2), characteristic_);
  }

  /**
   * @brief Multiplies the first element with the second one and adds the third one. Returns the result in the field.
   *
   * @warning Not overflow safe.
   * 
   * @param e First element.
   * @param m Second element.
   * @param a Third element.
   * @return `(e * m + a) % characteristic`, such that the result is positive.
   */
  Element multiply_and_add(Element e, Element m, Element a) const { return get_value(e * m + a); }

  /**
   * @brief Multiplies the first element with the second one and adds the third one, that is
   * `(e * m + a) % characteristic`, such that the result is positive. Stores the result in the first element.
   *
   * @warning Not overflow safe.
   * 
   * @param e First element.
   * @param m Second element.
   * @param a Third element.
   */
  void multiply_and_add_inplace_front(Element& e, Element m, Element a) const {
    e = get_value(e * m + a);
  }
  /**
   * @brief Multiplies the first element with the second one and adds the third one, that is
   * `(e * m + a) % characteristic`, such that the result is positive. Stores the result in the third element.
   *
   * @warning Not overflow safe.
   * 
   * @param e First element.
   * @param m Second element.
   * @param a Third element.
   */
  void multiply_and_add_inplace_back(Element e, Element m, Element& a) const {
    a = get_value(e * m + a);
  }

  /**
   * @brief Adds the first element to the second one and multiplies the third one with it.
   * Returns the result in the field.
   *
   * @warning Not overflow safe.
   * 
   * @param e First element.
   * @param a Second element.
   * @param m Third element.
   * @return `((e + a) * m) % characteristic`, such that the result is positive.
   */
  Element add_and_multiply(Element e, Element a, Element m) const { return get_value((e + a) * m); }

  /**
   * @brief Adds the first element to the second one and multiplies the third one with it, that is
   * `((e + a) * m) % characteristic`, such that the result is positive. Stores the result in the first element.
   *
   * @warning Not overflow safe.
   * 
   * @param e First element.
   * @param a Second element.
   * @param m Third element.
   */
  void add_and_multiply_inplace_front(Element& e, Element a, Element m) const {
    e = get_value((e + a) * m);
  }
  /**
   * @brief Adds the first element to the second one and multiplies the third one with it, that is
   * `((e + a) * m) % characteristic`, such that the result is positive. Stores the result in the third element.
   *
   * @warning Not overflow safe.
   * 
   * @param e First element.
   * @param a Second element.
   * @param m Third element.
   */
  void add_and_multiply_inplace_back(Element e, Element a, Element& m) const {
    m = get_value((e + a) * m);
  }

  /**
   * @brief Returns true if the two given elements are equal in the field, false otherwise.
   * 
   * @param e1 First element to compare.
   * @param e2 Second element to compare.
   * @return true If `e1 % characteristic == e2 % characteristic`.
   * @return false Otherwise.
   */
  bool are_equal(Element e1, Element e2) const { return get_value(e1) == get_value(e2); }

  /**
   * @brief Returns the inverse of the given element in the field.
   * 
   * @param e Element to get the inverse from.
   * @return Inverse in the current field of `e % characteristic`.
   */
  Element get_inverse(Element e) const { return inverse_[get_value(e)]; }
  /**
   * @brief For interface purposes with multi-fields. Returns the inverse together with the second argument.
   * 
   * @param e Element to get the inverse from.
   * @param productOfCharacteristics Some value.
   * @return Pair whose first element is the inverse of @p e and the second element is @p productOfCharacteristics.
   */
  std::pair<Element, Characteristic> get_partial_inverse(Element e,
                                                                   Characteristic productOfCharacteristics) const {
    return {get_inverse(e), productOfCharacteristics};
  }

  /**
   * @brief Returns the additive identity of the field.
   * 
   * @return 0.
   */
  static constexpr Element get_additive_identity() { return 0; }
  /**
   * @brief Returns the multiplicative identity of the field.
   * 
   * @return 1.
   */
  static constexpr Element get_multiplicative_identity() { return 1; }
  /**
   * @brief For interface purposes with multi-fields. Returns the multiplicative identity of the field.
   * 
   * @param productOfCharacteristics Some value.
   * @return 1.
   */
  static constexpr Element get_partial_multiplicative_identity(
      [[maybe_unused]] Characteristic productOfCharacteristics) {
    return 1;
  }

  // static constexpr bool handles_only_z2() { return false; }

  /**
   * @brief Assign operator.
   */
  Zp_field_operators& operator=(Zp_field_operators other) {
    std::swap(characteristic_, other.characteristic_);
    inverse_.swap(other.inverse_);
    return *this;
  }
  /**
   * @brief Swap operator.
   */
  friend void swap(Zp_field_operators& f1, Zp_field_operators& f2) {
    std::swap(f1.characteristic_, f2.characteristic_);
    f1.inverse_.swap(f2.inverse_);
  }

 private:
  Characteristic characteristic_;  /**< Current characteristic of the field. */
  std::vector<Element> inverse_;   /**< All inverse elements. */

  static Element _add(Element e1, Element e2, Characteristic characteristic) {
    if (UINT_MAX - e1 < e2) {
      // automatic unsigned integer overflow behaviour will make it work
      e1 += e2;
      e1 -= characteristic;
      return e1;
    }

    e1 += e2;
    if (e1 >= characteristic) e1 -= characteristic;

    return e1;
  }
  static Element _subtract(Element e1, Element e2, Characteristic characteristic) {
    if (e1 < e2) {
      e1 += characteristic;
    }
    e1 -= e2;

    return e1;
  }
  static Element _multiply(Element e1, Element e2, Characteristic characteristic) {
    unsigned int a = e1;
    e1 = 0;
    unsigned int temp_b;

    while (a != 0) {
      if (a & 1) {
        if (e2 >= characteristic - e1) e1 -= characteristic;
        e1 += e2;
      }
      a >>= 1;

      temp_b = e2;
      if (e2 >= characteristic - e2) temp_b -= characteristic;
      e2 += temp_b;
    }

    return e1;
  }
};

}  // namespace persistence_fields
}  // namespace Gudhi

#endif  // MATRIX_FIELD_ZP_OPERATOR_H_
